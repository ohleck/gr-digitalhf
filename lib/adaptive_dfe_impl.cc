/* -*- c++ -*- */
/*
 * Copyright 2018 hcab14@mail.com.
 *
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <boost/format.hpp>

#include <gnuradio/expj.h>
#include <gnuradio/io_signature.h>
#include <gnuradio/logger.h>

#include <volk/volk.h>

#include "adaptive_dfe_impl.h"

#define VOLK_SAFE_DELETE(x) \
  volk_free(x);             \
  x = nullptr

namespace gr {
namespace digitalhf {

namespace {
class GILLock {
  PyGILState_STATE _state;
public:
  GILLock()
    :_state(PyGILState_Ensure()) {}
  ~GILLock() {
    PyGILState_Release(_state);
  }
} ;

boost::python::numpy::ndarray
complex_vector_to_ndarray(std::vector<gr_complex> const& v) {
  return  boost::python::numpy::from_data
    (&v.front(),
     boost::python::numpy::dtype::get_builtin<gr_complex>(),
     boost::python::make_tuple(v.size()),
     boost::python::make_tuple(sizeof(gr_complex)),
     boost::python::object());
}
} // anonymous namespace

adaptive_dfe::sptr
adaptive_dfe::make(int sps, // samples per symbol
                   int nB,  // number of forward FIR taps
                   int nF,  // number of backward FIR taps
                   int nW,  // number of feedback taps
                   float mu,
                   float alpha,
                   std::string python_module_name)
{
  return gnuradio::get_initial_sptr
    (new adaptive_dfe_impl(sps, nB, nF, nW, mu, alpha, python_module_name));
}

/*
 * The private constructor
 */
adaptive_dfe_impl::adaptive_dfe_impl(int sps, // samples per symbol
                                     int nB,  // number of forward FIR taps
                                     int nF,  // number of backward FIR taps
                                     int nW,  // number of feedback taps
                                     float mu,
                                     float alpha,
                                     std::string python_module_name)
  : gr::block("adaptive_dfe",
              gr::io_signature::make(1, 1, sizeof(gr_complex)),
              gr::io_signature::make(1, 1, sizeof(gr_complex)))
  , _sps(sps)
  , _nB(nB*sps)
  , _nF(nF*sps)
  , _nW(nW)
  , _mu(mu)
  , _alpha(alpha)
  , _use_symbol_taps(true)
  , _py_module_name(python_module_name)
  , _physicalLayer()
  , _taps_samples(nullptr)
  , _taps_symbols(nullptr)
  , _hist_samples(nullptr)
  , _hist_symbols(nullptr)
  , _hist_sample_index(0)
  , _hist_symbol_index(0)
  , _ignore_filter_updates(0)
  , _saved_samples()
  , _sample_counter(0)
  , _constellations()
  , _npwr()
  , _npwr_counter()
  , _npwr_max_time_constant(10)
  , _constellation_index()
  , _samples()
  , _symbols()
  , _scramble()
  , _descrambled_symbols()
  , _symbol_counter(0)
  , _need_samples(false)
  , _save_soft_decisions(false)
  , _vec_soft_decisions()
  , _msg_port_name(pmt::mp("soft_dec"))
  , _msg_metadata(pmt::make_dict())
  , _df(0)
  , _phase(0)
  , _b{0.338187046465954, -0.288839024460507}
  , _ud(0)
  , _state(WAIT_FOR_PREAMBLE)
{
  GR_LOG_DECLARE_LOGPTR(d_logger);
  GR_LOG_ASSIGN_LOGPTR(d_logger, "adaptive_dfe");
  message_port_register_out(_msg_port_name);
}

/*
 * Our virtual destructor.
 */
adaptive_dfe_impl::~adaptive_dfe_impl()
{
  _msg_port_name = pmt::PMT_NIL;
  _msg_metadata  = pmt::PMT_NIL;
  VOLK_SAFE_DELETE(_taps_samples);
  VOLK_SAFE_DELETE(_taps_symbols);
  VOLK_SAFE_DELETE(_hist_samples);
  VOLK_SAFE_DELETE(_hist_symbols);
}

void
adaptive_dfe_impl::forecast (int noutput_items, gr_vector_int &ninput_items_required)
{
  ninput_items_required[0] = _sps*noutput_items;
}

int
adaptive_dfe_impl::general_work(int noutput_items,
                                gr_vector_int &ninput_items,
                                gr_vector_const_void_star &input_items,
                                gr_vector_void_star &output_items)
{
  gr::thread::scoped_lock lock(d_setlock);
  gr_complex const* in = (gr_complex const *)input_items[0];
  gr_complex *out = (gr_complex *)output_items[0];

  int nout = 0; // counter for produced output items
  int i    = 0; // counter for consumed input items
  for (; i<ninput_items[0] && nout < noutput_items;) {
    assert(nout < noutput_items);
    switch (_state) {
      case WAIT_FOR_PREAMBLE: {
        insert_sample(in[i++]);
        uint64_t offset = 0;
        float phase_est = 0;
        if (get_correlation_tag(i, offset, phase_est)) {
          GR_LOG_DEBUG(d_logger, "next state > INITIAL_DOPPLER_ESTIMATE");
          _state = INITIAL_DOPPLER_ESTIMATE;
          _sample_counter = 0;
          _symbol_counter = 0;
          // _symbols.clear();
          // _scramble.clear();
          _descrambled_symbols.clear();
          // _hist_sample_index = 0;
          _hist_symbol_index = 0;
          _ignore_filter_updates = 0;
          _saved_samples.clear();
          std::fill_n(_hist_symbols, 2*_nW, gr_complex(0));
          std::fill_n(_taps_samples, _nB+_nF+1, gr_complex(0));
          std::fill_n(_taps_symbols, _nW, gr_complex(0));
          _use_symbol_taps = true;
          _samples.clear();
          _phase = -phase_est;
          _taps_samples[_nB+1] = 0.01;
          _taps_symbols[0] = 1;
          GILLock gil_lock;
          try {
            update_frame_information(_physicalLayer.attr("get_frame")());
          } catch (boost::python::error_already_set const&) {
            PyErr_Print();
          }
        }
        break;
      } // WAIT_FOR_PREAMBLE

      case INITIAL_DOPPLER_ESTIMATE: {
        _samples.push_back(in[i++]);
        // buffer samples and replay them later once the initial doppler estimate is there
        if (_samples.size() == _sps * _symbols.size()) {
          GILLock gil_lock;
          try {
            std::vector<gr_complex> const empty_vec;
            // initial doppler estimate
            if (!update_doppler_information(_physicalLayer.attr("get_doppler")
                                            (complex_vector_to_ndarray(empty_vec),
                                             complex_vector_to_ndarray(_samples)))) {
              GR_LOG_DEBUG(d_logger, "next state > WAIT_FOR_PREAMBLE");
              _state = WAIT_FOR_PREAMBLE;
              break;
            }
          } catch (boost::python::error_already_set const&) {
            PyErr_Print();
          }
          // (1) correct all samples in the circular buffer with the inital doppler estimate
          for (int j=_nB+1; j<_nB+_nF+1; ++j) {
            assert(_hist_sample_index+j < 2*(_nB+_nF+1));
            _hist_samples[_hist_sample_index+j] *= gr_expj(-_phase);
            update_local_oscillator();
          }
          // (2) insert all buffered samples and run the adaptive filter for them
          //      instead of pop_front() we first reverse _samples and then insert back() + pop_back()
          //      O(N) instead of O(N^2)
          std::reverse(_samples.begin(), _samples.end());
          while (!_samples.empty() && nout < noutput_items) {
            insert_sample(_samples.back());
            _sample_counter += 1;
            _samples.pop_back();
            if ((_sample_counter%_sps) == 0)
              out[nout++] = filter();
          }
          if (_samples.empty()) {
            GR_LOG_DEBUG(d_logger,"next state > DO_FILTER");
            _state = DO_FILTER;
            break;
          } else {
            GR_LOG_DEBUG(d_logger, "next state > INITIAL_DOPPLER_ESTIMATE_CONTINUE");
            _state = INITIAL_DOPPLER_ESTIMATE_CONTINUE;
            break;
          }
        }
      } // INITIAL_DOPPLER_ESTIMATE_CONTINUE

      case INITIAL_DOPPLER_ESTIMATE_CONTINUE: {
        GR_LOG_DEBUG(d_logger, "INITIAL_DOPPLER_ESTIMATE_CONTINUE");
        while (!_samples.empty() && nout < noutput_items) {
          insert_sample(_samples.back());
          _sample_counter += 1;
          _samples.pop_back();
          if ((_sample_counter%_sps) == 0)
            out[nout++] = filter();
        }
        if (_samples.empty()) {
          GR_LOG_DEBUG(d_logger, "next state > DO_FILTER");
          _state = DO_FILTER;
        } else {
          GR_LOG_DEBUG(d_logger, "next state > INITIAL_DOPPLER_ESTIMATE_CONTINUE");
          _state = INITIAL_DOPPLER_ESTIMATE_CONTINUE;
        }
        break;
      } // INITIAL_DOPPLER_ESTIMATE_CONTINUE

      case DO_FILTER: {
        if ((_sample_counter%_sps) == 0) {
          if (_symbol_counter == _symbols.size()) { // frame is ready
            _symbol_counter = 0;
            GILLock gil_lock;
            try {
              // update doppler estimate
              if (!update_doppler_information(_physicalLayer.attr("get_doppler")
                                              (complex_vector_to_ndarray(_descrambled_symbols),
                                               complex_vector_to_ndarray(_samples)))) {
                GR_LOG_DEBUG(d_logger, "next state > WAIT_FOR_PREAMBLE");
                _state = WAIT_FOR_PREAMBLE;
                break;
              }
              // publish soft decisions
              if (!_vec_soft_decisions.empty()) {
                std::cout << "soft_dec " << _vec_soft_decisions.size() << "\n";
                unsigned int const bits_per_symbol = _constellations[_constellation_index]->bits_per_symbol();
                _msg_metadata = pmt::dict_add(_msg_metadata, pmt::mp("bits_per_symbol"), pmt::from_long(bits_per_symbol));
                _msg_metadata = pmt::dict_add(_msg_metadata, pmt::mp("packet_len"), pmt::mp(_vec_soft_decisions.size()));
                message_port_pub(_msg_port_name,
                                 pmt::cons(_msg_metadata,
                                           pmt::init_f32vector(_vec_soft_decisions.size(), _vec_soft_decisions)));
                _vec_soft_decisions.clear();
              }
              _samples.clear();
              // get information about the following frame
              update_frame_information(_physicalLayer.attr("get_frame")());
            } catch (boost::python::error_already_set const&) {
              PyErr_Print();
            }
          } // frame is ready
          if (_ignore_filter_updates == 0) {
            out[nout++] = filter();
            if (_symbol_counter+1 == _symbols.size())
              recenter_filter_taps();
          } else {
            _ignore_filter_updates -= 1;
          }
        } // (_sample_counter%_sps) == 0
        if (_need_samples) {
          _samples.push_back(_hist_samples[_hist_sample_index+_nB+1]);
        }
        if (_saved_samples.empty()) {
          insert_sample(in[i++]);
        } else {
          insert_sample(_saved_samples.back());
          _saved_samples.pop_back();
        }
        _sample_counter += 1;
      } // DO_FILTER
    } // switch _state
  } // next input sample

  consume(0, i);

  // Tell runtime system how many output items we produced.
  return nout;
}

bool adaptive_dfe_impl::start()
{
  gr::thread::scoped_lock lock(d_setlock);
  // make sure python is ready for threading
  if( Py_IsInitialized() ){
    GILLock gil_lock;
    if(PyEval_ThreadsInitialized() != 1 ){
      PyEval_InitThreads();
    }
    boost::python::numpy::initialize();
  } else {
    throw std::runtime_error("dont use adaptive_dfe without python!");
  }
  _taps_samples = (gr_complex*)(volk_malloc(  (_nB+_nF+1)*sizeof(gr_complex), volk_get_alignment()));
  _taps_symbols = (gr_complex*)(volk_malloc(          _nW*sizeof(gr_complex), volk_get_alignment()));
  _hist_samples = (gr_complex*)(volk_malloc(2*(_nB+_nF+1)*sizeof(gr_complex), volk_get_alignment()));
  _hist_symbols = (gr_complex*)(volk_malloc(        2*_nW*sizeof(gr_complex), volk_get_alignment()));

  _samples.clear();
  std::fill_n(_hist_samples, 2*(_nB+_nF+1), gr_complex(0));
  std::fill_n(_hist_symbols,         2*_nW, gr_complex(0));
  std::fill_n(_taps_samples,   (_nB+_nF+1), gr_complex(0));
  std::fill_n(_taps_symbols,           _nW, gr_complex(0));

  _taps_samples[_nB+1] = 0.01;
  _taps_symbols[0]     = 1;

  GR_LOG_DEBUG(d_logger,str(boost::format("adaptive_dfe_impl::start() nB=%d nF=%d mu=%f alpha=%f")
                             % _nB % _nF % _mu % _alpha));
  GILLock gil_lock;
  try {
    boost::python::object module        = boost::python::import(boost::python::str("digitalhf.physical_layer." + _py_module_name));
    boost::python::object PhysicalLayer = module.attr("PhysicalLayer");
    _physicalLayer = PhysicalLayer(_sps);
    update_constellations(_physicalLayer.attr("get_constellations")());
  } catch (boost::python::error_already_set const&) {
    PyErr_Print();
    return false;
  }
  return true;
}
bool adaptive_dfe_impl::stop()
{
  gr::thread::scoped_lock lock(d_setlock);
  GR_LOG_DEBUG(d_logger, "adaptive_dfe_impl::stop()");
  GILLock gil_lock;
  _physicalLayer = boost::python::object();
  VOLK_SAFE_DELETE(_taps_samples);
  VOLK_SAFE_DELETE(_taps_symbols);
  VOLK_SAFE_DELETE(_hist_samples);
  VOLK_SAFE_DELETE(_hist_symbols);
  return true;
}

gr_complex adaptive_dfe_impl::filter() {
  gr_complex filter_output = 0;
  volk_32fc_x2_dot_prod_32fc(&filter_output,
                             _hist_samples+_hist_sample_index,
                             _taps_samples,
                             _nB+_nF+1);
  gr_complex dot_symbols=0;
  gr::digital::constellation_sptr constell = _constellations[_constellation_index];
  bool const update_taps = constell->bits_per_symbol() <= 3;
  if (constell->bits_per_symbol() > 3)
    _use_symbol_taps = false;
  if (_use_symbol_taps) {
    for (int l=0; l<_nW; ++l) {
      assert(_hist_symbol_index+l < 2*_nW);
      dot_symbols += _hist_symbols[_hist_symbol_index+l]*_taps_symbols[l];
    }
    filter_output += dot_symbols;
  }
  gr_complex known_symbol = _symbols[_symbol_counter];
  bool const is_known     = std::abs(known_symbol) > 1e-5;
  if (not is_known) { // not known
    gr_complex const descrambled_filter_output = std::conj(_scramble[_symbol_counter]) * filter_output;
    unsigned int jc = constell->decision_maker(&descrambled_filter_output);
    gr_complex descrambled_symbol = 0;
    constell->map_to_points(jc, &descrambled_symbol);

    if (_save_soft_decisions) {
      float const err = std::abs(descrambled_filter_output - descrambled_symbol);
      _npwr_counter[_constellation_index] += (_npwr_counter[_constellation_index] < _npwr_max_time_constant);
      float const alpha = 1.0f/_npwr_counter[_constellation_index];
      _npwr[_constellation_index] = (1-alpha)*_npwr[_constellation_index] + alpha*err;
      std::vector<float> const soft_dec = constell->calc_soft_dec(descrambled_filter_output, _npwr[_constellation_index]);
      std::copy(soft_dec.begin(), soft_dec.end(), std::back_inserter<std::vector<float> >(_vec_soft_decisions));
    }
    known_symbol = _scramble[_symbol_counter] * descrambled_symbol;
  }
  if (is_known || update_taps) {
    gr_complex const err =  filter_output - known_symbol;
    for (int j=0; j<_nB+_nF+1; ++j) {
      assert(_hist_sample_index+j < 2*(_nB+_nF+1));
      _taps_samples[j] -= _mu*err*std::conj(_hist_samples[_hist_sample_index+j]);
    }
    if (_use_symbol_taps) {
      for (int j=0; j<_nW; ++j) {
        assert(_hist_symbol_index+j < 2*_nW);
        _taps_symbols[j] -= _mu*err*std::conj(_hist_symbols[_hist_symbol_index+j]) + _alpha*_taps_symbols[j];
      }
      _hist_symbols[_hist_symbol_index] = _hist_symbols[_hist_symbol_index + _nW] = known_symbol;
      if (++_hist_symbol_index == _nW)
        _hist_symbol_index = 0;
    }
  }
  _descrambled_symbols[_symbol_counter] = filter_output*std::conj(_scramble[_symbol_counter]);
  return filter_output*std::conj(_scramble[_symbol_counter++]);
}

void adaptive_dfe_impl::recenter_filter_taps() {
  // get max(abs(taps))
  ssize_t const idx_max = std::distance(_taps_samples,
                                        std::max_element(_taps_samples+_nB+1-3*_sps, _taps_samples+_nB+1+3*_sps,
                                                         [](gr_complex a, gr_complex b) {
                                                           return std::norm(a) < std::norm(b);
                                                         }));
  GR_LOG_DEBUG(d_logger, str(boost::format("idx_max=%2d abs(tap_max)=%f") % idx_max % std::abs(_taps_samples[idx_max])));

  if (idx_max-_nB-1 >= 2*_sps && _saved_samples.empty() && _ignore_filter_updates==0) {
    // maximum is right of the center tap
    //   -> shift taps to the left left
    GR_LOG_DEBUG(d_logger, "shift left");
    std::copy(_taps_samples+2*_sps, _taps_samples+_nB+_nF+1, _taps_samples);
    std::fill_n(_taps_samples+_nB+_nF+1-2*_sps, 2*_sps, gr_complex(0));
    //      and omit the next two calls to filter in order to keep the alignment between samples and taps
    _ignore_filter_updates = 2;

  } else if (idx_max-_nB-1 <= -2*_sps && _saved_samples.empty() && _ignore_filter_updates==0) {
    // maximum is left of the center tap
    //   -> shift taps to the right
    GR_LOG_DEBUG(d_logger, "shift right");
    std::copy_backward(_taps_samples, _taps_samples+_nB+_nF+1-2*_sps,
                       _taps_samples+_nB+_nF+1);
    std::fill_n(_taps_samples, 2*_sps, gr_complex(0));
    //      save the last 2*_sps samples (will be reinserted)
    _saved_samples.resize(2*_sps);
    std::reverse_copy(_hist_samples+_hist_sample_index+(_nB+_nF+1)-2*_sps,
                      _hist_samples+_hist_sample_index+(_nB+_nF+1),
                      _saved_samples.begin());
    //      shift samples index
    _hist_sample_index += (_nB+_nF+1)-2*_sps;
    _hist_sample_index %= (_nB+_nF+1);
    //      set the 1st 2*_sps unknown old samples to zero
    for (int l=_hist_sample_index; l<_hist_sample_index+2*_sps; ++l) {
      int const k = (l+_nB+_nF+1)%(2*(_nB+_nF+1));
      _hist_samples[l] = _hist_samples[k] = gr_complex(0);
    }
  }
}
void adaptive_dfe_impl::set_mode(std::string mode) {
  gr::thread::scoped_lock lock(d_setlock);
  GR_LOG_DEBUG(d_logger, "adaptive_dfe_impl::set_mode "+ mode);
  GILLock gil_lock;
  try {
    _physicalLayer.attr("set_mode")(mode);
  } catch (boost::python::error_already_set const&) {
    PyErr_Print();
    return;
  }
}

void adaptive_dfe_impl::update_constellations(boost::python::object obj)
{
  int const n = boost::python::extract<int>(obj.attr("__len__")());
  _constellations.resize(n);
  _npwr.resize(n);
  _npwr_counter.resize(n);
  for (int i=0; i<n; ++i) {
    boost::python::numpy::ndarray const& array = boost::python::numpy::array(obj[i]);
    char const* data = array.get_data();
    int  const     m = array.shape(0);
    std::vector<gr_complex> constell(m);
    std::vector<int>   pre_diff_code(m);
    for (int j=0; j<m; ++j) {
      std::memcpy(&constell[j], data+9*j, sizeof(gr_complex));
      pre_diff_code[j] = (data+9*j)[8];
    }
    unsigned int const rotational_symmetry = 0;
    unsigned int const dimensionality = 1;
    _constellations[i] = gr::digital::constellation_calcdist::make(constell, pre_diff_code, rotational_symmetry, dimensionality);
    _npwr[i] = 0.0f;
    _npwr_counter[i] = 0;
  }
}
bool adaptive_dfe_impl::update_frame_information(boost::python::object obj)
{
  int const n = boost::python::extract<int>(obj.attr("__len__")());
  assert(n==4);
  boost::python::numpy::ndarray array = boost::python::numpy::array(obj[0]);
  char const* data = array.get_data();
  int  const     m = array.shape(0);
  _symbols.resize(m);
  _scramble.resize(m);
  _descrambled_symbols.resize(m);
  _samples.clear();
  for (int i=0; i<m; ++i) {
    std::memcpy(&_symbols[i],  data+16*i,   sizeof(gr_complex));
    std::memcpy(&_scramble[i], data+16*i+8, sizeof(gr_complex));
    // std::cout << "get_frame " << i << " " << _symbols[i] << " " << _scramble[i] << std::endl;
  }
  _constellation_index = boost::python::extract<int> (obj[1]);
  _need_samples        = boost::python::extract<bool>(obj[2]);
  _save_soft_decisions = boost::python::extract<bool>(obj[3]);
  return true;
}
bool adaptive_dfe_impl::update_doppler_information(boost::python::object obj)
{
  int const n = boost::python::extract<int>(obj.attr("__len__")());
  assert(n==2);
  bool  const do_continue = boost::python::extract<bool>(obj[0]);
  if (!do_continue) {
    _phase = 0;
    _df    = 0;
    std::fill_n(_hist_samples, 2*(_nB+_nF+1), gr_complex(0));
    _hist_sample_index = 0;
    _sample_counter    = 0;
    return false;
  }
  float const doppler    = boost::python::extract<float>(obj[1]);
  update_pll(doppler);
  return true;
}

void adaptive_dfe_impl::update_pll(float doppler) {
  if (doppler == 0)
    return;
  float const delta_f = doppler/_sps;
  if (_df == 0) { // init
    _ud = _df = delta_f;
  } else {
    float const ud_old = _ud;
    _ud  = delta_f;
    _df +=_b[0]*_ud + _b[1]*ud_old;
  }
  GR_LOG_DEBUG(d_logger, str(boost::format("PLL: df=%f delta_f=%f (rad/sample)") % _df % delta_f));
}
void adaptive_dfe_impl::insert_sample(gr_complex z) {
  // insert sample into the circular buffer
  _hist_samples[_hist_sample_index] = _hist_samples[_hist_sample_index+_nB+_nF+1] = z * gr_expj(-_phase);
  if (++_hist_sample_index == _nB+_nF+1)
    _hist_sample_index = 0;
  if (z != gr_complex(0))
    update_local_oscillator();
}
void adaptive_dfe_impl:: update_local_oscillator() {
  _phase += _df;
  if (_phase > M_PI)
    _phase -= 2*M_PI;
  if (_phase < -M_PI)
    _phase += 2*M_PI;
}
bool adaptive_dfe_impl::get_correlation_tag(uint64_t i, uint64_t& offset, float& phase_est) {
  std::vector<tag_t> v;
  get_tags_in_window(v, 0, i,i+1);
  for (int j=0; j<v.size(); ++j) {
    std::cout << "tag " << v[j].key << " " << v[j].offset-nitems_read(0) << " " << v[j].value << std::endl;
    if (v[j].key == pmt::mp("phase_est")) {
      phase_est = pmt::to_double(v[j].value);
      std::cout << "phase_est " << v[j].offset <<" " << nitems_read(0) << " " << phase_est << std::endl;
      offset = v[j].offset - nitems_read(0);
    }
    if (v[j].key == pmt::mp("corr_est")) {
      double const corr_est = pmt::to_double(v[j].value);
      if (v[j].offset - nitems_read(0) == offset)// && corr_est > 10e3)
        return true;
    }
  }
  return false;
}

} /* namespace digitalhf */
} /* namespace gr */
