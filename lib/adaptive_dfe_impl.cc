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

#include <gnuradio/io_signature.h>
#include <volk/volk.h>
#include "adaptive_dfe_impl.h"

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
}

adaptive_dfe::sptr
adaptive_dfe::make(int sps, // samples per symbol
                   int nB,  // number of forward FIR taps
                   int nF,  // number of backward FIR taps
                   int nW,  // number of feedback taps
                   std::string python_module_name)
{
  return gnuradio::get_initial_sptr
    (new adaptive_dfe_impl(sps, nB, nF, nW, python_module_name));
}

/*
 * The private constructor
 */
adaptive_dfe_impl::adaptive_dfe_impl(int sps, // samples per symbol
                                     int nB,  // number of forward FIR taps
                                     int nF,  // number of backward FIR taps
                                     int nW,  // number of feedback taps
                                     std::string python_module_name)
  : gr::block("adaptive_dfe",
              gr::io_signature::make(1, 1, sizeof(gr_complex)),
              gr::io_signature::make(1, 1, sizeof(gr_complex)))
  , _sps(sps)
  , _nB(nB)
  , _nF(nF)
  , _nW(nW)
  , _mu(0.01)
  , _alpha(0.0005)
  , _py_module_name(python_module_name)
  , _physicalLayer()
  , _taps_samples(nullptr)
  , _taps_symbols(nullptr)
  , _hist_samples(nullptr)
  , _hist_symbols(nullptr)
  , _hist_sample_index(0)
  , _hist_symbol_index(0)
  , _sample_counter(0)
  , _constellations()
  , _constellation_index()
  , _symbols()
  , _scramble()
  , _descrambled_symbols()
  , _symbol_counter(0)
  , _sum_phase_diff(0)
  , _df(0)
  , _phase(0)
  , _b{0.338187046465954, -0.288839024460507}
  , _ud(0)
  , _state(WAIT_FOR_PREAMBLE)
{
}

/*
 * Our virtual destructor.
 */
adaptive_dfe_impl::~adaptive_dfe_impl()
{
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

  int nout = 0;
  int i = 0;
  for (; i<ninput_items[0] && nout < noutput_items; ++i) {
    assert(nout < noutput_items);
    _phase += _df;
    if (_phase > M_PI)
      _phase -= 2*M_PI;
    if (_phase < -M_PI)
      _phase += 2*M_PI;

    _hist_samples[_hist_sample_index] = _hist_samples[_hist_sample_index+_nB+_nF+1] = in[i] * std::exp(gr_complex(0,_phase));
    if (++_hist_sample_index == _nB+_nF+1)
      _hist_sample_index = 0;

    if (_state == WAIT_FOR_PREAMBLE) {
      std::vector<tag_t> v;
      get_tags_in_window(v, 0, i,i+1);
      float phase_est = 0;
      float  corr_est = 0;
      uint64_t offset = 0;
      for (int j=0; j<v.size(); ++j) {
        std::cout << "tag " << v[j].key << " " << v[j].offset-nitems_read(0) << std::endl;
        if (v[j].key == pmt::mp("phase_est")) {
          phase_est = pmt::to_double(v[j].value);
          std::cout << "phase_est " << v[j].offset <<" " << nitems_read(0) << " " << phase_est << std::endl;
        }
        if (v[j].key == pmt::mp("corr_est")) {
          corr_est = pmt::to_double(v[j].value);
          std::cout << "corr_est " << v[j].offset << " " << nitems_read(0) << " "
                    << pmt::is_number(v[j].value) << " "
                    << pmt::is_integer(v[j].value) << " "
                    << pmt::is_real(v[j].value) << " "
                    << pmt::to_double(v[j].value)
                    << std::endl;
          if (corr_est > 130e3) {
            offset = v[j].offset - nitems_read(0);
            break;
          }
        }
      }
      if (corr_est > 130e3) {
        _state = DO_FILTER;
        _sample_counter = 0;
        _symbol_counter = 0;
        // _symbols.clear();
        // _scramble.clear();
        _descrambled_symbols.clear();
        // _hist_sample_index = 0;
        _hist_symbol_index = 0;
        std::fill_n(_hist_symbols, 2*_nW, gr_complex(0));
        std::fill_n(_taps_samples, _nB+_nF+1, gr_complex(0));
        std::fill_n(_taps_symbols, _nW, gr_complex(0));
        //_phase = -phase_est;
        _taps_samples[_nB+1] = std::exp(gr_complex(0, -phase_est));
        _taps_symbols[0] = 1;
        GILLock lock;
        try {
          update_frame_information(_physicalLayer.attr("get_frame")());
        } catch (boost::python::error_already_set const&) {
          PyErr_Print();
        }
      }
    }
    if (_state == DO_FILTER) {
      gr_complex dot_samples = 0;
      // volk_32fc_x2_dot_prod_32fc(&dot_samples,
      //                            &_hist_samples.front()+_hist_sample_index,
      //                            &_taps_samples.front(),
      //                            _taps_samples.size());
      gr_complex filter_output = dot_samples;
      // if (_sample_counter < 80*5)
      //   std::cout << "SAMPLE " << _sample_counter << " " << dot_samples << std::endl;
      if ((_sample_counter%_sps) == 0) {
        if (_symbol_counter == _symbols.size()) {
          _symbol_counter = 0;
          GILLock lock;
          try {
            boost::python::numpy::ndarray s = boost::python::numpy::from_data(&_descrambled_symbols.front(),
                                                                              boost::python::numpy::dtype::get_builtin<gr_complex>(),
                                                                              boost::python::make_tuple(_descrambled_symbols.size()),
                                                                              boost::python::make_tuple(sizeof(gr_complex)),
                                                                              boost::python::object());
            update_doppler_information(_physicalLayer.attr("get_doppler")(s));
            update_frame_information(_physicalLayer.attr("get_frame")());
          } catch (boost::python::error_already_set const&) {
            PyErr_Print();
          }
        }
        gr_complex known_symbol = _symbols[_symbol_counter];
        bool is_known = true;
        for (int k=0; k<1; ++k) {
          filter_output = 0;
#if 1
          volk_32fc_x2_dot_prod_32fc(&filter_output,
                                     _hist_samples+_hist_sample_index,
                                     _taps_samples,
                                     _nB+_nF+1);
#else
          for (int l=0; l<_nB+_nF+1; ++l) {
            assert(_hist_sample_index+l < 2*(_nB+_nF+1));
            filter_output += _hist_samples[_hist_sample_index+l]*_taps_samples[l];
          }
#endif
          gr_complex dot_symbols=0;
          for (int l=0; l<_nW; ++l) {
            assert(_hist_symbol_index+l < 2*_nW);
            dot_symbols += _hist_symbols[_hist_symbol_index+l]*_taps_symbols[l];
          }
          filter_output += dot_symbols;
          if (std::abs(known_symbol) < 1e-5) { // not known
            is_known = false;
            gr_complex descrambled_filter_output = std::conj(_scramble[_symbol_counter]) * filter_output;
            gr::digital::constellation_sptr constell = _constellations[_constellation_index];
            unsigned int jc = constell->decision_maker(&descrambled_filter_output);
            constell->map_to_points(jc, &descrambled_filter_output);
            known_symbol = _scramble[_symbol_counter] * descrambled_filter_output;
          }
          gr_complex err =  filter_output - known_symbol;
          if (_symbol_counter >= 0) {
            for (int j=0; j<_nB+_nF+1; ++j) {
              _taps_samples[j] -= _mu*err*std::conj(_hist_samples[_hist_sample_index+j]);
            }
            for (int j=0; j<_nW; ++j) {
              assert(_hist_symbol_index+j < 2*_nW);
              _taps_symbols[j] -= _mu*err*std::conj(_hist_symbols[_hist_symbol_index+j]) + _alpha*_taps_symbols[j];
            }
          }
          // std::cout << "filter: " << _symbol_counter << " " << _sample_counter << " " << filter_output << " " << known_symbol << " " << std::abs(err) << std::endl;
        }
        if (is_known) {
          _taps_symbols[_hist_symbol_index] = _taps_symbols[_hist_symbol_index + _nW] = known_symbol;
          if (++_hist_symbol_index == _nW)
            _hist_symbol_index = 0;
        }
        _descrambled_symbols[_symbol_counter] = filter_output*std::conj(_scramble[_symbol_counter]);
        out[nout++] = filter_output;
        ++_symbol_counter;
      }
      _sample_counter += 1;
    }
  }

  consume(0, i);

  // Tell runtime system how many output items we produced.
  return nout;
}

bool adaptive_dfe_impl::start()
{
  // make sure python is ready for threading
  if( Py_IsInitialized() ){
    GILLock lock;
    if(PyEval_ThreadsInitialized() != 1 ){
      PyEval_InitThreads();
    }
    boost::python::numpy::initialize();
  } else {
    throw std::runtime_error("dont use es_pyhandler without python!");
  }
  _taps_samples = (gr_complex*)(volk_malloc(  (_nB+_nF+1)*sizeof(gr_complex), volk_get_alignment()));
  _taps_symbols = (gr_complex*)(volk_malloc(          _nW*sizeof(gr_complex), volk_get_alignment()));
  _hist_samples = (gr_complex*)(volk_malloc(2*(_nB+_nF+1)*sizeof(gr_complex), volk_get_alignment()));
  _hist_symbols = (gr_complex*)(volk_malloc(        2*_nW*sizeof(gr_complex), volk_get_alignment()));
  _taps_samples[_nB+1] = 1;
  _taps_symbols[0]     = 1;

  std::cout << "adaptive_dfe_impl::start()" << std::endl;
  GILLock lock;
  try {
    boost::python::object module        = boost::python::import(boost::python::str("digitalhf.physical_layer." + _py_module_name));
    boost::python::object PhysicalLayer = module.attr("PhysicalLayer");
    _physicalLayer = PhysicalLayer();
    update_constellations(_physicalLayer.attr("get_constellations")());
  } catch (boost::python::error_already_set const&) {
    PyErr_Print();
    return false;
  }
  return true;
}
bool adaptive_dfe_impl::stop()
{
  std::cout << "adaptive_dfe_impl::stop()" << std::endl;
  GILLock lock;
  _physicalLayer = boost::python::object();
  volk_free(_taps_samples);
  volk_free(_taps_symbols);
  volk_free(_hist_samples);
  volk_free(_hist_symbols);
  return true;
}

void adaptive_dfe_impl::update_constellations(boost::python::object obj)
{
  int const n = boost::python::extract<int>(obj.attr("__len__")());
  _constellations.resize(n);
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
  }
}
void adaptive_dfe_impl::update_frame_information(boost::python::object obj)
{
  int const n = boost::python::extract<int>(obj.attr("__len__")());
  assert(n==2);
  boost::python::numpy::ndarray array = boost::python::numpy::array(obj[0]);
  char const* data = array.get_data();
  int  const     m = array.shape(0);
  _symbols.resize(m);
  _scramble.resize(m);
  _descrambled_symbols.resize(m);
  for (int i=0; i<m; ++i) {
    std::memcpy(&_symbols[i],  data+16*i,   sizeof(gr_complex));
    std::memcpy(&_scramble[i], data+16*i+8, sizeof(gr_complex));
    // std::cout << "get_frame " << i << " " << _symbols[i] << " " << _scramble[i] << std::endl;
  }
  _constellation_index = boost::python::extract<int>(obj[1]);
}
void adaptive_dfe_impl::update_doppler_information(boost::python::object obj)
{
  int const n = boost::python::extract<int>(obj.attr("__len__")());
  assert(n==2);
  double const do_continue = boost::python::extract<bool>(obj[0]);
  double const doppler = boost::python::extract<float>(obj[1]);

  float delta_f = doppler/_sps;
  if (_df == 0) { // init
    _ud = _df = -delta_f;
  } else {
    const float ud_old = _ud;
    _ud  = -delta_f;
    _df +=_b[0]*_ud + _b[1]*ud_old;
  }
  std::cout << "PLL: " << _df << " " << delta_f << std::endl;
  _sum_phase_diff = 0;
}

} /* namespace digitalhf */
} /* namespace gr */
