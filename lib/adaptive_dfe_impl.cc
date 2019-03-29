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

#include <gnuradio/math.h>
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

adaptive_dfe::sptr
adaptive_dfe::make(int sps, // samples per symbol
                   int nB,  // number of forward FIR taps
                   int nF,  // number of backward FIR taps
                   int nW,  // number of feedback taps
                   float mu,
                   float alpha)
{
  return gnuradio::get_initial_sptr
    (new adaptive_dfe_impl(sps, nB, nF, nW, mu, alpha));
}

adaptive_dfe_impl::adaptive_dfe_impl(int sps, // samples per symbol
                                     int nB,  // number of forward FIR taps
                                     int nF,  // number of backward FIR taps
                                     int nW,  // number of feedback taps
                                     float mu,
                                     float alpha)
  : gr::block("adaptive_dfe",
              gr::io_signature::make(1, 1, sizeof(gr_complex)),
              gr::io_signature::make2(2, 2,
                                      sizeof(gr_complex),
                                      sizeof(gr_complex)*(sps*(nF+nB)+1)))
  , _sps(sps)
  , _nB(nB*sps)
  , _nF(nF*sps)
  , _nW(nW)
  , _nGuard(2*sps)
  , _mu(mu)
  , _alpha(alpha)
  , _use_symbol_taps(true)
  // , _py_module_name(python_module_name)
  // , _physicalLayer()
  , _taps_samples(nullptr)
  , _taps_symbols(nullptr)
  , _hist_symbols(nullptr)
  , _hist_symbol_index(0)
  , _constellations()
  , _npwr()
  , _npwr_max_time_constant(10)
  , _constellation_index()
  , _symbols()
  , _scramble()
  , _descrambled_symbols()
  , _symbol_counter(0)
  , _save_soft_decisions(false)
  , _vec_soft_decisions()
  , _msg_ports{{"soft_dec",   pmt::intern("soft_dec")},
               {"frame_info", pmt::intern("frame_info")}}
  , _msg_metadata(pmt::make_dict())
  , _state(WAIT_FOR_PREAMBLE)
  , _num_samples_since_filter_update(0)
  , _rotated_samples()
  , _rotator()
  , _control_loop(2*M_PI/100, 5e-3, -5e-3)
{
  GR_LOG_DECLARE_LOGPTR(d_logger);
  GR_LOG_ASSIGN_LOGPTR(d_logger, "adaptive_dfe");

  set_history(_nGuard+_nB+1);

  message_port_register_out(_msg_ports["soft_dec"]);

  pmt::pmt_t constellations_port = pmt::intern("constellations");
  message_port_register_in(constellations_port);
  set_msg_handler(constellations_port, boost::bind(&adaptive_dfe_impl::update_constellations, this, _1));

  pmt::pmt_t frame_info_port = _msg_ports["frame_info"];
  message_port_register_in(frame_info_port);
  message_port_register_out(frame_info_port);
  set_msg_handler(frame_info_port, boost::bind(&adaptive_dfe_impl::update_frame_info, this, _1));
}

adaptive_dfe_impl::~adaptive_dfe_impl()
{
  _msg_metadata = pmt::PMT_NIL;
  VOLK_SAFE_DELETE(_taps_samples);
  VOLK_SAFE_DELETE(_taps_symbols);
  VOLK_SAFE_DELETE(_hist_symbols);
}

void
adaptive_dfe_impl::forecast(int noutput_items, gr_vector_int &ninput_items_required)
{
  // [guard | nB | 1 | nF | guard ]
  ninput_items_required[0] = _sps*noutput_items + 2*_nGuard + _nB + _nF + 1;
}

int
adaptive_dfe_impl::general_work(int noutput_items,
                                gr_vector_int &ninput_items,
                                gr_vector_const_void_star &input_items,
                                gr_vector_void_star &output_items)
{
  gr::thread::scoped_lock lock(d_setlock);
  //GR_LOG_DEBUG(d_logger, str(boost::format("work: %d") % noutput_items));
  gr_complex const* in = (gr_complex const *)input_items[0];
  gr_complex *out_symb = (gr_complex *)output_items[0];
  gr_complex *out_taps = (gr_complex *)output_items[1];

  const int nin = ninput_items[0];

  assert(ninput_items[0] >= 2*_nGuard + _nB + _nF + 1);
  if (ninput_items[0] < 2*_nGuard + _nB + _nF + 1)
    return 0;
  int const ninput = ninput_items[0] - _nGuard - _nF;

  int nout = 0; // counter for produced output items
  switch (_state) {
    case WAIT_FOR_PREAMBLE: {
      std::vector<tag_t> v;
      get_tags_in_window(v, 0, history()-1, ninput, pmt::intern("preamble_start"));
      if (v.empty()) {
        consume(0, ninput - history()+1);
      } else {
        tag_t const& tag = v.front();
        // uint64_t const offset = tag.offset - nitems_read(0) + history() - 1;
        // std::cout << "========= offset= " << offset
        //           << " tag.offset= " << tag.offset
        //           << " nitems_read(0)= " << nitems_read(0)
        //           << " tag.offset-nitems_read(0)= " << tag.offset - nitems_read(0) << " ==========" << std::endl;
        // for (int k=0; k<ninput; ++k)
        //   std::cout << "SAMPLE: " << k << " " << k-int(offset) << " " << k+nitems_read(0) << " " << in[k] << std::endl;
        reset_filter();
        _descrambled_symbols.clear();
        publish_frame_info();
        consume(0, tag.offset - nitems_read(0));
        _state = WAIT_FOR_FRAME_INFO;
        GR_LOG_DEBUG(d_logger, "got preamble tag > wait for frame info");
      }
      break;
    } // WAIT_FOR_PREAMBLE
    case WAIT_FOR_FRAME_INFO: {
      //GR_LOG_DEBUG(d_logger, "WAIT_FOR_FRAME_INFO");
      //update_frame_info(delete_head_blocking(_msg_ports["frame_info"]));
      break;
    } // WAIT_FOR_FRAME_INFO
    case DO_FILTER: {
      // std::cout << "========= offset (DO_FILTER) nitems_read(0)= " << nitems_read(0) << " ==========" << std::endl;
      _rotated_samples.resize(ninput+_nF+1);
      int ninput_processed = 0;
      for (int i=history()-1; i<ninput && nout<noutput_items; i+=_sps, ninput_processed+=_sps) {
        // rotate samples
        if (i == history()-1) {
          _rotator.rotateN(&_rotated_samples[0] + i - _nB,
                           in      + i - _nB,
                           _nB+_nF+1);
        } else {
          _rotator.rotateN(&_rotated_samples[0] + i + _nF+1 - _sps,
                           in      + i + _nF+1 - _sps,
                           _sps);
        }
        if (_symbol_counter == _symbols.size()) {
          publish_frame_info();
          publish_soft_dec();
          _symbol_counter = 0;

          int const shift = recenter_filter_taps();
          if (shift != 0)
            ninput_processed += shift;

          _state = WAIT_FOR_FRAME_INFO;
          break;
        }
        // std::cout << "FILTER_CHECK: " << i << " " << i-1-_nB << " " << i+_nF << " " << in[i] << std::endl;
        assert(i+_nF < nin && i-1-_nB >= 0);
        out_symb[nout] = filter(&_rotated_samples[0] + i - _nB,
                                &_rotated_samples[0] + i + _nF+1);
        std::memcpy(&out_taps[(_nB+_nF+1)*nout], _taps_samples, (_nB+_nF+1)*sizeof(gr_complex));
        ++nout;
      } // next sample
      consume(0, ninput_processed);
      break;
    } // DO_FILTER
  }
  return nout;
}

bool adaptive_dfe_impl::start()
{
  gr::thread::scoped_lock lock(d_setlock);
  _taps_samples      = (gr_complex*)(volk_malloc((_nB+_nF+1)*sizeof(gr_complex), volk_get_alignment()));
  _last_taps_samples = (gr_complex*)(volk_malloc((_nB+_nF+1)*sizeof(gr_complex), volk_get_alignment()));
  _taps_symbols      = (gr_complex*)(volk_malloc(        _nW*sizeof(gr_complex), volk_get_alignment()));
  _hist_symbols      = (gr_complex*)(volk_malloc(      2*_nW*sizeof(gr_complex), volk_get_alignment()));
  reset_filter();
  GR_LOG_DEBUG(d_logger,str(boost::format("adaptive_dfe_impl::start() nB=%d nF=%d mu=%f alpha=%f")
                             % _nB % _nF % _mu % _alpha));
  return true;
}
bool adaptive_dfe_impl::stop()
{
  gr::thread::scoped_lock lock(d_setlock);
  GR_LOG_DEBUG(d_logger, "adaptive_dfe_impl::stop()");
  VOLK_SAFE_DELETE(_taps_samples);
  VOLK_SAFE_DELETE(_last_taps_samples);
  VOLK_SAFE_DELETE(_taps_symbols);
  VOLK_SAFE_DELETE(_hist_symbols);
  return true;
}

gr_complex adaptive_dfe_impl::filter(gr_complex const* start, gr_complex const* end) {
  assert(end-start == _nB + _nF + 1);

  _num_samples_since_filter_update += _sps;

  // (1) run the filter filter
  gr_complex filter_output(0);
  // (1a) taps_samples
  volk_32fc_x2_dot_prod_32fc(&filter_output,
                             start,
                             _taps_samples,
                             _nB+_nF+1);
  // (1b) taps_symbols
  gr_complex dot_symbols(0);
  gr::digital::constellation_sptr constell = _constellations[_constellation_index];
  bool const update_taps = true; //constell->bits_per_symbol() <= 3;
  if (constell->bits_per_symbol() > 3)
    _use_symbol_taps = false;
  if (_use_symbol_taps) {
    for (int l=0; l<_nW; ++l) {
      assert(_hist_symbol_index+l < 2*_nW);
      dot_symbols += _hist_symbols[_hist_symbol_index+l]*_taps_symbols[l];
    }
    filter_output += dot_symbols;
  }
  assert(_symbol_counter < _symbols.size());

  gr_complex known_symbol = _symbols[_symbol_counter];
  bool const is_known     = std::abs(known_symbol) > 1e-5;
  // (2)  unknown symbols (=data): compute soft decisions
  if (not is_known) {
    gr_complex const descrambled_filter_output = std::conj(_scramble[_symbol_counter]) * filter_output;
    unsigned int const jc = constell->decision_maker(&descrambled_filter_output);
    gr_complex descrambled_symbol = 0;
    constell->map_to_points(jc, &descrambled_symbol);
    if (_save_soft_decisions) {
      float const err = std::abs(descrambled_filter_output - descrambled_symbol);
      std::vector<float> const soft_dec = constell->calc_soft_dec
        (descrambled_filter_output, _npwr[_constellation_index].filter(err));
      std::copy(soft_dec.begin(), soft_dec.end(), std::back_inserter<std::vector<float> >(_vec_soft_decisions));
    }
    known_symbol = _scramble[_symbol_counter] * descrambled_symbol;
  }
  // std::cout << "FILTER: " << filter_output <<" " << known_symbol << " " << start[_nB+1] << std::endl;
  // (3) filter update
  if (is_known || update_taps) {
    // (3a) control loop update for doppler correction using the adaptibve filter taps
    gr_complex acc(0);
    for (int j=_nB+1-2*_sps; j<_nB+1+2*_sps+1; ++j)
      acc += std::conj(_last_taps_samples[j]) * _taps_samples[j];
    float const frequency_err = gr::fast_atan2f(acc)/_num_samples_since_filter_update; // frequency error (rad/sample)
    _control_loop.advance_loop(frequency_err);
    _control_loop.phase_wrap();
    _control_loop.frequency_limit();
    _rotator.set_phase_incr(gr_expj(_control_loop.get_frequency()));

    // (3b) update of adaptive filter taps
    gr_complex const err =  filter_output - known_symbol;
    //       taps_samples
    for (int j=0; j<_nB+_nF+1; ++j) {
      _last_taps_samples[j] = _taps_samples[j];
      _taps_samples[j]     -= _mu*err*std::conj(start[j]);
      _num_samples_since_filter_update = 0;
    }
    //       taps_symbols
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
  // (4) save the descrambled symbol (-> frame_info)
  _descrambled_symbols[_symbol_counter] = filter_output*std::conj(_scramble[_symbol_counter]);
  return _descrambled_symbols[_symbol_counter++];
}

int
adaptive_dfe_impl::recenter_filter_taps() {
  // get max(abs(taps))
  ssize_t const idx_max = std::distance(_taps_samples,
                                        std::max_element(_taps_samples+_nB+1-3*_sps, _taps_samples+_nB+1+3*_sps,
                                                         [](gr_complex a, gr_complex b) {
                                                           return std::norm(a) < std::norm(b);
                                                         }));
  // GR_LOG_DEBUG(d_logger, str(boost::format("idx_max=%2d abs(tap_max)=%f") % idx_max % std::abs(_taps_samples[idx_max])));
  if (idx_max-_nB-1 > +2*_sps) {
    // maximum is right of the center tap
    //   -> shift taps to the left left
    GR_LOG_DEBUG(d_logger, "shift left");
    std::copy(_taps_samples+2*_sps, _taps_samples+_nB+_nF+1, _taps_samples);
    std::fill_n(_taps_samples+_nB+_nF+1-2*_sps, 2*_sps, gr_complex(0));
    return +2*_sps;
  }
  if (idx_max-_nB-1 < -2*_sps) {
    // maximum is left of the center tap
    //   -> shift taps to the right
    GR_LOG_DEBUG(d_logger, "shift right");
    std::copy_backward(_taps_samples, _taps_samples+_nB+_nF+1-2*_sps,
                       _taps_samples+_nB+_nF+1);
    std::fill_n(_taps_samples, 2*_sps, gr_complex(0));
    return -2*_sps;
  }
  return 0;
}

void adaptive_dfe_impl::reset_filter()
{
  std::fill_n(_taps_samples,              _nB+_nF+1, gr_complex(0));
  std::fill_n(_last_taps_samples,         _nB+_nF+1, gr_complex(0));
  std::fill_n(_taps_symbols,                    _nW, gr_complex(0));
  std::fill_n(_hist_symbols,                  2*_nW, gr_complex(0));
  _taps_symbols[0]     = 1;
  _hist_symbol_index   = 0;
  _num_samples_since_filter_update = 0;
}

void adaptive_dfe_impl::publish_frame_info()
{
  pmt::pmt_t data = pmt::make_dict();
  GR_LOG_DEBUG(d_logger, str(boost::format("publish_frame_info %d == %d") % _descrambled_symbols.size() % _symbols.size()));
  data = pmt::dict_add(data, pmt::intern("symbols"), pmt::init_c32vector(_descrambled_symbols.size(), &_descrambled_symbols.front()));
  // for (int i=0; i<_vec_soft_decisions.size(); ++i)
  //   _vec_soft_decisions[i] = std::max(-1.0f, std::min(1.0f, _vec_soft_decisions[i]));
  data = pmt::dict_add(data, pmt::intern("soft_dec"), pmt::init_f32vector(_vec_soft_decisions.size(), &_vec_soft_decisions.front()));
  message_port_pub(_msg_ports["frame_info"], data);
  _descrambled_symbols.clear();
}

void adaptive_dfe_impl::publish_soft_dec()
{
  if (_vec_soft_decisions.empty())
    return;
  message_port_pub(_msg_ports["soft_dec"],
                   pmt::cons(pmt::dict_add(_msg_metadata, pmt::intern("packet_len"), pmt::mp(_vec_soft_decisions.size())),
                             pmt::init_f32vector(_vec_soft_decisions.size(), _vec_soft_decisions)));
  _vec_soft_decisions.clear();
}

void adaptive_dfe_impl::update_constellations(pmt::pmt_t data) {
  int const n = pmt::length(data);
  _constellations.resize(n);
  _npwr.resize(n);
  std::cout << "adaptive_dfe_impl::update_constellations " << data << std::endl;
  std::cout << "adaptive_dfe_impl::update_constellations n=" << n << std::endl;
  unsigned int const rotational_symmetry = 0;
  unsigned int const dimensionality = 1;

  for (int i=0; i<n; ++i) {
    pmt::pmt_t c = pmt::vector_ref(data, i);
    int const idx = pmt::to_long(pmt::dict_ref(c, pmt::intern("idx"), pmt::from_long(-1)));
    assert(idx>=0 && idx < n);
    _constellations[idx] = gr::digital::constellation_calcdist::make
      (pmt::c32vector_elements(pmt::dict_ref(c, pmt::intern("points"),  pmt::PMT_NIL)),
       pmt::s32vector_elements(pmt::dict_ref(c, pmt::intern("symbols"), pmt::PMT_NIL)),
       rotational_symmetry, dimensionality);
    _npwr[i].reset(_npwr_max_time_constant);
  }
}

void adaptive_dfe_impl::update_frame_info(pmt::pmt_t data)
{
  //GR_LOG_DEBUG(d_logger,str(boost::format("adaptive_dfe_impl::update_frame_info() %s") % data));
  _symbols    = pmt::c32vector_elements(pmt::dict_ref(data, pmt::intern("symb"),     pmt::PMT_NIL));
  _scramble   = pmt::c32vector_elements(pmt::dict_ref(data, pmt::intern("scramble"), pmt::PMT_NIL));
  _constellation_index   = pmt::to_long(pmt::dict_ref(data, pmt::intern("constellation_idx"), pmt::PMT_NIL));
  _save_soft_decisions   = pmt::to_bool(pmt::dict_ref(data, pmt::intern("save_soft_dec"),     pmt::PMT_F));
  bool const do_continue = pmt::to_bool(pmt::dict_ref(data, pmt::intern("do_continue"),       pmt::PMT_F));
  assert(_symbols.size() == _scramble.size());
  _descrambled_symbols.resize(_symbols.size());
  _vec_soft_decisions.clear();
  _symbol_counter = 0;
  _state = (do_continue ? DO_FILTER : WAIT_FOR_PREAMBLE);
}

} /* namespace digitalhf */
} /* namespace gr */
