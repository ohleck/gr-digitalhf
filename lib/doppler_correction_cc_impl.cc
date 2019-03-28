/* -*- c++ -*- */
/*
 * Copyright 2018 hcab14@gmail.com.
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
#include <gnuradio/expj.h>
#include <gnuradio/logger.h>

#include "doppler_correction_cc_impl.h"

namespace gr {
namespace digitalhf {

doppler_correction_cc::sptr
doppler_correction_cc::make(unsigned int preamble_length, unsigned int preamble_length_cc)
{
  return gnuradio::get_initial_sptr
    (new doppler_correction_cc_impl(preamble_length, preamble_length_cc));
}

/*
 * The private constructor
 */
doppler_correction_cc_impl::doppler_correction_cc_impl(unsigned int preamble_length,
                                                       unsigned int preamble_length_cc)
  : gr::sync_block("doppler_correction_cc",
                   gr::io_signature::make(1, 1, sizeof(gr_complex)),
                   gr::io_signature::make(1, 1, sizeof(gr_complex)))
  , _preamble_length(preamble_length)
  , _preamble_length_cc(preamble_length_cc)
  , _rotator()
  , _state(WAIT_FOR_PHASE_EST_TAG)
  , _msg_metadata(pmt::make_dict())
  , _port_name(pmt::intern("doppler"))
  , _phase_est(0)
{
  GR_LOG_DECLARE_LOGPTR(d_logger);
  GR_LOG_ASSIGN_LOGPTR(d_logger, "doppler_correction_cc");
  message_port_register_out(_port_name);
  message_port_register_in (_port_name);
  set_msg_handler(_port_name, boost::bind(&doppler_correction_cc_impl::handle_message, this, _1));
  set_tag_propagation_policy(TPP_DONT);
}

doppler_correction_cc_impl::~doppler_correction_cc_impl()
{
}

void
doppler_correction_cc_impl::handle_message(pmt::pmt_t msg)
{
  gr::thread::scoped_lock lock(d_setlock);
  bool const success = pmt::to_bool(pmt::dict_ref(msg, pmt::intern("success"), pmt::get_PMT_F()));
  if (!success) {
    // GR_LOG_DEBUG(d_logger, "next state > CONSUME_AND_SKIP success=false");
    if (_state == WAIT_FOR_MSG)
      _state = CONSUME_AND_SKIP;
    return;
  }
  float const doppler = pmt::to_float(pmt::dict_ref(msg, pmt::intern("doppler"), pmt::from_float(0)));
  _rotator.set_phase_incr(gr_expj(-doppler));
  if (_state == WAIT_FOR_MSG) {
    _rotator.set_phase(gr_expj(-_phase_est + 0.5*doppler*_preamble_length_cc));
    _state = CONSUME_AND_INSERT_PREAMBLE_TAG;
    GR_LOG_DEBUG(d_logger, str(boost::format("next state > CONSUME_AND_INSERT_PREAMBLE_TAG phase_est=%f doppler=%f")
                               % (_phase_est - 0.5*doppler*_preamble_length_cc)
                               % doppler));
  }
}

void
doppler_correction_cc_impl::forecast(int noutput_items, gr_vector_int &ninput_items_required)
{
  ninput_items_required[0] = _preamble_length + 1;
//  GR_LOG_DEBUG(d_logger, str(boost::format("forecast: %d %d %d") % noutput_items % ninput_items_required[0] % _preamble_length));
}

int
doppler_correction_cc_impl::work(int noutput_items,
                                 gr_vector_const_void_star &input_items,
                                 gr_vector_void_star &output_items)
{
  gr::thread::scoped_lock lock(d_setlock);
  gr_complex const* in = (gr_complex const*)input_items[0];
  gr_complex *out = (gr_complex *) output_items[0];
  // GR_LOG_DEBUG(d_logger, str(boost::format("work: %d %d") % noutput_items % _preamble_length));
  if (noutput_items < _preamble_length)
    return 0;
  noutput_items -= _preamble_length;
  int nout = 0;
  switch (_state) {
    case WAIT_FOR_PHASE_EST_TAG: {
      std::vector<tag_t> v;
      get_tags_in_window(v, 0, 0, noutput_items, pmt::intern("phase_est"));
      if (v.empty()) {
        nout = noutput_items;
      } else {
        tag_t const& tag = v.front();
        uint64_t const offset = tag.offset - nitems_read(0);
        nout = offset;
        _phase_est = pmt::to_double(tag.value);
        _msg_metadata = pmt::dict_add(_msg_metadata, pmt::intern("packet_len"), pmt::from_long(_preamble_length));
        message_port_pub(_port_name,
                         pmt::cons(_msg_metadata,
                                   pmt::init_c32vector(_preamble_length, in+nout)));
        // GR_LOG_DEBUG(d_logger, str(boost::format("next state > WAIT_FOR_MSG %lld phase_est=%f") % tag.offset % _phase_est));
        _state = WAIT_FOR_MSG;
      }
      break;
    } // WAIT_FOR_PHASE_EST_TAG
    case WAIT_FOR_MSG: {
      // GR_LOG_DEBUG(d_logger, "WAIT_FOR_MSG");
      // handle_message(delete_head_nowait(_port_name));
      break;
    } // WAIT_FOR_MSG
    case CONSUME_AND_INSERT_PREAMBLE_TAG: {
      add_item_tag(0, nitems_read(0), pmt::intern("preamble_start"), pmt::from_long(0));
      nout = _preamble_length;
      // GR_LOG_DEBUG(d_logger, str(boost::format("next state > WAIT_FOR_PHASE_EST_TAG %lld") % nitems_read(0)));
      _state = WAIT_FOR_PHASE_EST_TAG;
      break;
    } // CONSUME_AND_INSERT_PREAMBLE_TAG
    case CONSUME_AND_SKIP: {
      nout = _preamble_length;
      // GR_LOG_DEBUG(d_logger, str(boost::format("next state > WAIT_FOR_PHASE_EST_TAG %lld") % nitems_read(0)));
      _state = WAIT_FOR_PHASE_EST_TAG;
      break;
    } // CONSUME_AND_SKIP
  }
  // apply current doppler correction to all produced samples
  _rotator.rotateN(out, in, nout);
  // Tell runtime system how many output items we produced.
  return nout;
}

} /* namespace digitalhf */
} /* namespace gr */
