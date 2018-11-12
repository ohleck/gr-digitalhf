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

#ifndef INCLUDED_DIGITALHF_DOPPLER_CORRECTION_CC_IMPL_H
#define INCLUDED_DIGITALHF_DOPPLER_CORRECTION_CC_IMPL_H

#include <gnuradio/blocks/rotator.h>

#include <digitalhf/doppler_correction_cc.h>

namespace gr {
namespace digitalhf {

class doppler_correction_cc_impl : public doppler_correction_cc
{
private:
  unsigned int    _preamble_length;    // length of preamble (in samples)
  unsigned int    _preamble_length_cc; // length of the part of the preamble used for cross correlation (in samples)
  blocks::rotator _rotator;
  enum {
    WAIT_FOR_PHASE_EST_TAG,   // wait for a tag from corr_est_cc
    WAIT_FOR_MSG,             // wait for response from msg_proxy
    CONSUME_AND_INSERT_PREAMBLE_TAG, //              insert a preamble tag (=doppler calculation was successful)
    CONSUME_AND_SKIP          // do not intsert a tag and skip the samples (=doppler calculation was not successful)
  } _state;
  pmt::pmt_t _msg_metadata;
  pmt::pmt_t _port_name;

  float _phase_est;

public:
  doppler_correction_cc_impl(unsigned int preamble_length, unsigned int preamble_length_cc);
  virtual ~doppler_correction_cc_impl();

  void forecast(int noutput_items, gr_vector_int &ninput_items_required);

  int work(int noutput_items,
           gr_vector_const_void_star &input_items,
           gr_vector_void_star &output_items);
protected:
  void handle_message(pmt::pmt_t msg);

};

} // namespace digitalhf
} // namespace gr

#endif /* INCLUDED_DIGITALHF_DOPPLER_CORRECTION_CC_IMPL_H */
