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

#include "viterbi48_impl.h"

namespace gr {
namespace digitalhf {

viterbi48::sptr
viterbi48::make(std::uint32_t pol0, std::uint32_t pol1,
                std::uint32_t pol2, std::uint32_t pol3)
{
  return sptr(new viterbi48_impl(pol0, pol1, pol2, pol3));
}

viterbi48_impl::viterbi48_impl(std::uint32_t pol0, std::uint32_t pol1,
                               std::uint32_t pol2, std::uint32_t pol3)
  : _v({pol0, pol1, pol2, pol3})
  , _quality(0)
  , _num_bits(0)
  , _bits()
{
}

viterbi48_impl::~viterbi48_impl()
{
}

void viterbi48_impl::reset() {
  _v.reset();
  _quality = 0;
}


const std::vector<std::uint8_t>& viterbi48_impl::udpate(const std::vector<float>& soft_dec) {
  size_t symb_len = soft_dec.size();
  float const *sd = &soft_dec[0];

  size_t const num_bits = symb_len/4;
  _bits.resize(num_bits);
  std::vector<uint8_t>::iterator iterator_bits = _bits.begin();

  size_t bits_per_frame = 5*8; // 5 times constraint length
  _v.resize(bits_per_frame);

  size_t num_frames = num_bits/bits_per_frame;

  _quality = 0;
  int i = 0;
  for (; i<num_frames; ++i) {
    for (int j=0; j<bits_per_frame; ++j) {
      std::array<std::uint8_t, 4> const symb = {
        table_lookup(sd[4*(i*bits_per_frame+j)  ]),
        table_lookup(sd[4*(i*bits_per_frame+j)+1]),
        table_lookup(sd[4*(i*bits_per_frame+j)+2]),
        table_lookup(sd[4*(i*bits_per_frame+j)+3])
      };
      _v.update(j, symb);
    }
    _quality      += _v.chainback(iterator_bits, iterator_bits+bits_per_frame);
    iterator_bits += bits_per_frame;
  }

  size_t remaining_bits = num_bits - bits_per_frame*num_frames;
  _v.resize(remaining_bits);
  for (int j=0; j<remaining_bits; ++j) {
    std::array<std::uint8_t, 4> const symb = {
      table_lookup(sd[4*(i*bits_per_frame+j)  ]),
      table_lookup(sd[4*(i*bits_per_frame+j)+1]),
      table_lookup(sd[4*(i*bits_per_frame+j)+2]),
      table_lookup(sd[4*(i*bits_per_frame+j)+3])
    };
    _v.update(j, symb);
  }
  _quality += _v.chainback(iterator_bits, iterator_bits+remaining_bits);
  return _bits;
}

float viterbi48_impl::quality() {
  return _quality;
}

} // namespace digitalhf
} // namespace gr
