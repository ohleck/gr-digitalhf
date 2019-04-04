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

#ifndef INCLUDED_DIGITALHF_VITERBI29IMPL_H
#define INCLUDED_DIGITALHF_VITERBI29IMPL_H

#include <digitalhf/viterbi29.h>

#include "viterbi.h"
#include "llr_to_prob.h"

namespace gr {
namespace digitalhf {

class viterbi29_impl : public viterbi29, public llr_to_prob<1024> {
private:
  viterbi<2,9> _v;
  float  _quality;
  size_t _num_bits;
  std::vector<std::uint8_t> _bits;

public:
  viterbi29_impl(std::uint32_t pol0, std::uint32_t pol1);
  virtual ~viterbi29_impl();

  virtual void       reset();
  virtual const std::vector<std::uint8_t>& udpate(const std::vector<float>& soft_dec);
  virtual float      quality();
} ;

} // namespace digitalhf
} // namespace gr

#endif /* INCLUDED_DIGITALHF_VITERBI29IMPL_H */
