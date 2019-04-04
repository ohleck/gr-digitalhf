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


#ifndef INCLUDED_DIGITALHF_VITERBI48_H
#define INCLUDED_DIGITALHF_VITERBI48_H

#include <digitalhf/api.h>
#include <pmt/pmt.h>

namespace gr {
namespace digitalhf {

/*!
 * \brief <+description of block+>
 * \ingroup digitalhf
 *
 */
class DIGITALHF_API viterbi48
{
  public:
  typedef boost::shared_ptr<viterbi48> sptr;

  virtual ~viterbi48() {}
  /*!
   * \brief Return a shared_ptr to a new instance of digitalhf::viterbi48.
   *
   * To avoid accidental use of raw pointers, digitalhf::viterbi48's
   * constructor is in a private implementation
   * class. digitalhf::viterbi48::make is the public interface for
   * creating new instances.
   */
  static sptr make(std::uint32_t pol0=0xB9,
                   std::uint32_t pol1=0x9D,
                   std::uint32_t pol2=0xD3,
                   std::uint32_t pol3=0xF7);

  virtual void  reset() = 0;
  virtual const std::vector<std::uint8_t>& udpate(const std::vector<float>& soft_dec) = 0;
  virtual float quality() = 0;
} ;

} // namespace digitalhf
} // namespace gr

#endif /* INCLUDED_DIGITALHF_VITERBI48_H */
