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

namespace gr {
namespace digitalhf {

#include <array>
#include <cmath>
#include <cstdint>

template<size_t N>
class llr_to_prob {
public:
  llr_to_prob()
    : _table() {
    for (int i=0; i<N; ++i) { // llr -> probability(1)
      float const x = -7.0f + 14.0f*float(i)/float(N-1);
      _table[i] = std::uint8_t(0.5f + 255.0f/(1.0f + std::exp(x)));
    }
  }
  virtual ~llr_to_prob() {}

  std::uint8_t table_lookup(float llr) const {
    float const x = (7.0f + std::max(-7.0f, std::min(+7.0f, llr)))/14.0f; // \in [0,1]
    return _table[int(0.5 + (N-1)*x)];
  }
protected:
private:
  std::array<std::uint8_t, N> _table;
} ;

} // namespace digitalhf
} // namespace gr
