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

#ifndef INCLUDED_DIGITALHF_ADAPTIVE_DFE_IMPL_H
#define INCLUDED_DIGITALHF_ADAPTIVE_DFE_IMPL_H

#include <digitalhf/adaptive_dfe.h>

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

namespace gr {
namespace digitalhf {

class adaptive_dfe_impl : public adaptive_dfe {
private:
  int _sps;
  int _nB, _nF, _nW;

  std::string           _py_file_name;
  boost::python::object _physicalLayer;

  std::vector<gr_complex> _taps_samples;
  std::vector<gr_complex> _taps_symbols;

//  boost::python::numpy::ndarray _symbols;
  boost::python::object _constellation;

  void get_next_frame();

public:
  adaptive_dfe_impl(int sps, // samples per symbol
                    int nB,  // number of forward FIR taps
                    int nF,  // number of backward FIR taps
                    int nW,  // number of symbol taps
                    std::string python_file_name);
  virtual ~adaptive_dfe_impl();

  // Where all the action really happens
  void forecast (int noutput_items, gr_vector_int &ninput_items_required);

  virtual bool start();
  virtual bool stop();

  virtual int general_work(int noutput_items,
                           gr_vector_int &ninput_items,
                           gr_vector_const_void_star &input_items,
                           gr_vector_void_star &output_items);

} ;

} // namespace digitalhf
} // namespace gr

#endif /* INCLUDED_DIGITALHF_ADAPTIVE_DFE_IMPL_H */
