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

#include <gnuradio/digital/constellation.h>
#include <digitalhf/adaptive_dfe.h>

namespace gr {
namespace digitalhf {

class constellation_distance_filter {
public:
  constellation_distance_filter(int max_time_constant=10)
    : _max_time_constant(max_time_constant)
    , _counter(0)
    , _pwr(0) {}

  void reset(int max_time_constant) {
    _max_time_constant = max_time_constant;
    _counter = 0;
    _pwr = 0;
  }
  float filter(float x) {
    _counter += (_counter < _max_time_constant);
    float const alpha = 1.0f/_counter;
    _pwr = (1-alpha)*_pwr + alpha*x;
    return _pwr;
  }
protected:
private:
  int   _max_time_constant;
  int   _counter;
  float _pwr; // filtered distance to constellation point
} ;

class adaptive_dfe_impl : public adaptive_dfe {
private:
  int _sps;
  int _nB, _nF, _nW;
  int _nGuard;

  float _mu;
  float _alpha;

  bool _use_symbol_taps;

  // module name w.r.t. digitalhf.physical_layer containing a PhysicalLayer class
  // std::string           _py_module_name;
  // boost::python::object _physicalLayer; // class instance of physical layer description

  gr_complex* _taps_samples;
  gr_complex* _taps_symbols;

  gr_complex* _hist_symbols;
  int _hist_symbol_index;

  std::vector<gr::digital::constellation_sptr> _constellations;
  std::vector<constellation_distance_filter> _npwr;
  int _npwr_max_time_constant;
  int _constellation_index;
  std::vector<gr_complex> _symbols;
  std::vector<gr_complex> _scramble;
  std::vector<gr_complex> _descrambled_symbols;
  int _symbol_counter;

  bool _save_soft_decisions;
  std::vector<float> _vec_soft_decisions;
  std::map<std::string, pmt::pmt_t> _msg_ports;
  pmt::pmt_t _msg_metadata;

  enum state {
    WAIT_FOR_PREAMBLE,
    WAIT_FOR_FRAME_INFO,
    DO_FILTER
  } _state;

//  void update_constellations(boost::python::object obj);
  void update_constellations(pmt::pmt_t );
  void update_frame_info(pmt::pmt_t );

  gr_complex filter(gr_complex const* start, gr_complex const* end);
  int recenter_filter_taps();
  void reset_filter();

  void publish_frame_info();
  void publish_soft_dec();

public:
  adaptive_dfe_impl(int sps, // samples per symbol
                    int nB,  // number of forward FIR taps
                    int nF,  // number of backward FIR taps
                    int nW,  // number of symbol taps
                    float mu,
                    float alpha);
  virtual ~adaptive_dfe_impl();

  void forecast (int noutput_items, gr_vector_int &ninput_items_required);

  virtual bool start();
  virtual bool stop();

  virtual int general_work(int noutput_items,
                           gr_vector_int &ninput_items,
                           gr_vector_const_void_star &input_items,
                           gr_vector_void_star &output_items);

  virtual void set_mu(float mu) { _mu = mu; }
  virtual void set_alpha(float alpha) { _alpha = alpha; }
} ;

} // namespace digitalhf
} // namespace gr

#endif /* INCLUDED_DIGITALHF_ADAPTIVE_DFE_IMPL_H */
