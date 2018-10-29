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

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <gnuradio/digital/constellation.h>
#include <digitalhf/adaptive_dfe.h>

namespace gr {
namespace digitalhf {

class adaptive_dfe_impl : public adaptive_dfe {
private:
  int _sps;
  int _nB, _nF, _nW;

  float _mu;
  float _alpha;

  // module name w.r.t. digitalhf.physical_layer containing a PhysicalLayer class
  std::string           _py_module_name;
  boost::python::object _physicalLayer; // class instance of physical layer description

  gr_complex* _taps_samples;
  gr_complex* _taps_symbols;

  gr_complex* _hist_samples;
  gr_complex* _hist_symbols;

  int _hist_sample_index;
  int _hist_symbol_index;

  std::size_t _sample_counter;

  std::vector<gr::digital::constellation_sptr> _constellations;
  std::vector<float> _npwr;
  std::vector<int>   _npwr_counter;
  int _npwr_max_time_constant;
  int _constellation_index;
  std::vector<gr_complex> _samples;
  std::vector<gr_complex> _symbols;
  std::vector<gr_complex> _scramble;
  std::vector<gr_complex> _descrambled_symbols;
  int _symbol_counter;

  bool _need_samples;
  bool _save_soft_decisions;
  std::vector<float> _vec_soft_decisions;
  pmt::pmt_t         _msg_port_name;
  pmt::pmt_t         _msg_metadata;

  // PLL for doppler tracking
  float _df;     // frequency offset in radians per sample
  float _phase;  // accumulated phase for frequency correction
  const float _b[2];
  float _ud;

  enum state {
    WAIT_FOR_PREAMBLE,
    INITIAL_DOPPLER_ESTIMATE,
    INITIAL_DOPPLER_ESTIMATE_CONTINUE,
    DO_FILTER
  } _state;

  void update_constellations(boost::python::object obj);
  bool update_frame_information(boost::python::object obj);
  bool update_doppler_information(boost::python::object obj);

  void update_local_oscillator();
  gr_complex filter();

  void insert_sample(gr_complex z);
  void update_pll(float doppler);
  bool get_correlation_tag(uint64_t i, uint64_t& offset, float& phase_est);

public:
  adaptive_dfe_impl(int sps, // samples per symbol
                    int nB,  // number of forward FIR taps
                    int nF,  // number of backward FIR taps
                    int nW,  // number of symbol taps
                    float mu,
                    float alpha,
                    std::string physical_layer_description);
  virtual ~adaptive_dfe_impl();

  void forecast (int noutput_items, gr_vector_int &ninput_items_required);

  virtual bool start();
  virtual bool stop();

  virtual int general_work(int noutput_items,
                           gr_vector_int &ninput_items,
                           gr_vector_const_void_star &input_items,
                           gr_vector_void_star &output_items);

  virtual void set_mu(float mu) { _mu = mu; }
  virtual void set_mode(std::string);

} ;

} // namespace digitalhf
} // namespace gr

#endif /* INCLUDED_DIGITALHF_ADAPTIVE_DFE_IMPL_H */
