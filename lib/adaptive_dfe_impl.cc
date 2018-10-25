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

#include <gnuradio/io_signature.h>
#include "adaptive_dfe_impl.h"

namespace gr {
namespace digitalhf {

namespace {
class GILLock {
  PyGILState_STATE _state;
  public:
  GILLock()
    :_state(PyGILState_Ensure()) {}
  ~GILLock() {
    PyGILState_Release(_state);
  }
} ;
}
adaptive_dfe::sptr
adaptive_dfe::make(int sps, // samples per symbol
                   int nB,  // number of forward FIR taps
                   int nF,  // number of backward FIR taps
                   int nW,  // number of feedback taps
                   std::string python_module_name)
{
  return gnuradio::get_initial_sptr
    (new adaptive_dfe_impl(sps, nB, nF, nW, python_module_name));
}

/*
 * The private constructor
 */
adaptive_dfe_impl::adaptive_dfe_impl(int sps, // samples per symbol
                                     int nB,  // number of forward FIR taps
                                     int nF,  // number of backward FIR taps
                                     int nW,  // number of feedback taps
                                     std::string python_module_name)
  : gr::block("adaptive_dfe",
              gr::io_signature::make(1, 1, sizeof(gr_complex)),
              gr::io_signature::make(0, 0, sizeof(gr_complex)))
  , _sps(sps)
  , _nB(nB)
  , _nF(nF)
  , _nW(nW)
  , _py_module_name(python_module_name)
  , _physicalLayer()
  , _taps_samples(nB+nF+1)
  , _taps_symbols(nW)
  , _constellations()
  , _constellation_index()
  , _symbols()
  , _scramble()
{
  // make sure python is ready for threading
  if( Py_IsInitialized() ){
    if(PyEval_ThreadsInitialized() != 1 ){
      PyEval_InitThreads();
    }
    boost::python::numpy::initialize();
  } else {
    throw std::runtime_error("dont use es_pyhandler without python!");
  }
}

/*
 * Our virtual destructor.
 */
adaptive_dfe_impl::~adaptive_dfe_impl()
{
}

void
adaptive_dfe_impl::forecast (int noutput_items, gr_vector_int &ninput_items_required)
{
  /* <+forecast+> e.g. ninput_items_required[0] = noutput_items */
}

int
adaptive_dfe_impl::general_work(int noutput_items,
                                gr_vector_int &ninput_items,
                                gr_vector_const_void_star &input_items,
                                gr_vector_void_star &output_items)
{
  gr_complex const* in = (gr_complex const *)input_items[0];

  GILLock lock;
  // TODO: wait for preamble correlation tag etc...
  update_frame_information(_physicalLayer.attr("get_frame")());
  update_doppler_information(_physicalLayer.attr("get_doppler")()); // symbols

  consume_each (noutput_items);

  // Tell runtime system how many output items we produced.
  return noutput_items;
}

bool adaptive_dfe_impl::start()
{
  std::cout << "adaptive_dfe_impl::start()" << std::endl;
  GILLock lock;
  try {
    boost::python::object module        = boost::python::import(boost::python::str("digitalhf.physical_layer." + _py_module_name));
    boost::python::object PhysicalLayer = module.attr("PhysicalLayer");
    _physicalLayer = PhysicalLayer();
    update_constellations(_physicalLayer.attr("get_constellations")());
  } catch (boost::python::error_already_set const&) {
    PyErr_Print();
    return false;
  }
  return true;
}
bool adaptive_dfe_impl::stop()
{
  std::cout << "adaptive_dfe_impl::stop()" << std::endl;
  GILLock lock;
  _physicalLayer = boost::python::object();
  return true;
}

void adaptive_dfe_impl::update_constellations(boost::python::object obj)
{
  int const n = boost::python::extract<int>(obj.attr("__len__")());
  _constellations.resize(n);
  for (int i=0; i<n; ++i) {
    boost::python::numpy::ndarray const& array = boost::python::numpy::array(obj[i]);
    char const* data = array.get_data();
    int  const     m = array.shape(0);
    std::vector<gr_complex> constell(m);
    std::vector<int>   pre_diff_code(m);
    for (int j=0; j<m; ++j) {
      std::memcpy(&constell[j], data+9*j, sizeof(gr_complex));
      pre_diff_code[j] = (data+9*j)[8];
    }
    unsigned int const rotational_symmetry = 0;
    unsigned int const dimensionality = 1;
    _constellations[i] = gr::digital::constellation_calcdist::make(constell, pre_diff_code, rotational_symmetry, dimensionality);
  }
}
void adaptive_dfe_impl::update_frame_information(boost::python::object obj)
{
  int const n = boost::python::extract<int>(obj.attr("__len__")());
  assert(n==2);
  boost::python::numpy::ndarray array = boost::python::numpy::array(obj[0]);
  char const* data = array.get_data();
  int  const     m = array.shape(0);
  _symbols.resize(m);
  _scramble.resize(m);
  for (int i=0; i<m; ++i) {
    std::memcpy(&_symbols[i],  data+16*i,   sizeof(gr_complex));
    std::memcpy(&_scramble[i], data+16*i+8, sizeof(gr_complex));
  }
  _constellation_index = boost::python::extract<int>(obj[1]);
}
void adaptive_dfe_impl::update_doppler_information(boost::python::object obj)
{
  int const n = boost::python::extract<int>(obj.attr("__len__")());
  assert(n==2);
  double const do_continue = boost::python::extract<bool>(obj[0]);
  double const doppler = boost::python::extract<float>(obj[1]);
  // TODO
}

} /* namespace digitalhf */
} /* namespace gr */
