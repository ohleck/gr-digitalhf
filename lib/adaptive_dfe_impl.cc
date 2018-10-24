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
#include <gnuradio/digital/constellation.h>
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
                   std::string python_file_name)
{
  return gnuradio::get_initial_sptr
    (new adaptive_dfe_impl(sps, nB, nF, nW, python_file_name));
}

/*
 * The private constructor
 */
adaptive_dfe_impl::adaptive_dfe_impl(int sps, // samples per symbol
                                     int nB,  // number of forward FIR taps
                                     int nF,  // number of backward FIR taps
                                     int nW,  // number of feedback taps
                                     std::string python_file_name)
  : gr::block("adaptive_dfe",
              gr::io_signature::make(1, 1, sizeof(gr_complex)),
              gr::io_signature::make(0, 0, sizeof(gr_complex)))
  , _sps(sps)
  , _nB(nB)
  , _nF(nF)
  , _nW(nW)
  , _py_file_name(python_file_name)
  , _physicalLayer()
  , _taps_samples(nB+nF+1)
  , _taps_symbols(nW) {
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
  const gr_complex *in = (const gr_complex *) input_items[0];

  get_next_frame();
  GILLock lock;
  std::cout << "bits_per_symbol: " << boost::python::extract<int>(_constellation.attr("bits_per_symbol")()) << std::endl;

  consume_each (noutput_items);

  // Tell runtime system how many output items we produced.
  return noutput_items;
}

void adaptive_dfe_impl::get_next_frame()
{
  GILLock lock;
  boost::python::object const& obj = _physicalLayer.attr("get_frame")();
  std::cout << "get_frame" << std::endl;
  boost::python::numpy::ndarray  symbols = boost::python::numpy::array(obj[0]);
  _constellation = obj[1];
}
boost::python::object import(const std::string& module, const std::string& path, boost::python::object& globals)
{
    boost::python::dict locals;
    locals["module_name"] = module;
    locals["path"]        = path;

    boost::python::exec("import imp\n"
             "new_module = imp.load_module(module_name, open(path), path, ('py', 'U', imp.PY_SOURCE))\n",
             globals,
             locals);
    return locals["new_module"];
}
bool adaptive_dfe_impl::start()
{
  std::cout << "adaptive_dfe_impl::start()" << std::endl;
  GILLock lock;
  try {
    boost::python::object main          = boost::python::import("__main__");
    boost::python::object globals       = main.attr("__dict__");
    boost::python::object module        = import("physicalLayer", _py_file_name, globals);
    boost::python::object PhysicalLayer = module.attr("PhysicalLayer");
    _physicalLayer = PhysicalLayer();

    _physicalLayer.attr("get_frame")();
  } catch (const boost::python::error_already_set& ) {
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

} /* namespace digitalhf */
} /* namespace gr */
