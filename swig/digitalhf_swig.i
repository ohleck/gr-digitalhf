/* -*- c++ -*- */

#define DIGITALHF_API

%include "gnuradio.i"			// the common stuff

//load generated python docstrings
%include "digitalhf_swig_doc.i"

%{
#include "digitalhf/adaptive_dfe.h"
%}

%include "digitalhf/adaptive_dfe.h"
GR_SWIG_BLOCK_MAGIC2(digitalhf, adaptive_dfe);
