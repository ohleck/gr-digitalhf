#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2018 hcab14@gmail.com.
#
# This is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3, or (at your option)
# any later version.
#
# This software is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this software; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin Street,
# Boston, MA 02110-1301, USA.
#

import importlib
from gnuradio import blocks
from gnuradio import digital
from gnuradio import filter
from gnuradio import gr
import pmt
import digitalhf
import digitalhf.physical_layer

class physical_layer_driver(gr.hier_block2):
    """
    docstring for block physical_layer_driver
    """
    def __init__(self, samp_rate,sps,alpha,mu,nB,nF,nW,description_name,mode):
        gr.hier_block2.__init__(self,
                                "physical_layer_driver",
                                gr.io_signature(1, 1, gr.sizeof_gr_complex), # Input signature
                                gr.io_signature3(3, 3,
                                                 gr.sizeof_gr_complex,
                                                 gr.sizeof_gr_complex,
                                                 gr.sizeof_gr_complex*(1+sps*(nB+nF)))) # Output signature

        self._sps   = sps
        self._alpha = alpha
        self._mu    = mu
        self._nB    = nB
        self._nF    = nF
        self._nW    = nW

        m = importlib.import_module('digitalhf.physical_layer.'+description_name)
        self._physical_layer_driver_description = m.PhysicalLayer(sps)
        self._physical_layer_driver_description.set_mode(mode)

        ## TODO: get rrc tap information from physical layer description
        self._rrc_taps = filter.firdes.root_raised_cosine(1.0, samp_rate, samp_rate/sps, 0.35, 11*sps)
        preamble_offset,preamble_samples = self._physical_layer_driver_description.get_preamble_z()
        preamble_length          = sps*len(self._physical_layer_driver_description.get_preamble())
        self._rrc_filter         = filter.fir_filter_ccc(1, (self._rrc_taps))
        self._corr_est           = digital.corr_est_cc(symbols    = (preamble_samples.tolist()),
                                                       sps        = sps,
                                                       mark_delay = preamble_offset,
                                                       threshold  = 0.5,
                                                       threshold_method = 1)
        self._doppler_correction = digitalhf.doppler_correction_cc(preamble_length, len(preamble_samples))
        self._adaptive_filter    = digitalhf.adaptive_dfe(sps, nB, nF, nW, mu, alpha)
        self._msg_proxy          = digitalhf.msg_proxy(self._physical_layer_driver_description)
        self.connect((self, 0),
                     (self._rrc_filter, 0),
                     (self._corr_est, 0),
                     (self._doppler_correction, 0),
                     (self._adaptive_filter, 0),
                     (self, 0))
        self.connect((self._corr_est, 1),        ## correlation
                     (self, 1))
        self.connect((self._adaptive_filter, 1), ## taps
                     (self, 2))

        self.msg_connect((self._doppler_correction, 'doppler'), (self._msg_proxy, 'doppler'))
        self.msg_connect((self._msg_proxy, 'doppler'), (self._doppler_correction, 'doppler'))

        self.msg_connect((self._adaptive_filter, 'frame_info'), (self._msg_proxy, 'frame_info'))
        self.msg_connect((self._msg_proxy, 'frame_info'), (self._adaptive_filter, 'frame_info'))

        constellations_data = self._physical_layer_driver_description.get_constellations()
        constellations_msg  = pmt.to_pmt([{'idx': idx, 'points': c['points'], 'symbols': c['symbols']}
                                          for (idx,c) in enumerate(constellations_data)])
        self._adaptive_filter.to_basic_block()._post(pmt.intern('constellations'), constellations_msg)

        self.message_port_register_hier_out('soft_dec')
        self.msg_connect((self._adaptive_filter, 'soft_dec'), (self, 'soft_dec'))
        self.msg_connect((self._msg_proxy, 'soft_dec'), (self, 'soft_dec'))

    def set_mu(self, mu):
        self._adaptive_filter.set_mu(mu)

    def set_alpha(self, alpha):
        self._adaptive_filter.set_alpha(alpha)

    def set_mode(self, mode):
        self._physical_layer_driver_description.set_mode(mode)
