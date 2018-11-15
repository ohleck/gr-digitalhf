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

import numpy as np
from gnuradio import gr
import pmt

class msg_proxy(gr.basic_block):
    """
    docstring for block msg_proxy
    """
    def __init__(self, physical_layer_object):
        gr.basic_block.__init__(self,
                                name="msg_proxy",
                                in_sig=[],
                                out_sig=[])
        self._obj = physical_layer_object

        self._port_doppler = pmt.intern("doppler")
        self.message_port_register_in(self._port_doppler)
        self.message_port_register_out(self._port_doppler)
        self.set_msg_handler(self._port_doppler, self.msg_handler_doppler)

        self._port_frame_info = pmt.intern("frame_info")
        self.message_port_register_in(self._port_frame_info)
        self.message_port_register_out(self._port_frame_info)
        self.set_msg_handler(self._port_frame_info, self.msg_handler_frame)

    def msg_handler_doppler(self, msg_in):
        ## print('-------------------- msg_handler_doppler --------------------')
        iq_samples = pmt.to_python(pmt.cdr(msg_in))
        success,doppler = self._obj.get_doppler(iq_samples)
        msg_out = pmt.make_dict()
        msg_out = pmt.dict_add(msg_out, pmt.intern('success'), pmt.to_pmt(np.bool(success)))
        msg_out = pmt.dict_add(msg_out, pmt.intern('doppler'), pmt.to_pmt(doppler))
        ## print(msg_out)
        self.message_port_pub(self._port_doppler, msg_out)

    def msg_handler_frame(self, msg_in):
        ## print('-------------------- msg_handler_frame --------------------')
        ## print(msg_in)
        symbols = pmt.to_python(pmt.dict_ref(msg_in, pmt.intern('symbols'), pmt.PMT_NIL))
        symb,constellation_idx,do_continue,save_soft_dec = self._obj.get_next_frame(symbols)
        ##print('symb=', symb, symb['symb'], symb['scramble'])
        msg_out = pmt.make_dict()
        msg_out = pmt.dict_add(msg_out, pmt.intern('symb'), pmt.to_pmt(symb['symb']))
        msg_out = pmt.dict_add(msg_out, pmt.intern('scramble'), pmt.to_pmt(symb['scramble']))
        msg_out = pmt.dict_add(msg_out, pmt.intern('constellation_idx'), pmt.to_pmt(constellation_idx))
        msg_out = pmt.dict_add(msg_out, pmt.intern('do_continue'), pmt.to_pmt(np.bool(do_continue)))
        msg_out = pmt.dict_add(msg_out, pmt.intern('save_soft_dec'), pmt.to_pmt(np.bool(save_soft_dec)))
        ## print(msg_out)
        self.message_port_pub(self._port_frame_info, msg_out)
