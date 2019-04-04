#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2018 hcab14@mail.com.
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

import digitalhf.digitalhf_swig as digitalhf
import random
import pmt
def main():
    v = digitalhf.viterbi48(0xB9, 0x9D, 0xD3, 0xF7)
    N = 40*500+5
    bits = [random.randint(0,1) for _ in range(N+8)]
    a = [1,0,0,1,1,1,0,1]
    b = [1,0,1,1,1,0,0,1]
    c = [1,1,0,0,1,0,1,1]
    d = [1,1,1,0,1,1,1,1]
    llr_encoded_bits = []
    for i in range(N):
        t1=t2=t3=t4=0
        for j in range(8):
            t1 += a[j]*bits[i+8-j]
            t2 += b[j]*bits[i+8-j]
            t3 += c[j]*bits[i+8-j]
            t4 += d[j]*bits[i+8-j]
        llr_encoded_bits.extend([7*(1-2*(t1%2)), 7*(1-2*(t2%2)), 7*(1-2*(t3%2)), 7*(1-2*(t4%2))])
    v.reset()
    decoded_bits = v.udpate(llr_encoded_bits)
    print(bits[8:38])
    print(decoded_bits[0:30])
    print('quality:', v.quality())
    test = [all([decoded_bits[i] == bits[i+8] for i in range(N)]), abs(v.quality()-4*N)<1]
    print('test:', test)
    if not all(test):
        raise Exception(test)

if __name__ == '__main__':
    main()
