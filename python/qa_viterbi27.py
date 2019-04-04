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
    v = digitalhf.viterbi27(0x6d, 0x4f)
    N = 35*500+5
    bits = [random.randint(0,1) for _ in range(N+7)]
    a = [1,0,1,1,0,1,1]
    b = [1,1,1,1,0,0,1]
    llr_encoded_bits = []
    for i in range(N):
        t1=t2=0
        for j in range(7):
            t1 += a[j]*bits[i+7-j]
            t2 += b[j]*bits[i+7-j]
        llr_encoded_bits.extend([7*(1-2*(t1%2)), 7*(1-2*(t2%2))])
    v.reset()
    decoded_bits = v.udpate(llr_encoded_bits)
    print(bits[7:37])
    print(decoded_bits[0:30])
    print('quality:', v.quality())
    test = [all([decoded_bits[i] == bits[i+7] for i in range(N)]), abs(v.quality()-2*N)<1]
    print('test:', test)
    if not all(test):
        raise Exception(test)

if __name__ == '__main__':
    main()
