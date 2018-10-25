## -*- python -*-

import numpy as np
from gnuradio import digital

class PhysicalLayer(object):
    """Physical layer description for STANAG 4285"""

    def __init__(self, mode=0):
        """For STANAG 4258 the mode has to be set manually: mode=0 -> BPSK, mode=1 -> QPSK, mode=2 -> 8PSK"""
        self._constellations = [PhysicalLayer.make_psk(2, [0,1]),
                                PhysicalLayer.make_psk(4, [0,1,3,2]),
                                PhysicalLayer.make_psk(8, [1,0,2,3,6,7,5,4])]
        self._preamble = [PhysicalLayer.get_preamble(), 0] ## BPSK
        self._data     = [PhysicalLayer.get_data(),  mode] ## according to the mode
        self._counter  = 0

    def set_mode(self):
        """For STANAG 4258 the mode has to be set manually: mode=0 -> BPSK, mode=1 -> QPSK, mode=2 -> 8PSK"""
        self._data[1] = mode

    def get_constellations(self):
        return self._constellations

    def get_frame(self):
        """returns the known+unknown symbols and scrambling"""
        if self._counter == 0:
            return self._preamble
        else:
            return self._data

    def get_doppler(self): ## symbols
        """used for doppler shift update, for determining which frame to provide next,
        and for stopping at end of data/when the signal quality is too low"""
        self._counter = (self._counter+1)&1
        ## TODO: doppler calculations
        doppler = 0.1234
        return [True, doppler]

    @staticmethod
    def get_preamble():
        """preamble symbols + scrambler(=1)"""
        state = np.array([1,1,0,1,0], dtype=np.bool)
        taps  = np.array([0,0,1,0,1], dtype=np.bool)
        p = np.zeros(80, dtype=np.uint8)
        for i in range(80):
            p[i]      = state[-1]
            state     = np.concatenate(([np.sum(state&taps)&1], state[0:-1]))
        a = np.zeros(80, dtype=[('symb',np.complex64), ('scramble', np.complex64)])
        ## BPSK modulation
        constellation = PhysicalLayer.make_psk(2,range(2))['points']
        a['symb']     = constellation[p,]
        a['scramble'] = 1
        return a

    @staticmethod
    def get_data():
        """data symbols + scrambler; for unknown symbols 'symb'=0"""
        state = np.array([1,1,1,1,1,1,1,1,1], dtype=np.bool)
        taps =  np.array([0,0,0,0,1,0,0,0,1], dtype=np.bool)
        p = np.zeros(176, dtype=np.uint8)
        for i in range(176):
            p[i] = np.sum(state[-3:]*[4,2,1])
            for j in range(3):
                state = np.concatenate(([np.sum(state&taps)&1], state[0:-1]))
        a=np.zeros(176, dtype=[('symb',np.complex64), ('scramble', np.complex64)])
        ## PSK-8 modulation
        constellation = PhysicalLayer.make_psk(8,range(8))['points']
        a['scramble'] = constellation[p,]
        a['symb'][ 32: 48] = 1 ## mini-probe 1
        a['symb'][ 80: 96] = 1 ## mini-probe 2
        a['symb'][128:144] = 1 ## mini-probe 3
        return a

    @staticmethod
    def make_psk(n, gray_code):
        c = np.zeros(n, dtype=[('points', np.complex64), ('symbols', np.uint8)])
        c['points']  = np.exp(2*np.pi*1j*np.array(range(n))/n)
        c['symbols'] = gray_code
        return c
