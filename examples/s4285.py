## -*- python -*-

import numpy as np
from gnuradio import digital

class PhysicalLayer(object):
    """Physical layer description for STANAG 4285"""

    def __init__(self, mode=0):
        """For STANAG 4258 the mode has to be set manually: mode=0 -> BPSK, mode=1 -> QPSK, mode=2 -> 8PSK"""
        print('Hello from PhysicalLayer!\n')
        self._preamble = [get_preamble(), digital.constellation_bpsk()]
        self._constellations = [digital.constellation_bpsk(),
                                digital.constellation_qpsk(),
                                digital.constellation_8psk()]
        self._data     = [get_data(), self._constellations[mode]]
        self._counter  = 0

    def __del__(self):
        print('Bye from PhysicalLayer!\n')

    def set_mode(self):
        """For STANAG 4258 the mode has to be set manually: mode=0 -> BPSK, mode=1 -> QPSK, mode=2 -> 8PSK"""
        self._data[1] = self._constellations[mode];

    def get_frame(self):
        """returns the known+unknown symbols and scrambling"""
        print('get_frame', self._counter, self._preamble, self._preamble[1].__deref__())
        if self._counter == 0:
            return self._preamble
        else:
            return self._data

    def get_doppler(self, symbols):
        """used for doppler shift update, for determining which frame to provide next,
        and for stopping at end of data/when the signal quality is too low"""
        self._counter = (self._counter+1)&2
        ## TODO: doppler calculations
        doppler = 0.0
        return [True, doppler]

def get_preamble():
    """preamble symbols + scrambler(=1)"""
    state = np.array([1,1,0,1,0], dtype='b')
    taps  = np.array([0,0,1,0,1], dtype='b')
    p = np.zeros(80, dtype='f')
    for i in range(80):
        p[i]      = state[-1];
        state     = np.concatenate(([np.sum(state&taps)&1], state[0:-1]))
    a = np.zeros(80, dtype=[('symb','c'), ('scramble', 'c')])
    ## BPSK modulation
    a['symb']     = np.exp(np.pi*1j*p)
    a['scramble'] = 1;
    return a

def get_data():
    """data symbols + scrambler; for unknown symbols 'symb'=0"""
    state = np.array([1,1,1,1,1,1,1,1,1], dtype='b')
    taps =  np.array([0,0,0,0,1,0,0,0,1], dtype='b')
    p = np.zeros(176, dtype='f')
    for i in range(176):
        p[i] = np.sum(state[-3:]*[4,2,1]);
        for j in range(3):
            state = np.concatenate(([np.sum(state&taps)&1], state[0:-1]))
    a=np.zeros(176, dtype=[('symb','c'), ('scramble', 'c')])
    ## PSK-8 modulation
    a['scramble'] = np.exp(np.pi*1j*p/4)
    a['symb'][ 32: 48] = 1; ## mini-probe 1
    a['symb'][ 80: 96] = 1; ## mini-probe 2
    a['symb'][128:144] = 1; ## mini-probe 3
    return a
