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
        self._preamble_phases = []

    def set_mode(self, mode):
        """For STANAG 4258 the mode has to be set manually: mode=0 -> BPSK, mode=1 -> QPSK, mode=2 -> 8PSK"""
        self._data[1] = mode

    def get_constellations(self):
        return self._constellations

    def get_frame(self):
        """returns the known+unknown symbols and scrambling"""
        print('-------------------- get_frame --------------------',self._counter)
        if self._counter == 0:
            x= self._preamble
        else:
            x=self._data
        print('get_frame end\n')
        return x;

    def get_doppler(self, s):
        """used for doppler shift update, for determining which frame to provide next,
        and for stopping at end of data/when the signal quality is too low"""
        print('-------------------- get_doppler --------------------',self._counter)
        doppler = 0
        if self._counter == 0: ## preamble
            doppler = PhysicalLayer.data_aided_frequency_estimation(s, self._preamble[0]['symb'])
        self._counter = (self._counter+1)&1
        return [True, 2*doppler]

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
        a['symb'][ 32: 48] = a['scramble'][ 32: 48] ## mini-probe 1
        a['symb'][ 80: 96] = a['scramble'][ 80: 96] ## mini-probe 2
        a['symb'][128:144] = a['scramble'][128:144] ## mini-probe 3
        return a

    @staticmethod
    def make_psk(n, gray_code):
        c = np.zeros(n, dtype=[('points', np.complex64), ('symbols', np.uint8)])
        c['points']  = np.exp(2*np.pi*1j*np.array(range(n))/n)
        c['symbols'] = gray_code
        return c

    @staticmethod
    def data_aided_frequency_estimation(x,c):
        """Data-Aided Frequency Estimation for Burst Digital Transmission,
        Umberto Mengali and M. Morelli, IEEE TRANSACTIONS ON COMMUNICATIONS,
        VOL. 45, NO. 1, JANUARY 1997"""
        z  = x*np.conj(c)                                      ## eq (2)
        L0 = len(z)
        N  = L0//2
        R  = np.zeros(N, dtype=np.complex64)
        for i in range(N):
            R[i] = 1.0/(L0-i)*np.sum(z[i:]*np.conj(z[0:L0-i])) ## eq (3)
        m  = np.array(range(N), dtype=np.float)
        w  = 3*((L0-m)*(L0-m+1)-N*(L0-N))/(N*(4*N*N - 6*N*L0 + 3*L0*L0-1)) ## eq (9)
        mod_2pi = lambda x : np.mod(x-np.pi, 2*np.pi) - np.pi
        fd = np.sum(w[1:] * mod_2pi(np.diff(np.angle(R))))     ## eq (8)
        return fd
