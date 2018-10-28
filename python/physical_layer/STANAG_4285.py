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
        self._is_first_frame  = True

    def set_mode(self, mode):
        """For STANAG 4258 the mode has to be set manually: mode=0 -> BPSK, mode=1 -> QPSK, mode=2 -> 8PSK"""
        print('set_mode', mode)
        self._data[1] = int(mode)

    def get_constellations(self):
        return self._constellations

    def get_frame(self):
        """returns the known+unknown symbols and scrambling"""
        print('-------------------- get_frame --------------------',self._counter)
        return self._preamble if self._counter == 0 else self._data

    def get_doppler(self, sy, sa):
        """used for doppler shift update, for determining which frame to provide next,
        and for stopping at end of data/when the signal quality is too low
        sy ... equalized symbols; sa ... samples"""
        print('-------------------- get_doppler --------------------',self._counter,len(sy),len(sa))
        success,doppler = self.quality_preamble(sy,sa) if self._counter == 0 else self.quality_data(sy)
        self._counter = (self._counter+1)&1 if success else 0
        self._is_first_frame = not success
        return success,doppler

    def quality_preamble(self, sy, sa):
        sps  = 5
        zp   = [x for x in self._preamble[0]['symb'][9:40] for i in range(sps)]
        cc   = np.array([np.sum(sa[ i*5:(31+i)*5]*zp) for i in range(49)])
        imax = np.argmax(np.abs(cc[0:18]))
        pks  = cc[(imax,imax+15,imax+16,imax+31),]
        apks = np.abs(pks)
        test = np.mean(apks[(0,3),]) > 2*np.mean(apks[(1,2),])
        doppler = np.diff(np.unwrap(np.angle(pks[(0,3),])))[0]/31 if test else 0
        idx = range(80)
        if self._is_first_frame:
           idx = range(30,80)
        z = sy[idx]*np.conj(self._preamble[0]['symb'][idx])
        success = np.sum(np.real(z)<0) < 30
        return success,doppler

    def quality_data(self, s):
        known_symbols = np.mod(range(176),48)>=32
        success = np.sum(np.real(s[known_symbols])<0) < 20
        return success,0 ## no doppler estimate for data frames

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
        known_symbols = np.mod(range(176),48)>=32
        a['symb'][known_symbols] = a['scramble'][known_symbols]
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
