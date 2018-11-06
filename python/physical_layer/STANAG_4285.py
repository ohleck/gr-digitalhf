## -*- python -*-

import numpy as np

class PhysicalLayer(object):
    """Physical layer description for STANAG 4285"""

    MODE_BPSK=0
    MODE_QPSK=1
    MODE_8PSK=2

    def __init__(self, sps):
        """intialization"""
        self._sps     = sps
        self._mode    = self.MODE_BPSK
        self._frame_counter = 0
        self._is_first_frame = True
        self._constellations = [self.make_psk(2, [0,1]),
                                self.make_psk(4, [0,1,3,2]),
                                self.make_psk(8, [1,0,2,3,6,7,5,4])]
        self._preamble = self.get_preamble()
        self._data     = self.get_data()

    def set_mode(self, mode):
        """set phase modultation type"""
        print('set_mode', mode)
        self._mode = int(mode)

    def get_constellations(self):
        return self._constellations

    def get_frame(self):
        """returns a tuple describing the frame:
        [0] ... known+unknown symbols and scrambling
        [1] ... modulation type after descrambling
        [2] ... a boolean indicating whethere or not raw IQ samples needed
        [3] ... a boolean indicating if the soft decision for the unknown symbols are saved"""
        print('-------------------- get_frame --------------------', self._frame_counter)
        return [self._preamble,self.MODE_BPSK,True,False] if self.is_preamble() else [self._data,self._mode,False,True]

    def get_doppler(self, symbols, iq_samples):
        """returns a tuple
        [0] ... quality flag
        [1] ... doppler estimate (rad/symbol) if available"""
        print('-------------------- get_doppler --------------------', self._frame_counter,len(symbols),len(iq_samples))
        success,doppler = self.quality_preamble(symbols,iq_samples) if self.is_preamble() else self.quality_data(symbols)
        if len(symbols) != 0:
            self._frame_counter = (self._frame_counter+1)&1 if success else 0
            self._is_first_frame = not success
        return success,doppler

    def is_preamble(self):
        return self._frame_counter == 0

    def quality_preamble(self, symbols, iq_samples):
        """quality check and doppler estimation for preamble"""
        success = True
        doppler = 0
        if len(iq_samples) != 0:
            sps  = self._sps
            zp   = np.array([x for x in self._preamble['symb'][9:40]
                             for _ in range(sps)], dtype=np.complex64)
            cc   = np.correlate(iq_samples, zp)
            imax = np.argmax(np.abs(cc[0:18*sps]))
            pks  = cc[(imax,imax+31*sps),]
            tpks = cc[imax+15*sps:imax+16*sps]
            success = np.mean(np.abs(pks)) > 2*np.mean(np.abs(tpks))
            doppler = np.diff(np.unwrap(np.angle(pks)))[0]/31 if success else 0
        if len(symbols) != 0:
            idx = range(30,80) if self._is_first_frame else range(80)
            z = symbols[idx]*np.conj(self._preamble['symb'][idx])
            print('quality_preamble',np.sum(np.real(z)<0))
            success = np.sum(np.real(z)<0) < 30
        return success,doppler

    def quality_data(self, s):
        """quality check for the data frame"""
        known_symbols = np.mod(range(176),48)>=32
        print('quality_data',np.sum(np.real(s[known_symbols])<0))
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
            for _ in range(3):
                state = np.concatenate(([np.sum(state&taps)&1], state[0:-1]))
        a=np.zeros(176, dtype=[('symb',np.complex64), ('scramble', np.complex64)])
        ## 8PSK modulation
        constellation = PhysicalLayer.make_psk(8,range(8))['points']
        a['scramble'] = constellation[p,]
        known_symbols = np.mod(range(176),48)>=32
        a['symb'][known_symbols] = a['scramble'][known_symbols]
        return a

    @staticmethod
    def make_psk(n, gray_code):
        """generates n-PSK constellation data"""
        c = np.zeros(n, dtype=[('points', np.complex64), ('symbols', np.uint8)])
        c['points']  = np.exp(2*np.pi*1j*np.arange(n)/n)
        c['symbols'] = gray_code
        return c
