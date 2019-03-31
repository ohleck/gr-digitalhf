## -*- python -*-

import numpy as np

class Deinterleaver(object):
    "S4285 deinterleaver"
    def __init__(self, incr):
        ## incr = 12 -> L
        ## incr =  1 -> S
        self._buf = [np.zeros(incr*(31-i) + 1) for i in range(32)]

    def push(self, a):
        assert(len(a) == 32)
        for i in range(32):
            self._buf[i][0] = a[i]
            self._buf[i] = np.roll(self._buf[i],1)

    def fetch(self):
        return np.array([self._buf[(9*i)%32][0] for i in range(32)])

class PhysicalLayer(object):
    """Physical layer description for STANAG 4285"""

    MODE_BPSK=0
    MODE_QPSK=1
    MODE_8PSK=2

    def __init__(self, sps):
        """intialization"""
        self._sps     = sps
        self._mode    = self.MODE_QPSK
        self._frame_counter = 0
        self._is_first_frame = True
        self._constellations = [self.make_psk(2, [0,1]),
                                self.make_psk(4, [0,1,3,2]),
                                self.make_psk(8, [1,0,2,3,6,7,5,4])]
        self._preamble = self.get_preamble()
        self._data     = self.get_data()
        self._deinterleaver = Deinterleaver(12) ## for now BPSK L fixed

    def set_mode(self, mode):
        """set phase modultation type"""
        print('set_mode', mode)
        self._mode = int(mode)

    def get_constellations(self):
        return self._constellations

    def get_next_frame(self, symbols):
        """returns a tuple describing the frame:
        [0] ... known+unknown symbols and scrambling
        [1] ... modulation type after descrambling
        [2] ... a boolean indicating if the processing should continue
        [3] ... a boolean indicating if the soft decision for the unknown symbols are saved"""
        ## print('-------------------- get_frame --------------------', self._frame_counter, len(symbols))
        if len(symbols) == 0: ## 1st preamble
            self._frame_counter = 0

        success,frame_description = True,[]
        if (self._frame_counter%2) == 0:
            frame_description = [self._preamble,self.MODE_BPSK,success,False]
        else:
            idx = range(30,80)
            z = symbols[idx]*np.conj(self._preamble['symb'][idx])
            ## print('quality_preamble',np.sum(np.real(z)<0), symbols[idx])
            success = np.sum(np.real(z)<0) < 30
            frame_description = [self._data,self._mode,success,True]

        self._frame_counter += 1
        return frame_description

    def get_doppler(self, iq_samples):
        """returns a tuple
        [0] ... quality flag
        [1] ... doppler estimate (rad/symbol) if available"""
        ## print('-------------------- get_doppler --------------------', self._frame_counter,len(iq_samples))
        success,doppler = False,0
        if len(iq_samples) == 0:
            return success,doppler

        sps  = self._sps
        zp   = np.array([x for x in self._preamble['symb'][9:40]
                         for _ in range(sps)], dtype=np.complex64)
        cc   = np.correlate(iq_samples, zp)
        imax = np.argmax(np.abs(cc[0:18*sps]))
        pks  = cc[(imax,imax+31*sps),]
        tpks = cc[imax+15*sps:imax+16*sps]
        ## print('doppler: ', np.abs(pks), np.abs(tpks))
        success = np.mean(np.abs(pks)) > 5*np.mean(np.abs(tpks))
        doppler = np.diff(np.unwrap(np.angle(pks)))[0]/31/self._sps if success else 0
        return success,doppler

    def is_preamble(self):
        return self._frame_counter == 0

    def quality_data(self, s):
        """quality check for the data frame"""
        known_symbols = np.mod(range(176),48)>=32
        print('quality_data',np.sum(np.real(s[known_symbols])<0))
        success = np.sum(np.real(s[known_symbols])<0) < 20
        return success,0 ## no doppler estimate for data frames

    def get_preamble_z(self):
        """preamble symbols for preamble correlation"""
        a = PhysicalLayer.get_preamble()
        return 2,np.array([z for z in a['symb'][0:31] for _ in range(self._sps)])

    def decode_soft_dec(self, soft_dec):
        assert(len(soft_dec) == 128)
        print('decode_soft_dec: ', len(soft_dec))
        res = []
        for i in range(0,128,32):
            self._deinterleaver.push(soft_dec[i:i+32])
            res.extend(self._deinterleaver.fetch().tolist())
        return res

    @staticmethod
    def get_preamble():
        """preamble symbols + scrambler(=1)"""
        state = np.array([1,1,0,1,0], dtype=np.bool)
        taps  = np.array([0,0,1,0,1], dtype=np.bool)
        p = np.zeros(80, dtype=np.uint8)
        for i in range(80):
            p[i]      = state[-1]
            state     = np.concatenate(([np.sum(state&taps)&1], state[0:-1]))
        a = np.zeros(80, common.SYMB_SCRAMBLE_DTYPE)
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
        a = np.zeros(176, common.SYMB_SCRAMBLE_DTYPE)
        ## 8PSK modulation
        constellation = PhysicalLayer.make_psk(8,range(8))['points']
        a['scramble'] = constellation[p,]
        known_symbols = np.mod(range(176),48)>=32
        a['symb'][known_symbols] = a['scramble'][known_symbols]
        return a

    @staticmethod
    def make_psk(n, gray_code):
        """generates n-PSK constellation data"""
        c = np.zeros(n, dtype=[('points', np.complex64), ('symbols', np.int32)])
        c['points']  = np.exp(2*np.pi*1j*np.arange(n)/n)
        c['symbols'] = gray_code
        return c
