## -*- python -*-

from __future__ import print_function
import numpy as np
from common import *

## ---- constellations -----------------------------------------------------------
BPSK=np.array(zip(np.exp(2j*np.pi*np.arange(2)/2), [0,1]), CONST_DTYPE)
QPSK=np.array(zip(np.exp(2j*np.pi*np.arange(4)/4), [0,1,3,2]), CONST_DTYPE)
PSK8=np.array(zip(np.exp(2j*np.pi*np.arange(8)/8), [0,1,3,2,7,6,4,5]), CONST_DTYPE)
QAM16=np.array(
    zip([+0.866025+0.500000j,  0.500000+0.866025j,  1.000000+0.000000j,  0.258819+0.258819j,
         -0.500000+0.866025j,  0.000000+1.000000j, -0.866025+0.500000j, -0.258819+0.258819j,
         +0.500000-0.866025j,  0.000000-1.000000j,  0.866025-0.500000j,  0.258819-0.258819j,
         -0.866025-0.500000j, -0.500000-0.866025j, -1.000000+0.000000j, -0.258819-0.258819j],
        range(16)), CONST_DTYPE)
QAM32=np.array(
    zip([+0.866380+0.499386j,  0.984849+0.173415j,  0.499386+0.866380j,  0.173415+0.984849j,
         +0.520246+0.520246j,  0.520246+0.173415j,  0.173415+0.520246j,  0.173415+0.173415j,
         -0.866380+0.499386j, -0.984849+0.173415j, -0.499386+0.866380j, -0.173415+0.984849j,
         -0.520246+0.520246j, -0.520246+0.173415j, -0.173415+0.520246j, -0.173415+0.173415j,
         +0.866380-0.499386j,  0.984849-0.173415j,  0.499386-0.866380j,  0.173415-0.984849j,
         +0.520246-0.520246j,  0.520246-0.173415j,  0.173415-0.520246j,  0.173415-0.173415j,
         -0.866380-0.499386j, -0.984849-0.173415j, -0.499386-0.866380j, -0.173415-0.984849j,
         -0.520246-0.520246j, -0.520246-0.173415j, -0.173415-0.520246j, -0.173415-0.173415j],
        range(32)), CONST_DTYPE)
QAM64=np.array(
    zip([+1.000000+0.000000j,  0.822878+0.568218j,  0.821137+0.152996j,  0.932897+0.360142j,
         +0.000000-1.000000j,  0.822878-0.568218j,  0.821137-0.152996j,  0.932897-0.360142j,
         +0.568218+0.822878j,  0.588429+0.588429j,  0.588429+0.117686j,  0.588429+0.353057j,
         +0.568218-0.822878j,  0.588429-0.588429j,  0.588429-0.117686j,  0.588429-0.353057j,
         +0.152996+0.821137j,  0.117686+0.588429j,  0.117686+0.117686j,  0.117686+0.353057j,
         +0.152996-0.821137j,  0.117686-0.588429j,  0.117686-0.117686j,  0.117686-0.353057j,
         +0.360142+0.932897j,  0.353057+0.588429j,  0.353057+0.117686j,  0.353057+0.353057j,
         +0.360142-0.932897j,  0.353057-0.588429j,  0.353057-0.117686j,  0.353057-0.353057j,
         +0.000000+1.000000j, -0.822878+0.568218j, -0.821137+0.152996j, -0.932897+0.360142j,
         -1.000000+0.000000j, -0.822878-0.568218j, -0.821137-0.152996j, -0.932897-0.360142j,
         -0.568218+0.822878j, -0.588429+0.588429j, -0.588429+0.117686j, -0.588429+0.353057j,
         -0.568218-0.822878j, -0.588429-0.588429j, -0.588429-0.117686j, -0.588429-0.353057j,
         -0.152996+0.821137j, -0.117686+0.588429j, -0.117686+0.117686j, -0.117686+0.353057j,
         -0.152996-0.821137j, -0.117686-0.588429j, -0.117686-0.117686j, -0.117686-0.353057j,
         -0.360142+0.932897j, -0.353057+0.588429j, -0.353057+0.117686j, -0.353057+0.353057j,
         -0.360142-0.932897j, -0.353057-0.588429j, -0.353057-0.117686j, -0.353057-0.353057j],
        range(64)), CONST_DTYPE)

## for test
#QAM64 = QAM64[(7,3,24,56,35,39,60,28),]
#QAM64['symbols'] = [1, 0, 2, 6, 4, 5, 7, 3]

## ---- constellation indices ---------------------------------------------------
MODE_BPSK  = 0
MODE_QPSK  = 1
MODE_8PSK  = 2
MODE_16QAM = 3
MODE_32QAM = 4
MODE_64QAM = 5

## ---- data scrambler -----------------------------------------------------------
class ScrambleData(object):
    """data scrambling sequence generator"""
    def __init__(self):
        self.reset()

    def reset(self):
        self._state = np.array([0,0,0,0,0,0,0,0,1], dtype=np.bool)
        self._taps =  np.array([0,0,0,0,1,0,0,0,1], dtype=np.bool)

    def next(self, num_bits):
        r = np.packbits(self._state[1:])[0]&((1<<num_bits)-1)
        for _ in range(num_bits):
            self._advance()
        return r

    def _advance(self):
        self._state = np.concatenate(([np.sum(self._state&self._taps)&1],
                                      self._state[0:-1]))

## ---- preamble definitions  ---------------------------------------------------
## 184 = 8*23
PREAMBLE=n_psk(8, np.array(
    [1,5,1,3,6,1,3,1,1,6,3,7,7,3,5,4,3,6,6,4,5,4,0,
     2,2,2,6,0,7,5,7,4,0,7,5,7,1,6,1,0,5,2,2,6,2,3,
     6,0,0,5,1,4,2,2,2,3,4,0,6,2,7,4,3,3,7,2,0,2,6,
     4,4,1,7,6,2,0,6,2,3,6,7,4,3,6,1,3,7,4,6,5,7,2,
     0,1,1,1,4,4,0,0,5,7,7,4,7,3,5,4,1,6,5,6,6,4,6,
     3,4,3,0,7,1,3,4,7,0,1,4,3,3,3,5,1,1,1,4,6,1,0,
     6,0,1,3,1,4,1,7,7,6,3,0,0,7,2,7,2,0,2,6,1,1,1,
     2,7,7,5,3,3,6,0,5,3,3,1,0,7,1,1,0,3,0,4,0,7,3]))

## 103 = 31 + 1 + 3*13 + 1 + 31
REINSERTED_PREAMBLE=n_psk(8, np.array(
    [0,0,0,0,0,2,4,6,0,4,0,4,0,6,4,2,0,0,0,0,0,2,4,6,0,4,0,4,0,6,4,   ## MP+
     2,
     0,4,0,4,0,0,4,4,0,0,0,0,0, # + D0
     0,4,0,4,0,0,4,4,0,0,0,0,0, # + D1
     0,4,0,4,0,0,4,4,0,0,0,0,0, # + D2
     6,
     4,4,4,4,4,6,0,2,4,0,4,0,4,2,0,6,4,4,4,4,4,6,0,2,4,0,4,0,4,2,0])) ## MP-
## length 31
MINI_PROBE=[n_psk(8, np.array([0,0,0,0,0,2,4,6,0,4,0,4,0,6,4,2,0,0,0,0,0,2,4,6,0,4,0,4,0,6,4])), ## sign = + (0)
            n_psk(8, np.array([4,4,4,4,4,6,0,2,4,0,4,0,4,2,0,6,4,4,4,4,4,6,0,2,4,0,4,0,4,2,0]))] ## sign = - (1)

## ---- di-bits ----------------------------------------------------------------
TO_DIBIT=[(0,0),(0,1),(1,1),(1,0)]

## ---- rate -------------------------------------------------------------------
TO_RATE={(0,0,0): {'baud':    0, 'bits_per_symbol': 0},  ## reserved
         (0,0,1): {'baud': 3200, 'bits_per_symbol': 2, 'ci': MODE_QPSK},
         (0,1,0): {'baud': 4800, 'bits_per_symbol': 3, 'ci': MODE_8PSK},
         (0,1,1): {'baud': 6400, 'bits_per_symbol': 4, 'ci': MODE_16QAM},
         (1,0,0): {'baud': 8000, 'bits_per_symbol': 5, 'ci': MODE_32QAM},
         (1,0,1): {'baud': 9600, 'bits_per_symbol': 6, 'ci': MODE_64QAM},
         (1,1,0): {'baud':12800, 'bits_per_symbol': 6, 'ci': MODE_64QAM},
         (1,1,1): {'baud':    0, 'bits_per_symbol': 0}}  ## reserved

## ---- interleaver ------------------------------------------------------------
TO_INTERLEAVER={(0,0,0): {'frames': -1, 'name': 'illegal'},
                (0,0,1): {'frames':  1, 'name': 'Ultra Short (US)'},
                (0,1,0): {'frames':  3, 'name': 'Very Short (VS)'},
                (0,1,1): {'frames':  9, 'name': 'Short (S)'},
                (1,0,0): {'frames': 18, 'name': 'Medium (M)'},
                (1,0,1): {'frames': 36, 'name': 'Long (L)'},
                (1,1,0): {'frames': 72, 'name': 'very Long (VL)'},
                (1,1,1): {'frames': -1, 'name': 'illegal'}}

MP_COUNTER=[(0,0,1),
            (0,1,0),
            (0,1,1),
            (1,0,0)]

## ---- physcal layer class -----------------------------------------------------
class PhysicalLayer(object):
    """Physical layer description for MIL-STD-188-110 Appendix C = STANAG 4539"""

    def __init__(self, sps):
        """intialization"""
        self._sps = sps
        self._frame_counter = -1
        self._constellations = [BPSK, QPSK, PSK8, QAM16, QAM32, QAM64]
        self._preamble = self.get_preamble()

    def get_constellations(self):
        return self._constellations

    def get_next_frame(self, symbols):
        """returns a tuple describing the frame:
        [0] ... known+unknown symbols and scrambling
        [1] ... modulation type after descrambling
        [2] ... a boolean indicating if the processing should continue
        [3] ... a boolean indicating if the soft decision for the unknown
                symbols are saved"""
        print('-------------------- get_frame --------------------', self._frame_counter)
        success = True
        if self._frame_counter == -1: ## ---- preamble
            self._preamble_offset = 0
            self._frame_counter += 1
            return [self._preamble,MODE_BPSK,success,False]

        frame_counter_mod72 = self._frame_counter%72
        if frame_counter_mod72 == 0: ## --- re-inserted preamble
            self._frame_counter += 1
            success = self.get_preamble_quality(symbols)
            return [self.make_reinserted_preamble(self._preamble_offset,success),MODE_QPSK,success,False]

        if frame_counter_mod72 >= 1: ## ---- data frames
            got_reinserted_preamble = frame_counter_mod72 == 1
            self._frame_counter += 1
            if got_reinserted_preamble:
                success = self.decode_reinserted_preamble(symbols)
            else:
                success = self.get_data_frame_quality(symbols)
            return [self.make_data_frame(success),self._constellation_index,success,not got_reinserted_preamble]

    def get_doppler(self, iq_samples):
        """quality check and doppler estimation for preamble"""
        success,doppler = True,0
        if len(iq_samples) != 0:
            sps  = self._sps
            m    = 23*sps
            idx  = np.arange(m)
            idx2  = np.arange(m+23*sps)
            _,zp = self.get_preamble_z()
            n    = len(zp)
            cc   = np.correlate(iq_samples, zp)
            imax = np.argmax(np.abs(cc[0:23*sps]))
            print('imax=', imax, len(iq_samples))
            pks  = [np.correlate(iq_samples[imax+i*m+idx],
                                 zp[i*m+idx])[0]
                    for i in range(n//m)]
            val  = [np.mean(np.abs(np.correlate(iq_samples[imax+i*m+idx2],
                                                zp[i*m+idx])[11*sps+np.arange(-2*sps,2*sps)]))
                    for i in range((n//m)-1)]
            tests = np.abs(pks[0:-1])/val
            success = np.median(tests) > 2.0
            print('test:', np.abs(pks), tests)
            if success:
                print('doppler apks', np.abs(pks))
                print('doppler ppks', np.angle(pks),
                      np.diff(np.unwrap(np.angle(pks)))/m,
                      np.mean(np.diff(np.unwrap(np.angle(pks)))/m))
                doppler = freq_est(pks)/m;
            print('success=', success, 'doppler=', doppler)
        return success,doppler

    def set_mode(self, mode):
        pass

    def get_preamble_quality(self, symbols):
        return np.abs(np.mean(symbols[-40:])) > 0.5

    def get_data_frame_quality(self, symbols):
        return np.abs(np.mean(symbols[-31:])) > 0.5

    def decode_reinserted_preamble(self, symbols):
        ## decode D0,D1,D2
        idx    = np.arange(13)
        z      = np.array([np.mean(symbols[32+13*i+idx]) for i in range(3)])
        d0d1d2 = map(np.uint8, np.mod(np.round(np.angle(z)/np.pi*2),4))
        dibits = [TO_DIBIT[idx] for idx in d0d1d2]
        self._mode = {'rate':        tuple([x[0] for x in dibits]),
                      'interleaver': tuple([x[1] for x in dibits])}
        print('======== rate,interleaver:',
              TO_RATE[self._mode['rate']],
              TO_INTERLEAVER[self._mode['interleaver']])
        rate_info = TO_RATE[self._mode['rate']]
        print('rate_info', rate_info)
        self._constellation_index = rate_info['ci']
        print('constellation index', self._constellation_index)
        scr = ScrambleData()
        iscr = [scr.next(rate_info['bits_per_symbol']) for _ in range(256)]
        if rate_info['ci'] > MODE_8PSK:
            self._data_scramble = np.ones(256, dtype=np.complex64)
        else:
            constell = self._constellations[rate_info['ci']]
            self._data_scramble = constell[iscr]['points']
        success = True ## TODO
        return success

    def make_reinserted_preamble(self, offset, success):
        """ offset=  0 -> 1st reinserted preamble
            offset=-72 -> all following reinserted preambles"""
        print('make_reinserted_preamble', offset, success)
        a=np.array(zip(REINSERTED_PREAMBLE[offset:],
                       REINSERTED_PREAMBLE[offset:]),
                   dtype=[('symb',     np.complex64),
                          ('scramble', np.complex64)])
        a['symb'][-72:-72+3*13] = 0 ## D0,D1,D2
        if not success:
            self._frame_counter = -1
        return a

    def make_data_frame(self, success):
        self._preamble_offset = -72 ## all following reinserted preambles start at index -72
        a = np.zeros(256+31, dtype=[('symb',     np.complex64),
                                    ('scramble', np.complex64)])
        a['scramble'][:256] = self._data_scramble
        n = (self._frame_counter-2)%72
        m = n%18
        if m == 0:
            cnt = n//18
            self._mp = (1,1,1,1,1,1,1,0)+self._mode['rate']+self._mode['interleaver']+MP_COUNTER[cnt]+(0,)
            print('new mini-probe signs n=',n,'m=',m,self._mp)
        a['symb'][256:]     = MINI_PROBE[self._mp[m]]
        a['scramble'][256:] = MINI_PROBE[self._mp[m]]
        if not success:
            self._frame_counter = -1
        return a

    def decode_soft_dec(self, soft_dec):
        return soft_dec

    @staticmethod
    def get_preamble():
        """preamble symbols + scrambler"""
        return np.array(zip(PREAMBLE,
                            PREAMBLE),
                        dtype=[('symb',     np.complex64),
                               ('scramble', np.complex64)])

    def get_preamble_z(self):
        """preamble symbols for preamble correlation"""
        return 2,np.array([z for z in PREAMBLE for _ in range(self._sps)])

if __name__ == '__main__':
    print(PREAMBLE)
    z = n_psk(8,PREAMBLE)
    cc = [np.sum(z[0:23]*np.conj(z[23*i:23*i+23])) for i in range(6)]
    print(np.abs(cc))
    print(np.angle(cc)/np.pi*4)
    print(all(z==PhysicalLayer.get_preamble()['symb']))
    print(len(PhysicalLayer.get_preamble()['symb']))
    s = ScrambleData()
    print([s.next(1) for _ in range(511)])
    print([s.next(1) for _ in range(511)] ==
          [s.next(1) for _ in range(511)])
    #print(QAM64)
    #print(QAM32)
    #print(QAM16)
    #print(PSK8)
    #print(QPSK)
    #print(BPSK)
    #print(MINI_PROBE_PLUS)
    #print(MINI_PROBE_MINUS)
    #print(MINI_PROBE_PLUS*MINI_PROBE_MINUS)
    #for i in range(len(QAM64)):
    #    print(QAM64['points'][i])

    print([s.next(6) for _ in range(256)])
