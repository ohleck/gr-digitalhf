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
QAM64 = QAM64[(7,3,24,56,35,39,60,28),]
QAM64['symbols'] = [1, 0, 2, 6, 4, 5, 7, 3]

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
MINI_PROBE_PLUS=n_psk(8, np.array(
    [0,0,0,0,0,2,4,6,0,4,0,4,0,6,4,2,0,0,0,0,0,2,4,6,0,4,0,4,0,6,4]))
## length 31
MINI_PROBE_MINUS=n_psk(8, np.array(
    [4,4,4,4,4,6,0,2,4,0,4,0,4,2,0,6,4,4,4,4,4,6,0,2,4,0,4,0,4,2,0]))

## ---- physcal layer class -----------------------------------------------------
class PhysicalLayer(object):
    """Physical layer description for MIL-STD-188-110 Appendix C = STANAG 4539"""

    def __init__(self, sps):
        """intialization"""
        self._sps = sps
        self._frame_counter = -1
        self._constellations = [BPSK, QPSK, PSK8, QAM16, QAM32, QAM64]
        self._preamble = self.get_preamble()
        self._scr_data = ScrambleData()

    def get_constellations(self):
        return self._constellations

    def get_frame(self):
        """returns a tuple describing the frame:
        [0] ... known+unknown symbols and scrambling
        [1] ... modulation type after descrambling
        [2] ... a boolean indicating whethere or not raw IQ samples needed
        [3] ... a boolean indicating if the soft decision for the unknown
                symbols are saved"""
        print('-------------------- get_frame --------------------',
              self._frame_counter)
        ## --- preamble frame ----
        if self._frame_counter == -1:
            return [self._preamble,MODE_BPSK,True,False]
        ## ----- (re)inserted preamble ------
        if self._frame_counter == 0:
            self.a = self.make_reinserted_preamble()
            return [self.a, MODE_QPSK,False,False]
        if self._frame_counter >= 1:
            self.a = self.make_data_frame()
            return [self.a, MODE_64QAM,False,True]

    def get_doppler(self, symbols, iq_samples):
        """returns a tuple
        [0] ... quality flag
        [1] ... doppler estimate (rad/symbol) if available"""
        print('-------------------- get_doppler --------------------',
              self._frame_counter,len(symbols),len(iq_samples))
        #if len(symbols)!=0:
        #    print('symb=', symbols)
        success,doppler = False,0
        if self._frame_counter == -1: ## -- preamble ----
            success,doppler = self.get_doppler_from_preamble(symbols, iq_samples)
            if len(symbols) != 0:
                for s in symbols:
                    print(s)
                self._frame_counter = 0
        elif self._frame_counter >= 0 and self._frame_counter <5: ## -- reinserted preamble ----
            for s in symbols:
                print(s)
            self._frame_counter += 1
            success = True
        else: ## ------------------------ data frame ----
            for s in symbols:
                print(s)
            success = False
            self._frame_counter = -1
        return success,doppler

    def get_doppler_from_preamble(self, symbols, iq_samples):
        """quality check and doppler estimation for preamble"""
        success,doppler = True,0
        if len(iq_samples) != 0:
            sps  = self._sps
            idx  = np.arange(23*sps)
            zp   = self.get_preamble_z(self._sps)
            cc   = np.correlate(iq_samples, zp[idx])
            imax = np.argmax(np.abs(cc[0:23*sps]))
            pks  = [np.correlate(iq_samples[imax+i*23*sps+idx],
                                 zp[i*23*sps+idx])[0]
                    for i in range(7)]
            success = np.mean(np.abs(pks)) > 2*np.mean(np.abs(cc[imax+11*sps+range(-sps,sps)]))
            print('test:',imax, np.mean(np.abs(pks)), np.mean(np.abs(cc[imax+11*sps+range(-sps,sps)])))
            if success:
                print('doppler apks', np.abs(pks))
                print('doppler ppks', np.angle(pks),
                      np.diff(np.unwrap(np.angle(pks)))/23,
                      np.mean(np.diff(np.unwrap(np.angle(pks)))/23))
                doppler = freq_est(pks[1:])/23;
            print('success=', success, 'doppler=', doppler)
        return success,doppler

    def make_reinserted_preamble(self):
        a=np.array(zip(REINSERTED_PREAMBLE,
                       REINSERTED_PREAMBLE),
                   dtype=[('symb',     np.complex64),
                          ('scramble', np.complex64)])
        a['symb'][32:32+3*13] = 0 ## D0,D1,D2
        return a
    def make_data_frame(self):
        a=np.zeros(256+31,  dtype=[('symb',     np.complex64),
                                   ('scramble', np.complex64)])
        a['scramble'] = 1
        a['symb'][256:]     = MINI_PROBE_MINUS
        a['scramble'][256:] = MINI_PROBE_MINUS
        return a

    @staticmethod
    def get_preamble():
        """preamble symbols + scrambler"""
        return np.array(zip(PREAMBLE,
                            PREAMBLE),
                        dtype=[('symb',     np.complex64),
                               ('scramble', np.complex64)])
    @staticmethod
    def get_preamble_z(sps):
        """preamble symbols for preamble correlation"""
        return np.array([z for z in PREAMBLE for _ in range(sps)])

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
