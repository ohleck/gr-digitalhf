## -*- python -*-

from __future__ import print_function
import numpy as np

def n_psk(n,x):
    return np.complex64(np.exp(2j*np.pi*x/n))

## ---- constellations -----------------------------------------------------------
CONST_DTYPE=np.dtype([('points',  np.complex64),
                      ('symbols', np.uint8)])
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
        self._state = 1

    def next(self, num_bits):
        r = self._state & ((1<<num_bits)-1)
        for i in range(num_bits):
            self._advance()
        return r

    def _advance(self):
        lsb = self._state&1
        self._state = (self._state>>1)&511
        if lsb:
            self._state ^= 0x10B
        return self._state

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
        ## ----- data frame ------
        if self._frame_counter == 0:
            self.a = self.make_reinserted_preamble()
            return [self.a, MODE_QPSK,False,True]

    def get_doppler(self, symbols, iq_samples):
        """returns a tuple
        [0] ... quality flag
        [1] ... doppler estimate (rad/symbol) if available"""
        print('-------------------- get_doppler --------------------',
              self._frame_counter,len(symbols),len(iq_samples))
        #if len(symbols)!=0:
        #    print('symb=', symbols)
        success = False
        doppler = 0
        if self._frame_counter == -1: ## -- preamble ----
            success,doppler = self.get_doppler_from_preamble(symbols, iq_samples)
            if len(symbols) != 0:
                for s in symbols:
                    print(s)
                self._frame_counter = 0
        else: ## ------------------------ data frame ----
            if len(symbols) != 0:
                for s in symbols:
                    print(s)
            success = False
            self._frame_counter = -1
        return success,doppler

    def get_doppler_from_preamble(self, symbols, iq_samples):
        """quality check and doppler estimation for preamble"""
        success = True
        doppler = 0
        shift=9
        if len(iq_samples) != 0:
            zp   = np.conj(self.get_preamble_z(self._sps)[shift*self._sps:])
            cc   = np.array([np.sum(iq_samples[i:i+23*self._sps] *
                                    zp[0:23*self._sps])
                             for i in range(23*3*self._sps)])
            acc = np.abs(cc)
            for i in range(0,len(cc),23*self._sps):
                print('i=%3d: '%i,end='')
                for j in range(23*self._sps):
                    print('%3.0f ' % acc[i+j], end='')
                print()

            imax = np.argmax(np.abs(cc[0:3*23*self._sps]))
            print(imax)
            pks  = np.array([np.sum(iq_samples[(imax+23*i*self._sps):
                                               (imax+23*i*self._sps+23*self._sps)] *
                                    zp[(23*i*self._sps):
                                       (23*i*self._sps+23*self._sps)])
                             for i in range(1,5)])
            print('doppler apks', np.abs(pks))
            print('doppler ppks', np.angle(pks), np.diff(np.unwrap(np.angle(pks)))/23)
            doppler = np.mean(np.diff(np.unwrap(np.angle(pks))))/23
            success = True
            print('success=', success, 'doppler=', doppler)
        return success,doppler

    def make_reinserted_preamble(self):
        a=np.array(zip(REINSERTED_PREAMBLE,
                       REINSERTED_PREAMBLE),
                   dtype=[('symb',     np.complex64),
                          ('scramble', np.complex64)])
        a['symb'][32:32+3*13] = 0 ## D0,D1,D2
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
        return np.array([z for z in PREAMBLE for i in range(sps)])

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
    print(QAM64)
    print(QAM32)
    print(QAM16)
    print(PSK8)
    print(QPSK)
    print(BPSK)
    print(MINI_PROBE_PLUS)
    print(MINI_PROBE_MINUS)
    print(MINI_PROBE_PLUS*MINI_PROBE_MINUS)
