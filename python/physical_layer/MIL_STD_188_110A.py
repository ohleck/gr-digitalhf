## -*- python -*-

from __future__ import print_function
import numpy as np

## ---- Walsh-4 codes -----------------------------------------------------------
WALSH = np.array([[0,0,0,0, 0,0,0,0],
                  [0,1,0,1, 0,1,0,1],
                  [0,0,1,1, 0,0,1,1],
                  [0,1,1,0, 0,1,1,0],
                  [0,0,0,0, 1,1,1,1],
                  [0,1,0,1, 1,0,1,0],
                  [0,0,1,1, 1,1,0,0],
                  [0,1,1,0, 1,0,0,1]],
                 dtype=np.uint8)

def walsh_to_num(w):
    return sum(w*(1<<np.arange(8)[::-1]))

FROM_WALSH = -np.ones(256, dtype=np.int8)
for i in range(8):
    FROM_WALSH[walsh_to_num(WALSH[i][:])] = i

## ---- tri-bit codes -----------------------------------------------------------
TRIBIT = np.zeros((8,32), dtype=np.uint8)
for i in range(8):
    TRIBIT[i][:] = np.concatenate([WALSH[i][:] for j in range(4)])

## ---- tri-bit scramble sequence for preamble ----------------------------------
TRIBIT_SCRAMBLE = np.array(
    [7,4,3,0,5,1,5,0,2,2,1,1,5,7,4,3,5,0,2,6,2,1,6,2,0,0,5,0,5,2,6,6],
    dtype=np.uint8)

def n_psk(n,x):
    return np.complex64(np.exp(2j*np.pi*x/n))

## ---- preamble symbols ---------------------------------------------------------
D1=D2=C1=C2=C3=0 ## not known
PRE_SYMBOLS  = n_psk(2, np.concatenate(
    [TRIBIT[i][:] for i in [0,1,3,0,1,3,1,2,0,D1,D2,C1,C2,C3,0]]))
PRE_SYMBOLS[9*32:14*32] = 0

## ---- preamble scramble symbols ------------------------------------------------
PRE_SCRAMBLE = n_psk(8, np.concatenate([TRIBIT_SCRAMBLE for i in range(15)]))

## ---- data scrambler -----------------------------------------------------------
class ScrambleData(object):
    """data scrambling sequence generator"""
    def __init__(self):
        self.reset()

    def reset(self):
        self._state = 0xBAD
        self._counter = 0

    def next(self):
        if self._counter == 160:
            self.reset()
        for j in range(8):
            self._advance()
        self._counter += 1
        return self._state&7

    def _advance(self):
        msb = self._state>>11
        self._state = (self._state<<1)&4095
        if msb:
            self._state ^= 0x053
        return self._state

## ---- constellation indices ---------------------------------------------------
MODE_BPSK=0
MODE_QPSK=1
MODE_8PSK=2

## ---- mode definitions --------------------------------------------------------
MODE = [[{} for x in range(8)] for y  in range(8)]
MODE[7][6] = {'bit_rate':4800, 'ci':MODE_8PSK, 'interleaver':['N',  1,  1], 'unknown':32,'known':16, 'nsymb': 1, 'coding_rate': -1 }
MODE[7][7] = {'bit_rate':2400, 'ci':MODE_8PSK, 'interleaver':['N',  1,  1], 'unknown':32,'known':16, 'nsymb': 1, 'coding_rate':1./2}

MODE[6][4] = {'bit_rate':2400, 'ci':MODE_8PSK, 'interleaver':['S', 40, 72], 'unknown':32,'known':16, 'nsymb': 1, 'coding_rate':1./2}
MODE[4][4] = {'bit_rate':2400, 'ci':MODE_8PSK, 'interleaver':['L', 40,576], 'unknown':32,'known':16, 'nsymb': 1, 'coding_rate':1./2}

MODE[6][5] = {'bit_rate':1200, 'ci':MODE_QPSK, 'interleaver':['S', 40, 36], 'unknown':20,'known':20, 'nsymb': 1, 'coding_rate':1./2}
MODE[4][5] = {'bit_rate':1200, 'ci':MODE_QPSK, 'interleaver':['L', 40,288], 'unknown':20,'known':20, 'nsymb': 1, 'coding_rate':1./2}

MODE[6][6] = {'bit_rate': 600, 'ci':MODE_BPSK, 'interleaver':['S', 40, 18], 'unknown':20,'known':20, 'nsymb': 1, 'coding_rate':1./2}
MODE[4][6] = {'bit_rate': 600, 'ci':MODE_BPSK, 'interleaver':['L', 40,144], 'unknown':20,'known':20, 'nsymb': 1, 'coding_rate':1./2}

MODE[6][7] = {'bit_rate': 300, 'ci':MODE_BPSK, 'interleaver':['S', 40, 18], 'unknown':20,'known':20, 'nsymb': 1, 'coding_rate':1./4}
MODE[4][7] = {'bit_rate': 300, 'ci':MODE_BPSK, 'interleaver':['L', 40,144], 'unknown':20,'known':20, 'nsymb': 1, 'coding_rate':1./4}

MODE[7][4] = {'bit_rate': 150, 'ci':MODE_BPSK, 'interleaver':['S', 40, 18], 'unknown':20,'known':20, 'nsymb': 1, 'coding_rate':1./8}
MODE[5][4] = {'bit_rate': 150, 'ci':MODE_BPSK, 'interleaver':['L', 40,144], 'unknown':20,'known':20, 'nsymb': 1, 'coding_rate':1./8}

MODE[7][5] = {'bit_rate':  75, 'ci':MODE_QPSK, 'interleaver':['S', 10,  9], 'unknown':-1,'known': 0, 'nsymb':32, 'coding_rate':1./2}
MODE[5][4] = {'bit_rate':  75, 'ci':MODE_QPSK, 'interleaver':['L', 20, 36], 'unknown':-1,'known': 0, 'nsymb':32, 'coding_rate':1./2}

## ---- physcal layer class -----------------------------------------------------
class PhysicalLayer(object):
    """Physical layer description for MIL-STD-188-110 Appendix A"""

    def __init__(self, sps):
        """intialization"""
        self._sps     = sps
        self._frame_counter = 0
        self._is_first_frame = True
        self._constellations = [self.make_psk(2, [0,1]),
                                self.make_psk(4, [0,1,3,2]),
                                self.make_psk(8, [0,1,3,2,7,6,4,5])] ## TODO: check 8PSK gray code
        self._preamble    = self.get_preamble()
        self._pre_counter = -1
        self._d1d2        = [-1,-1] ## D1,D2
        self._mode        = {}
        self._scr_data    = ScrambleData()
        ##self._data     = self.get_data()

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
              self._pre_counter, self._frame_counter)
        if self._pre_counter != 0:
            self._scr_data.reset()
            return [self._preamble,MODE_BPSK,True,False]
        num_symb = 11520 if self._mode['interleaver'][0] == 'L' else 1440
        a = np.zeros(num_symb, dtype=[('symb',     np.complex64),
                                      ('scramble', np.complex64)])
        n_known = self._mode['known']
        n_unknown = self._mode['unknown']
        counter_d1d2 = 0
        for i in range(0,num_symb,n_known+n_unknown):
            a['symb'][i          :i+n_unknown        ] = 0
            a['symb'][i+n_unknown:i+n_unknown+n_known] = 1
            if i>=num_symb-2*(n_unknown+n_known):
                a['symb'][i+0:i+ 8] *= n_psk(2, WALSH[self._d1d2[counter_d1d2]][:])
                a['symb'][i+8:i+16] *= n_psk(2, WALSH[self._d1d2[counter_d1d2]][:])
                counter_d1d2 += 1

        a['scramble'] = n_psk(8, np.array([self._scr_data.next() for _ in range(num_symb)]))
        a['symb']    *= a['scramble']
        self._frame_counter += 1
        return [a, self._mode['ci'],False,True]

    def get_doppler(self, symbols, iq_samples):
        """returns a tuple
        [0] ... quality flag
        [1] ... doppler estimate (rad/symbol) if available"""
        print('-------------------- get_doppler --------------------',
              self._frame_counter,len(symbols),len(iq_samples))
        success = False
        doppler = 0
        if self._frame_counter == 0:
            success,doppler = self.quality_preamble(symbols,iq_samples)
            if len(symbols) != 0:
                data = [FROM_WALSH[walsh_to_num
                                    (np.real
                                     (np.sum
                                      (symbols[i:i+32].reshape((4,8)),0))<0)]
                        for i in range(0,15*32,32)]
                print('data=',data)
                self._pre_counter = sum((np.array(data[11:14])&3)
                                        *(1<<2*np.arange(3)[::-1]))
                self._d1d2 = data[9:11]
                self._mode = MODE[data[9]][data[10]]
                print('pre_counter', self._pre_counter, 'mode', self._mode)
                self._is_first_frame = not success
                success = True
        else:
            for i in range(0,len(symbols),40):
                print(i,symbols[i:i+40], np.mean(np.abs(symbols[i:i+40])))
            success = np.mean(np.abs(symbols[0:40])) > 0.5
            if not success:
                self._frame_counter = 0
                self._pre_counter = -1
        return success,doppler

    def is_preamble(self):
        return self._frame_counter == 0

    def quality_preamble(self, symbols, iq_samples):
        """quality check and doppler estimation for preamble"""
        success = True
        doppler = 0
        if len(iq_samples) != 0:
            zp   = np.conj(self.get_preamble_z(self._sps))[9*self._sps:]
            cc   = np.array([np.sum(iq_samples[i*self._sps:(3*32+i-9)*self._sps]*zp)
                             for i in range(4*32)])
            acc = np.abs(cc)
            for i in range(0,len(cc),32):
                print('i=%3d: '%i,end='')
                for j in range(32):
                    print('%3.0f ' % acc[i+j], end='')
                print()
            imax = np.argmax(np.abs(cc[0:2*32]))
            pks  = cc[(imax,imax+3*16,imax+3*16+1,imax+3*32),]
            apks = np.abs(pks)
            print('imax=', imax, 'apks=',apks)
            success = np.mean(apks[(0,3),]) > 2*np.mean(apks[(1,2),])
            doppler = np.diff(np.unwrap(np.angle(pks[(0,3),])))[0]/(3*32) if success else 0
            print('success=', success, 'doppler=', doppler)
        #if len(symbols) != 0:
        ## TODO: check the symbols
        return success,doppler

    @staticmethod
    def get_preamble():
        """preamble symbols + scrambler"""
        a=np.zeros(15*32, dtype=[('symb',     np.complex64),
                                 ('scramble', np.complex64)])
        a['symb']     = PRE_SCRAMBLE*PRE_SYMBOLS
        a['scramble'] = PRE_SCRAMBLE
        return a

    @staticmethod
    def get_preamble_z(sps):
        """preamble symbols for preamble correlation"""
        a = PhysicalLayer.get_preamble()
        return np.array([z for z in a['symb'][0:32*3]
                         for i in range(sps)])

    @staticmethod
    def make_psk(n, gray_code):
        """generates n-PSK constellation data"""
        c = np.zeros(n, dtype=[('points',  np.complex64),
                               ('symbols', np.uint8)])
        c['points']  = n_psk(n,np.arange(n))
        c['symbols'] = gray_code
        return c

if __name__ == '__main__':
    def gen_data_scramble():
        def advance(s):
            msb = s>>11
            s = (s<<1)&((1<<12)-1)
            if msb: s ^= 0x053
            return s
        a = np.zeros(160, dtype=np.uint8)
        s = 0xBAD
        for i in range(160):
            for j in range(8): s = advance(s)
            a[i] = s&7;
        return a

    p=PhysicalLayer(5)
    z1=np.array([x for x in PRE_SYMBOLS  for i in range(5)])
    z2=np.array([x for x in PRE_SCRAMBLE for i in range(5)])
    z=z1*z2

    for i in range(3):
        print(i, all(z[32*5*i:32*5*(i+1)] == z[32*5*(3+i):32*5*(3+i+1)]))

    print(np.sum(np.sum(z[0:32*5] * np.conj(z[32*5*3:32*5*4]))))
    print(WALSH[1][:])
    print(sum(WALSH[1][:]*(1<<np.array(range(7,-1,-1)))))
    print(FROM_WALSH)
    print(gen_data_scramble())

    s=ScrambleData()
    print(type(s))
    print([s.next() for _ in range(160)])
    print([s.next() for _ in range(160)])
