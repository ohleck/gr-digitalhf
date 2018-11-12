## -*- python -*-

import numpy as np

CONST_DTYPE=np.dtype([('points',  np.complex64),
                      ('symbols', np.int32)])

def n_psk(n,x):
    """n-ary PSK constellation"""
    return np.complex64(np.exp(2j*np.pi*x/n))

def freq_est(z):
    """Data-Aided Frequency Estimation for Burst Digital Transmission,
        Umberto Mengali and M. Morelli, IEEE TRANSACTIONS ON COMMUNICATIONS,
        VOL. 45, NO. 1, JANUARY 1997"""
    L0 = len(z)
    N  = L0//2
    R  = np.zeros(N+1, dtype=np.complex64)
    for i in range(N+1):
        R[i] = 1.0/(L0-i)*np.sum(z[i:]*np.conj(z[0:L0-i])) ## eq (3)
    m  = np.arange(N+1, dtype=np.float32)
    w  = 3*((L0-m)*(L0-m+1)-N*(L0-N))/(N*(4*N*N - 6*N*L0 + 3*L0*L0-1)) ## eq (9)
    mod_2pi = lambda x : np.mod(x-np.pi, 2*np.pi) - np.pi
    return np.sum(w[1:] * mod_2pi(np.diff(np.angle(R))))   ## eq (8)

if __name__ == '__main__':
    idx=np.arange(3)
    z=np.exp(1j*idx*0.056+1j)
    print(freq_est(z)/0.056)
