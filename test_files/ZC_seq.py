#@+leo-ver=5-thin
#@+node:michael.20120305092148.1314: * @thin ./test_files/ZC_seq.py
#@+others
#@+node:michael.20120305092148.1313: ** source_code
from numpy import *
import matplotlib.pyplot as plt

# time scale is in 1 ms
T_s = 1.0/30720 # in ms
f_0 = (2620+0.1*(2620-2750))*1000*1000  # in kHz

#@+others
#@+node:michael.20120305092148.1316: *3* Zadoff-Chu seq
def ZC( n, N_ZC, q ):
    '''
    give the n-th element of 0 cyclic shift Z-C sequence with root index q and length N_ZC.
    '''
    return exp(-1j*2*pi*q*n*(n+1)/2/N_ZC)
#@-others

l = 0
N_DL_RB = 110
N_RB_sc = 12
N_DL_CP = 0 # normal DL CP
N_ZC = 63

zc_seq_d = dict()
for root_index in (29,34,25):
    zc_seq_d[root_index] = array([0]*N_ZC, dtype=complex128)
    for i in arange(N_ZC):
        zc_seq_d[root_index][i] = ZC(i, N_ZC, root_index)

cs_corr = array([0]*N_ZC)
for p in (29,34,25):
    for q in (29,34,25):
        for cs in arange(N_ZC):
            cs_corr[cs] = abs(correlate(zc_seq_d[p],roll(zc_seq_d[q],cs))[0])/N_ZC
            plt.cla()
            plt.axis( [-1, N_ZC+1, -0.1, 1.1] )
            plt.plot(cs_corr, marker='+', linestyle='-')
            fn = "cyclic_shift_corr_ZC_%s_%s.png"%(p,q)
            plt.savefig(fn)

#@-others
#@-leo
