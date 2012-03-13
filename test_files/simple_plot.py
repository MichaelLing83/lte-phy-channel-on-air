#@+leo-ver=5-thin
#@+node:michael.20120305092148.1300: * @thin ./test_files/simple_plot.py
#@+others
#@+node:michael.20120305092148.1298: ** test
# time scale is in 1 ms
T_s = 1.0/30720 # in ms
f_0 = (2620+0.1*(2620-2750))*1000*1000  # in kHz

#@+others
#@+node:michael.20120305092148.1293: *3* 6.12 OFDM baseband signal gen
def s_p_l(symbol_array, l, N_DL_RB, N_RB_sc, N_DL_CP, delta_f=15000):
    '''
    Note: len(symbol_array)==N_DL_RB*N_RB_sc must be True.
    '''
    T_s = 1./30720/1000  # all time scale is in 1 s
    if N_DL_CP==0 and delta_f==15000:
        if l==0:
            N_CP_l = 160
        else:
            N_CP_l = 144
    elif N_DL_CP==1:    # extended CP
        if delta_f==15000:
            N_CP_l = 512
        else:   # delta_f == 7500
            N_CP_l = 1024
    if delta_f==15000:
        N = 2048
    else:   # delta_f == 7500
        N = 4096
    t = arange(0, (N_CP_l+N)*T_s, T_s)
    signal_pl = 0
    down_limit = int(floor(N_DL_RB*N_RB_sc/2))
    up_limit = int(ceil(N_DL_RB*N_RB_sc/2))
    for k in arange( -1*down_limit, 0, 1 ):
        signal_pl += symbol_array[k+down_limit]*exp(1j*2*pi*k*delta_f*(t-N_CP_l*T_s))
    for k in arange(1, up_limit+1, 1):
        signal_pl += symbol_array[k+down_limit-1]*exp(1j*2*pi*k*delta_f*(t-N_CP_l*T_s))
    return signal_pl
#@+node:michael.20120305092148.1296: *3* 6.13 Modulation&upconversion
def downlink_modulate(s_p_l, t, f_0):
    return cos(2*pi*f_0*t) * s_p_l.real - sin(2*pi*f_0*t) * imag(s_p_l)

def downlink_downconvert(signal, t, f_0):
    
    cutoff_freq = f_0
    nyq = 2 * f_0
    numtaps = 80
    lp_fir = firwin(numtaps, cutoff_freq, window=('kaiser',8), nyq=nyq)
    
    I = -2* convolve( signal * cos(2*pi*f_0*t), lp_fir )[numtaps/2:len(signal)+numtaps/2]
    Q = -2 * convolve( signal * sin(2*pi*f_0*t), lp_fir )[numtaps/2:len(signal)+numtaps/2]
    
    return I + 1j*Q
#@+node:michael.20120305092148.1292: *3* plot_symbols
from numpy import *
import matplotlib.pyplot as plt

def plot_symbol(symbol_array, l, N_DL_CP, delta_f=15000):
    T_s = 1./30720  # all time scale is in 1 ms
    if N_DL_CP==0 and delta_f==15000:
        if l==0:
            N_CP_l = 160
        else:
            N_CP_l = 144
    elif N_DL_CP==1:    # extended CP
        if delta_f==15000:
            N_CP_l = 512
        else:   # delta_f == 7500
            N_CP_l = 1024
    if delta_f==15000:
        N = 2048
    else:   # delta_f == 7500
        N = 4096
    t = arange(0, (N_CP_l+N)*T_s, T_s)
    # use gnuplot
    plt.plot(t, symbol_array)
    plt.show()

def myplot(sig, t):
    plt.cla()
    plt.plot(t, sig)
    plt.xlabel('time (ms)')
    plt.show()
#@-others

l = 0
N_DL_RB = 110
N_RB_sc = 12
N_DL_CP = 0

symbol_array = array([0]*N_DL_RB*N_RB_sc)
symbol_array[0] = 1
symbol_array[-1] = 1

baseband_sig = s_p_l(symbol_array, l, N_DL_RB, N_RB_sc, N_DL_CP, delta_f=15000)
t = arange(0, (160+2048)*T_s, T_s)
uu_sig = downlink_modulate(baseband_sig, t, f_0)
myplot(uu_sig, t)
#@-others
#@-leo
