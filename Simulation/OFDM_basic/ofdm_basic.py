#@+leo-ver=5-thin
#@+node:michael.20120315195140.1454: * @thin ./Simulation/OFDM_basic/ofdm_basic.py
#@+others
#@+node:michael.20120305092148.1310: ** source
from scipy.signal import *
from numpy import *
import matplotlib.pyplot as plt

# time scale is in 1 s
T_s = 1.0/30720/1000 # in s

# configuration for CSRS
n_s = 0
l = 0
antenna_port = 0
N_DL_RB = 110
N_maxDL_RB = N_DL_RB
N_RB_sc = 12
N_DL_CP = 0 # normal DL CP
DL_CP_type = 0
N_DL_symb = 7
N_ID_2_tuple = (0,1,2)
delta_f = 15000
subframe = 0
N_cell_ID = 0
f_0 = (2620+0.1*(2620-2750))*1000*1000  # in Hz

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

def find_max( a_list ):
    m = max(a_list)
    for i in arange(len(a_list)):
        if a_list[i] == m:
            return (i, m)

def find_min( a_array ):
    x, y = 0, 0
    for i in arange(len(a_array)):
        if a_array[i] < y:
            x, y = i, a_array[i]
    return (x,y)

def find_abs_max( a_array ):
    m = max(abs(a_array))
    for i in arange(len(a_array)):
        if abs(a_array[i]) == m:
            return (i, m)

            
#@+others
#@+node:michael.20120315195140.1455: *3* 01. OFDM baseband signal generation
def ofdm_baseband_signal_generation():
    
    symbol_array = array( [0.0 + 0.0 * 1j]*(N_DL_RB*N_RB_sc) )
    symbol_array[0] = 1
    
    ofdm_baseband_IQ_direct = s_p_l(symbol_array, l, N_DL_RB, N_RB_sc, N_DL_CP, delta_f=15000, gen_method='DIRECT')
    ofdm_baseband_IQ_ifft = s_p_l(symbol_array, l, N_DL_RB, N_RB_sc, N_DL_CP, delta_f=15000, gen_method='IFFT')
    
    plt.plot(real(ofdm_baseband_IQ_ifft))
    #plt.plot(real(ofdm_baseband_IQ_direct))
    #plt.plot(real(ofdm_baseband_IQ_ifft)-real(ofdm_baseband_IQ_direct))
    #plt.plot(real(ofdm_baseband_IQ_ifft))
    #plt.plot(fft.fft(ofdm_baseband_IQ_ifft[-1*N:], N))
    plt.show()
#@+node:michael.20120305092148.1293: *3* 6.12 OFDM baseband signal gen
def s_p_l(symbol_array, l, N_DL_RB, N_RB_sc, N_DL_CP, delta_f=15000, gen_method='DIRECT'):
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
    signal_pl =  array([0.0+0.0*1j] * (N_CP_l + N))
    
    down_limit = int(floor(N_DL_RB*N_RB_sc/2))
    up_limit = int(ceil(N_DL_RB*N_RB_sc/2))
    
    if gen_method == 'DIRECT':
        for k in arange( -1*down_limit, 0, 1 ):
            signal_pl += symbol_array[k+down_limit]*exp(1j*2*pi*k*delta_f*(t-N_CP_l*T_s))
        for k in arange(1, up_limit+1, 1):
            signal_pl += symbol_array[k+down_limit-1]*exp(1j*2*pi*k*delta_f*(t-N_CP_l*T_s))
    elif gen_method == 'IFFT':
        mapped_seq = array([0.0+0.0*1j] * N)
        mapped_seq[:down_limit] = symbol_array[:down_limit]
        mapped_seq[down_limit] = 0.0 + 0.0 * 1j
        mapped_seq[down_limit+1:down_limit+up_limit+1] = symbol_array[down_limit:]
        #for i in arange(down_limit+up_limit+1, N):
            #mapped_seq[i] = 0.0 + 0.0 * 1j
        signal_pl[N_CP_l:] = fft.ifft(mapped_seq, N) * N * exp(1j*2*pi*down_limit*delta_f*t[-1*N:])
        #signal_pl[:N_CP_l] = signal_pl[-1*N_CP_l:]
        
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

test_enabling_bits = 0b1

# 01. OFDM baseband signal generation
if test_enabling_bits & (1<<0):
    ofdm_baseband_signal_generation()
#@-others
#@-leo
