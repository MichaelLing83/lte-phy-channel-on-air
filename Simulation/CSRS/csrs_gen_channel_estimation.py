#@+leo-ver=5-thin
#@+node:Michael.20120315095133.1440: * @thin ./Simulation/CSRS/csrs_gen_channel_estimation.py
#@+others
#@+node:Michael.20120315095133.1439: ** source

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
#@+node:Michael.20120315095133.1442: *3* 01. CSRS baseband IQ time domain signal
def csrs_baseband_IQ_time_domain_signal():
    
    csrs_re_array = get_CSRS_in_symbol(n_s, l, antenna_port, N_cell_ID, N_maxDL_RB, N_DL_RB, N_RB_sc, N_DL_symb, DL_CP_type)
    csrs_baseband_IQ = ofdm_baseband_IQ_signal_generate(csrs_re_array, l, N_DL_RB, N_RB_sc, N_DL_CP, delta_f)
    
    subplot_pos_tupe = (131,132,133)
    title_tuple = ('CSRS baseband IQ OFDM magnitude','CSRS baseband IQ OFDM signal real part','CSRS baseband IQ OFDM signal imag part')
    y_label_tuple = ('IQ Magnitude', 'I part', 'Q part')
    func_tuple = (abs, real, imag)
    legend_list = list()
    for i in (0,1,2):
        plt.subplot(subplot_pos_tupe[i])
        plt.title(title_tuple[i])
        plt.plot(t*1000, func_tuple[i](csrs_baseband_IQ))
        legend_list.append( ('antenna_port=%s'%(antenna_port), ))
        plt.xlabel('Time (ms)')
        plt.ylabel(y_label_tuple[i])
        #plt.axis([-0.01, 0.075, 0, 15])
        #plt.legend( ('N_ID_cell=%s'%N_cell_ID,) )
            
    plt.show()
    
#@+node:Michael.20120315095133.1443: *3* 02. CSRS baseband IQ spectrum
def csrs_baseband_IQ_spectrum(to_draw=True):
    
    csrs_re_array = get_CSRS_in_symbol(n_s, l, antenna_port, N_cell_ID, N_maxDL_RB, N_DL_RB, N_RB_sc, N_DL_symb, DL_CP_type)
    csrs_baseband_IQ = ofdm_baseband_IQ_signal_generate(csrs_re_array, l, N_DL_RB, N_RB_sc, N_DL_CP, delta_f)[-1*N:]
    csrs_baseband_IQ_fft = fft.fft(csrs_baseband_IQ, N)
    
    if to_draw:
        legend_list = list()
        plt.title('CSRS baseband IQ spectrum')
        legend_list.append( 'Spectrum magnitude' )
        plt.plot(abs(csrs_baseband_IQ_fft), linestyle='-')
        plt.xlabel('n (FFT index)')
        plt.ylabel('Spectrum magnitude')
        plt.legend(legend_list)
        plt.show()
    
    return csrs_baseband_IQ_fft
#@+node:Michael.20120315095133.1444: *3* 03. CSRS Uu signal
def CSRS_signal_Uu_ref(n_s, l, antenna_port, N_cell_ID, N_maxDL_RB, N_DL_RB, N_RB_sc, N_DL_symb, DL_CP_type):
    
    csrs_re_array = get_CSRS_in_symbol(n_s, l, antenna_port, N_cell_ID, N_maxDL_RB, N_DL_RB, N_RB_sc, N_DL_symb, DL_CP_type)
    csrs_baseband_IQ = ofdm_baseband_IQ_signal_generate(csrs_re_array, l, N_DL_RB, N_RB_sc, N_DL_CP, delta_f)
    csrs_Uu_signal = downlink_modulate(csrs_baseband_IQ, t, f_0)
    
    legend_list = list()
    plt.plot(t*1000, csrs_Uu_signal)
    legend_list.append('CSRS Uu signal')
    plt.title('CSRS Uu signal')
    plt.xlabel('Time (ms)')
    plt.ylabel('Signal level')
    plt.legend(legend_list)
    #plt.axis( [-0.01, 0.075, -0.1, 14] )
    plt.show()
    #plt.savefig('PSS_signal_Uu_for_N_ID_2=%s.png'%N_ID_2, dpi=300)

def CSRS_signal_Uu():
    #CSRS_signal_Uu_ref(n_s, l, antenna_port, N_cell_ID, N_maxDL_RB, N_DL_RB, N_RB_sc, N_DL_symb, DL_CP_type)
    CSRS_signal_Uu_ref(0, 0, 0, 0, 110, 110, 12, 7, 0)
#@+node:Michael.20120315095133.1445: *3* 04. CSRS received IQ
def CSRS_received_IQ_ref(n_s, l, antenna_port, N_cell_ID, N_maxDL_RB, N_DL_RB, N_RB_sc, N_DL_symb, DL_CP_type, to_draw=False):
    
    csrs_re_array = get_CSRS_in_symbol(n_s, l, antenna_port, N_cell_ID, N_maxDL_RB, N_DL_RB, N_RB_sc, N_DL_symb, DL_CP_type)
    csrs_baseband_IQ = ofdm_baseband_IQ_signal_generate(csrs_re_array, l, N_DL_RB, N_RB_sc, N_DL_CP, delta_f)
    csrs_Uu_signal = downlink_modulate(csrs_baseband_IQ, t, f_0)
    csrs_Uu_signal_downconverted = downlink_downconvert(csrs_Uu_signal, t, f_0)
    
    if to_draw:
        legend_list = list()
        plt.plot(t*1000, real(csrs_Uu_signal_downconverted))
        legend_list.append('CSRS received IQ real part')
        plt.title('CSRS received IQ')
        plt.xlabel('Time (ms)')
        plt.ylabel('Signal level')
        plt.legend(legend_list)
        #plt.axis( [-0.01, 0.075, -0.1, 14] )
        plt.show()
    return csrs_Uu_signal_downconverted

def CSRS_received_IQ():
    #CSRS_received_IQ_ref(n_s, l, antenna_port, N_cell_ID, N_maxDL_RB, N_DL_RB, N_RB_sc, N_DL_symb, DL_CP_type)
    CSRS_received_IQ_ref(0, 0, 0, 0, 110, 110, 12, 7, 0, to_draw=True)
#@+node:Michael.20120315095133.1446: *3* 05. CSRS channel estimation in one symbol
def csrs_channel_estimation_in_one_symbol(received_baseband_IQ, ref_csrs_baseband_IQ, to_draw=False):
    
    received_baseband_IQ_fft = fft.fft( received_baseband_IQ[-1*N:], N )
    ref_csrs_baseband_IQ_fft = fft.fft(ref_csrs_baseband_IQ[-1*N:], N)
    channel_estimation_in_symbol = received_baseband_IQ_fft / ref_csrs_baseband_IQ_fft
    
    
    subplot_pos_tupe = (121,122)
    title_tuple = ('channel estimation spectrum magnitude','channel estimation spectrum phase')
    y_label_tuple = ('spectrum magnitude', 'spectrum phase')
    func_tuple = (abs, arctan)
    legend_list = list()
    for i in (0,1):
        plt.subplot(subplot_pos_tupe[i])
        plt.title(title_tuple[i])
        plt.plot(func_tuple[i](channel_estimation_in_symbol))
        #legend_list.append( ('antenna_port=%s'%(antenna_port), ))
        plt.xlabel('n (FFT index)')
        plt.ylabel(y_label_tuple[i])
        #plt.axis([-0.01, 0.075, 0, 15])
        #plt.legend( ('N_ID_cell=%s'%N_cell_ID,) )
            
    plt.show()

def test_csrs_channel_estimation_in_one_symbol():
    n_s, l, antenna_port, N_cell_ID, N_maxDL_RB, N_DL_RB, N_RB_sc, N_DL_symb, DL_CP_type = 0, 0, 0, 0, 110, 110, 12, 7, 0
    csrs_re_array = get_CSRS_in_symbol(n_s, l, antenna_port, N_cell_ID, N_maxDL_RB, N_DL_RB, N_RB_sc, N_DL_symb, DL_CP_type)
    csrs_baseband_IQ = ofdm_baseband_IQ_signal_generate(csrs_re_array, l, N_DL_RB, N_RB_sc, N_DL_CP, delta_f)
    csrs_Uu_signal = downlink_modulate(csrs_baseband_IQ, t, f_0)
    csrs_Uu_signal_downconverted = downlink_downconvert(csrs_Uu_signal, t, f_0)
    
    csrs_channel_estimation_in_one_symbol(csrs_Uu_signal_downconverted, csrs_baseband_IQ, True)
#@+node:michael.20120305092148.1285: *3* 6.10.1.1 Seq gen
def r_l_ns(n_s, l, N_cell_ID, N_maxDL_RB, DL_CP_type):
    '''
    r_l_ns(l, n_s, N_cell_ID, N_maxDL_RB, DL_CP_type): list of complex symbols for CSRS signal in symbol index l of given slot.
    l: symbol index in given slot
    n_s: slot index
    N_cell_ID:  cell ID
    N_maxDL_RB: 110 for 20MHz
    DL_CP_type: CP type for downlink, 0 for normal CP and 1 for extended CP
    
    profile:
import cProfile
from math import sqrt
def main():
    for i in range(100):
        tmp = r_l_ns(i%20, i%7, i, 110, i%2)

cProfile.runctx('main()', globals(), locals())

   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
    44000    0.227    0.000  109.216    0.002 <ipython console>:1(c)
        1    0.001    0.001  109.585  109.585 <ipython console>:1(main)
      100    0.255    0.003  109.583    1.096 <ipython console>:1(r_l_ns)
    44000   46.434    0.001   46.434    0.001 <ipython console>:1(x_1)
    44000   62.555    0.001   62.555    0.001 <ipython console>:1(x_2)
    '''
    if DL_CP_type == 0: # normal DL CP
        N_CP = 1
    else:
        N_CP = 0
    c_init = 2**10 * (7*(n_s+1)+l+1) * (2*N_cell_ID+1) + 2*N_cell_ID + N_CP
    csrs_symbol_list = list()
    for m in range(2*N_maxDL_RB):
        real_part = 1/sqrt(2) * (1-2*c(c_init,2*m))
        image_part = 1/sqrt(2) * (1-2*c(c_init,2*m+1))
        csrs_symbol_list.append( complex(real_part,image_part) )
    return tuple(csrs_symbol_list)
#@+node:michael.20120305092148.1279: *3* 6.10.1.2 Mapping to REs
def get_CSRS_REs_in_slot(n_s, antenna_port, N_cell_ID, N_maxDL_RB, N_DL_RB, N_DL_symb):
    '''
    get_CSRS_REs_in_slot(n_s, antenna_port, N_cell_ID, N_maxDL_RB, N_DL_RB, N_DL_symb): tuple of CSRS REs in the specified symbol of RB.
    n_s: slot index
    antenna_port: antenna port for CSRS
    N_cell_ID: cell ID
    N_maxDL_RB: 110 for 20MHz configured by higher layer
    N_DL_RB: PHY number of downlink RB
    N_DL_symb: maximum 110 for 20MHz
    
    profile:
def main():
    for i in range(1000):
        tmp = get_CSRS_in_slot(i%20, i%4, i, 110, 110, (7,3,2,4,5,1)[i%6])

cProfile.runctx('main()', globals(), locals())

   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
        1    0.012    0.012    0.903    0.903 <ipython console>:1(main)
     1000    0.536    0.001    0.891    0.001 <ipython console>:2(get_CSRS_in_slot)
    '''
    
    REs = list()
    # symbol indices for CSRS of this AP
    if antenna_port in (0,1):
        if N_DL_symb>3:
            l_list = (0, N_DL_symb-3)
        else:   # DwPTS that has only 3 DL symbols
            l_list = (0,)
    else:   # antenna_port in (2,3)
        l_list = (1,)
    # v_shift
    v_shift = N_cell_ID % 6
    for l in l_list:
        # v
        if antenna_port==0 and l==0:
            v = 0
        elif antenna_port==0 and l!=0:
            v = 3
        elif antenna_port==1 and l==0:
            v = 3
        elif antenna_port==1 and l!=0:
            v = 0
        elif antenna_port==2:
            v = 3 * (n_s%2)
        elif antenna_port==3:
            v = 3 + 3 * (n_s%2)
        for m in range(2*N_DL_RB-1):
            m_ = m + N_maxDL_RB - N_DL_RB   # m'
            k = 6*m + (v+v_shift)%6
            REs.append( (k,l) )
    return tuple(REs)

def get_CSRS_in_symbol(n_s, l, antenna_port, N_cell_ID, N_maxDL_RB, N_DL_RB, N_RB_sc, N_DL_symb, DL_CP_type):
    '''
    '''
    symbol_array = ndarray(shape=(N_maxDL_RB*N_RB_sc,),dtype=complex128)
    for i in arange(len(symbol_array)):
        symbol_array[i] = 0
    csrs_seq = r_l_ns(n_s, l, N_cell_ID, N_maxDL_RB, DL_CP_type)
    # symbol indices for CSRS of this AP
    if antenna_port in (0,1):
        if N_DL_symb>3:
            l_list = (0, N_DL_symb-3)
        else:   # DwPTS that has only 3 DL symbols
            l_list = (0,)
    else:   # antenna_port in (2,3)
        l_list = (1,)
    # v_shift
    v_shift = N_cell_ID % 6
    if l in l_list:
        # v
        if antenna_port==0 and l==0:
            v = 0
        elif antenna_port==0 and l!=0:
            v = 3
        elif antenna_port==1 and l==0:
            v = 3
        elif antenna_port==1 and l!=0:
            v = 0
        elif antenna_port==2:
            v = 3 * (n_s%2)
        elif antenna_port==3:
            v = 3 + 3 * (n_s%2)
        for m in range(2*N_DL_RB):
            m_ = m + N_maxDL_RB - N_DL_RB   # m'
            k = 6*m + (v+v_shift)%6
            symbol_array[k] = csrs_seq[m_]
    return symbol_array
#@+node:michael.20120305092148.1293: *3* 6.12 OFDM baseband signal gen
def ofdm_baseband_IQ_signal_generate(symbol_array, l, N_DL_RB, N_RB_sc, N_DL_CP, delta_f=15000, gen_method='IFFT'):
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
        mapped_seq = map_baseband_IQ_for_ifft(symbol_array)
        signal_pl[N_CP_l:] = fft.ifft(mapped_seq, N) * N
        signal_pl[:N_CP_l] = signal_pl[-1*N_CP_l:]
        
    return signal_pl

def map_baseband_IQ_for_ifft(baseband_IQ_array):
    '''
    Note: len(symbol_array)==N_DL_RB*N_RB_sc must be True.
    '''
    #T_s = 1./30720/1000  # all time scale is in 1 s
    if delta_f==15000:
        N = 2048
    else:   # delta_f == 7500
        N = 4096
    #t = arange(0, (N_CP_l+N)*T_s, T_s)
    #signal_pl =  array([0.0+0.0*1j] * (N_CP_l + N))
    
    down_limit = int(floor(len(baseband_IQ_array)/2))
    up_limit = int(ceil(len(baseband_IQ_array)/2))
    
    mapped_seq = array([0.0+0.0*1j] * N)
    # do the mapping before IFFT
    tmp_index = N-1
    for i in arange(down_limit-1, -1, -1):
        mapped_seq[tmp_index] = baseband_IQ_array[i]
        tmp_index -= 1
    tmp_index = 1
    for i in arange(down_limit, down_limit+up_limit):
        mapped_seq[tmp_index] = baseband_IQ_array[i]
        tmp_index += 1
    return mapped_seq

def map_fft_result_to_RE_IQ_array(fft_result_array):
    '''
    map_fft_result_to_baseband_IQ(fft_result_array): baseband_IQ_array
    Note: len(fft_result_array)==N must be True
            len(baseband_IQ_array) is N_DL_RB*N_RB_sc
    '''
    if delta_f==15000:
        N = 2048
    else:   # delta_f == 7500
        N = 4096
    
    mapped_seq = array([0.0+0.0*1j] * (N_DL_RB*N_RB_sc))
    
    down_limit = int(floor(len(mapped_seq)/2))
    up_limit = int(ceil(len(mapped_seq)/2))
    
    tmp_index = N-1
    for i in arange(down_limit-1, -1, -1):
        mapped_seq[i] = fft_result_array[tmp_index]
        tmp_index -= 1
    tmp_index = 1
    for i in arange(down_limit, down_limit+up_limit):
        mapped_seq[i] = fft_result_array[tmp_index]
        tmp_index += 1
        
    return mapped_seq
    
def ofdm_baseband_IQ_to_RE_IQ_array(baseband_IQ_array, N_DL_RB, N_RB_sc, delta_f=15000):
    '''
    Note: len(baseband_IQ_array)==N must be True.
    '''
    if delta_f==15000:
        N = 2048
    else:   # delta_f == 7500
        N = 4096

    re_IQ_array =  array([0.0+0.0*1j] * (N_DL_RB * N_RB_sc))
    re_IQ_array = 1.0/N * map_fft_result_to_RE_IQ_array(fft.fft(baseband_IQ_array, N))
        
    return re_IQ_array
#@+node:michael.20120305092148.1296: *3* 6.13 Modulation&upconversion
def downlink_modulate(s_p_l, t, f_0):
    modulated_signal = cos(2*pi*f_0*t) * s_p_l.real - sin(2*pi*f_0*t) * imag(s_p_l)
    cutoff_freq = f_0
    nyq = 2 * f_0
    numtaps = 80
    lp_fir = firwin(numtaps, cutoff_freq, window=('kaiser',8), nyq=nyq)
    filtered_modulated_signal = convolve( modulated_signal, lp_fir )[numtaps/2:len(modulated_signal)+numtaps/2]
    return modulated_signal

def downlink_downconvert(signal, t, f_0):
    
    cutoff_freq = f_0
    nyq = 2 * f_0
    numtaps = 80
    lp_fir = firwin(numtaps, cutoff_freq, window=('kaiser',8), nyq=nyq)
    
    I = -2* convolve( signal * cos(2*pi*f_0*t), lp_fir )[numtaps/2:len(signal)+numtaps/2]
    Q = -2 * convolve( signal * sin(2*pi*f_0*t), lp_fir )[numtaps/2:len(signal)+numtaps/2]
    
    return I + 1j*Q
#@+node:michael.20120305092148.1283: *3* 7.2 Pseudo-random seq gen
def x_1(i):
    x1_init = 1
    while i>30:
        tmp = (x1_init&1 ^ x1_init&8)%2
        x1_init = x1_init>>1 ^ tmp*(2**30)
        i -= 1
    return (x1_init >> i) & 1

def x_2(c_init, i):
    while i>30:
        tmp = (c_init&1 ^ c_init&2 ^ c_init&4 ^ c_init&8)%2
        c_init = c_init>>1 ^ tmp*(2**30)
        i -= 1
    return (c_init >> i) & 1

def c(c_init, i):
    '''
    profile:
def main():
    for i in range(10000):
        c(i,i)

cProfile.runctx('main()', globals(), locals())

   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
    10000    0.059    0.000   91.819    0.009 <ipython console>:1(c)
        1    0.023    0.023   91.843   91.843 <ipython console>:1(main)
    10000   39.318    0.004   39.318    0.004 <ipython console>:1(x_1)
    10000   52.441    0.005   52.441    0.005 <ipython console>:1(x_2)
    '''
    N_C = 1600
    return (x_1(i+N_C) + x_2(c_init,i+N_C)) %2
#@-others

test_enabling_bits = 0b10000

# 01. CSRS baseband IQ time domain signal
if test_enabling_bits & (1<<0):
    csrs_baseband_IQ_time_domain_signal()

# 02. CSRS baseband IQ spectrum
if test_enabling_bits & (1<<1):
    csrs_baseband_IQ_spectrum()

# 03. CSRS Uu signal
if test_enabling_bits & (1<<2):
    CSRS_signal_Uu()

# 04. CSRS received IQ
if test_enabling_bits & (1<<3):
    CSRS_received_IQ()

# 05. CSRS channel estimation in one symbol
if test_enabling_bits & (1<<4):
    test_csrs_channel_estimation_in_one_symbol()
#@-others
#@-leo
