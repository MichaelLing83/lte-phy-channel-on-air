#@+leo-ver=5-thin
#@+node:Michael.20120314113327.1414: * @thin ./Simulation/SSS/sss_gen_detect.py
#@+others
#@+node:Michael.20120314113327.1413: ** source

from scipy.signal import *
from numpy import *
import matplotlib.pyplot as plt

# time scale is in 1 s
T_s = 1.0/30720/1000 # in s

# configuration for SSS
l = 6
N_DL_RB = 110
N_RB_sc = 12
N_DL_CP = 0 # normal DL CP
N_ID_2_tuple = (0,1,2)
delta_f = 15000
subframe = 0
N_ID_cell = 0
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
#@+node:Michael.20120314113327.1415: *3* 01. SSS sequence generation
def SSS_sequence_generation(action='load'):
    import cPickle
    if action=='init':
        sss_dict = dict()   # key is (N_ID_cell, subframe)
        N_ID_cell_list = list()
        for N_ID_1 in range(168):
            for N_ID_2 in range(3):
                N_ID_cell = N_ID_1*3 + N_ID_2
                N_ID_cell_list.append(N_ID_cell)
                for subframe in (0,5):
                    sss_dict[(N_ID_cell,subframe)] = sss_seq(subframe,N_ID_cell)
        f = open('sss_dict.dump','w')
        cPickle.dump(sss_dict,f)
        f.close()
    elif action=='load':
        f = open('sss_dict.dump','r')
        sss_dict = cPickle.load(f)
    return sss_dict
    
#@+node:Michael.20120314113327.1417: *3* 02. SSS baseband IQ time domain signal
def sss_baseband_IQ():
    
    sss_re_array = sss_symbol_array(subframe, N_ID_cell, N_DL_RB, N_RB_sc)
    sss_baseband_IQ = s_p_l(sss_re_array, l, N_DL_RB, N_RB_sc, N_DL_CP, delta_f)
    
    subplot_pos_tupe = (131,132,133)
    title_tuple = ('SSS baseband IQ OFDM signal magnitude','SSS baseband IQ OFDM signal real part','SSS baseband IQ OFDM signal imag part')
    y_label_tuple = ('IQ Magnitude', 'I part', 'Q part')
    func_tuple = (abs, real, imag)
        
    for i in (0,1,2):
        plt.subplot(subplot_pos_tupe[i])
        plt.title(title_tuple[i])
        plt.plot(t*1000, func_tuple[i](sss_baseband_IQ))
        plt.xlabel('Time (ms)')
        plt.ylabel(y_label_tuple[i])
        #plt.axis([-0.01, 0.075, 0, 15])
        plt.legend( ('N_ID_cell=%s'%N_ID_cell,) )
            
    plt.show()
    
#@+node:Michael.20120314113327.1418: *3* 03. SSS baseband IQ spectrum
def sss_baseband_IQ_spectrum():
    
    sss_re_array = sss_symbol_array(subframe, N_ID_cell, N_DL_RB, N_RB_sc)
    sss_baseband_IQ = s_p_l(sss_re_array, l, N_DL_RB, N_RB_sc, N_DL_CP, delta_f)[-1*N:]
    sss_baseband_IQ_fft = fft.fft(sss_baseband_IQ, N)
    
    legend_list = list()
    plt.title('SSS baseband IQ spectrum for N_ID_cell=%s'%N_ID_cell)
    legend_list.append( 'Spectrum magnitude' )
    plt.plot(abs(sss_baseband_IQ_fft), linestyle='-')
    plt.xlabel('n (FFT index)')
    plt.ylabel('Spectrum magnitude')
    plt.legend(legend_list)
    plt.show()
#@+node:michael.20120314211632.1426: *3* 04. SSS baseband IQ correlation
def sss_baseband_IQ_correlation(to_draw=True):
    
    #sss_dict = SSS_sequence_generation()
    
    subframe = 0
    N_ID_cell = 0
    
    sss_re_array = sss_symbol_array(subframe, N_ID_cell, N_DL_RB, N_RB_sc)
    sss_baseband_IQ = s_p_l(sss_re_array, l, N_DL_RB, N_RB_sc, N_DL_CP, delta_f)[-1*N:]
    sss_baseband_IQ_conj = conjugate(sss_baseband_IQ)

    corr_dict = dict()
    cs_list = arange(-1*(N/2), N/2, 1)
    max_dict = dict()

    #print cs_corr.shape, cs_list.shape
    legend_list = list()
    y_offsets = dict()


    corr_dict[(subframe,N_ID_cell)] = array( [0] *N )
    for i in arange(len(cs_list)):
        corr_dict[(subframe,N_ID_cell)][i] = abs(correlate(sss_baseband_IQ, roll(sss_baseband_IQ_conj,cs_list[i]))[0])
    max_dict[(subframe,N_ID_cell)] = find_max(corr_dict[(subframe,N_ID_cell)])
    # normalize the correlation results
    overall_max_y = 0
    for k in max_dict.keys():
        x, y = max_dict[k]
        if y>overall_max_y:
            overall_max_y = y
    overall_max_y = float(overall_max_y)
    for k in max_dict.keys():
        x, y = max_dict[k]
        max_dict[k] = (x, y/overall_max_y)
    for k in corr_dict.keys():
        corr_dict[k] = corr_dict[k]/overall_max_y
    
    if to_draw:
        for subframe,N_ID_cell in max_dict.keys():
            y_offsets[(subframe,N_ID_cell)] = -30
        for subframe,N_ID_cell in corr_dict.keys():
            plt.plot(cs_list, corr_dict[(subframe,N_ID_cell)], marker='+', linestyle='-')
            legend_list.append( 'subframe=%s, N_ID_cell=%s'%(subframe,N_ID_cell) )
            x, y = max_dict[(subframe,N_ID_cell)]
            plt.annotate('Max of subframe=%s, N_ID_cell=%s: %4.4s @ cs=%s'%(subframe,N_ID_cell,y,cs_list[x]), xy=(cs_list[x], y), arrowprops=dict(facecolor='black', shrink=0.15), textcoords='offset points', xytext=(-60, y_offsets[(subframe,N_ID_cell)]))
        plt.title('SSS baseband IQ correlation')
        plt.legend(legend_list)
        plt.xlabel("Cyclic Shift")
        plt.ylabel("Correlation (normalized to peak)")
        plt.show()
    
    return corr_dict
    #plt.savefig('PSS_Uu_signal_inner_products.png', figsize=(1280,800), dpi=200, pad_inches=2)
#@+node:michael.20120314211632.1427: *3* 05. SSS baseband IQ spectrum correlation
def sss_baseband_IQ_spectrum_correlation_ref(ref_subframe, ref_N_ID_cell, to_draw=True):
    
    #sss_dict = SSS_sequence_generation()

    ref_sss_re_array = sss_symbol_array(ref_subframe, ref_N_ID_cell, N_DL_RB, N_RB_sc)
    ref_sss_baseband_IQ = s_p_l(ref_sss_re_array, l, N_DL_RB, N_RB_sc, N_DL_CP, delta_f)[-1*N:]
    ref_sss_baseband_IQ_fft = fft.fft(ref_sss_baseband_IQ, N)
    ref_sss_baseband_IQ_fft_conj = conjugate(ref_sss_baseband_IQ_fft)
    
    subframe_list = (0, 5)
    N_ID_cell_list = (0, 1, 2, 50, 167)
    sss_baseband_IQ_dict = dict()
    for subframe in subframe_list:
        for N_ID_cell in N_ID_cell_list:
            sss_re_array = sss_symbol_array(subframe, N_ID_cell, N_DL_RB, N_RB_sc)
            sss_baseband_IQ_dict[(subframe,N_ID_cell)] = s_p_l(sss_re_array, l, N_DL_RB, N_RB_sc, N_DL_CP, delta_f)[-1*N:]

    corr_dict = dict()
    cs_list = arange(-1*(N/2), N/2, 1)
    max_dict = dict()

    #print cs_corr.shape, cs_list.shape
    legend_list = list()
    y_offsets = dict()

    for subframe in subframe_list:
        for N_ID_cell in N_ID_cell_list:
            corr_dict[(subframe,N_ID_cell)] = array( [0] *N )
            for i in arange(len(cs_list)):
                corr_dict[(subframe,N_ID_cell)][i] = abs(correlate(ref_sss_baseband_IQ_fft_conj, fft.fft(roll(sss_baseband_IQ_dict[(subframe,N_ID_cell)],cs_list[i])))[0])
            max_dict[(subframe,N_ID_cell)] = find_max(corr_dict[(subframe,N_ID_cell)])
    # normalize the correlation results
    overall_max_y = 0
    for k in max_dict.keys():
        x, y = max_dict[k]
        if y>overall_max_y:
            overall_max_y = y
    overall_max_y = float(overall_max_y)
    for k in max_dict.keys():
        x, y = max_dict[k]
        max_dict[k] = (x, y/overall_max_y)
    for k in corr_dict.keys():
        corr_dict[k] = corr_dict[k]/overall_max_y
    
    if to_draw:
        for subframe,N_ID_cell in max_dict.keys():
            y_offsets[(subframe,N_ID_cell)] = 60
        y_offsets[(ref_subframe,ref_N_ID_cell)] = -80
        for subframe,N_ID_cell in corr_dict.keys():
            plt.plot(cs_list, corr_dict[(subframe,N_ID_cell)], marker='+', linestyle='-')
            legend_list.append( 'subframe=%s, N_ID_cell=%s'%(subframe,N_ID_cell) )
            x, y = max_dict[(subframe,N_ID_cell)]
            plt.annotate('Max of subframe=%s, N_ID_cell=%s: %4.4s @ cs=%s'%(subframe,N_ID_cell,y,cs_list[x]), xy=(cs_list[x], y), arrowprops=dict(facecolor='black', shrink=0.15), textcoords='offset points', xytext=(-90, y_offsets[(subframe,N_ID_cell)]))
        plt.title('SSS baseband IQ correlation reference subframe=%s N_ID_cell=%s'%(ref_subframe,ref_N_ID_cell))
        plt.legend(legend_list)
        plt.xlabel("Cyclic Shift")
        plt.ylabel("Correlation (normalized to peak)")
        plt.show()
    
    return corr_dict

def sss_baseband_IQ_spectrum_correlation():
    for ref_subframe in (0, 5):
        for ref_N_ID_cell in (0,):
            sss_baseband_IQ_spectrum_correlation_ref(ref_subframe, ref_N_ID_cell, to_draw=True)
#@+node:michael.20120314211632.1428: *3* 06. SSS Uu signal
def SSS_signal_Uu_ref(ref_subframe, ref_N_ID_cell):
    
    sss_re_array = sss_symbol_array(ref_subframe, ref_N_ID_cell, N_DL_RB, N_RB_sc)
    sss_baseband_IQ = s_p_l(sss_re_array, l, N_DL_RB, N_RB_sc, N_DL_CP, delta_f)
    sss_Uu_signal = downlink_modulate(sss_baseband_IQ, t, f_0)
    
    legend_list = list()
    plt.plot(t*1000, sss_Uu_signal)
    legend_list.append('SSS Uu signal')
    plt.title('SSS Uu signal for subframe=%s N_ID_cell=%s'%(ref_subframe, ref_N_ID_cell))
    plt.xlabel('Time (ms)')
    plt.ylabel('Signal level')
    plt.legend(legend_list)
    #plt.axis( [-0.01, 0.075, -0.1, 14] )
    plt.show()
    #plt.savefig('PSS_signal_Uu_for_N_ID_2=%s.png'%N_ID_2, dpi=300)

def SSS_signal_Uu():
    SSS_signal_Uu_ref(0, 0)
#@+node:michael.20120314211632.1429: *3* 07. SSS received IQ
def SSS_received_IQ_ref(ref_subframe, ref_N_ID_cell):
    
    sss_re_array = sss_symbol_array(ref_subframe, ref_N_ID_cell, N_DL_RB, N_RB_sc)
    sss_baseband_IQ = s_p_l(sss_re_array, l, N_DL_RB, N_RB_sc, N_DL_CP, delta_f)
    sss_Uu_signal = downlink_modulate(sss_baseband_IQ, t, f_0)
    sss_Uu_signal_downconverted = downlink_downconvert(sss_Uu_signal, t, f_0)
    
    legend_list = list()
    plt.plot(t*1000, sss_Uu_signal_downconverted)
    legend_list.append('SSS Uu signal downconverted')
    plt.title('SSS Uu signal downconverted for subframe=%s N_ID_cell=%s'%(ref_subframe, ref_N_ID_cell))
    plt.xlabel('Time (ms)')
    plt.ylabel('Signal level')
    plt.legend(legend_list)
    #plt.axis( [-0.01, 0.075, -0.1, 14] )
    plt.show()
    #plt.savefig('PSS_signal_Uu_for_N_ID_2=%s.png'%N_ID_2, dpi=300)

def SSS_received_IQ():
    SSS_received_IQ_ref(0, 0)
#@+node:michael.20120314211632.1430: *3* 08. SSS received IQ spectrum correlation
def SSS_received_IQ_spectrum_correlation_ref(ref_subframe, ref_N_ID_cell, to_draw=True):

    ref_sss_re_array = sss_symbol_array(ref_subframe, ref_N_ID_cell, N_DL_RB, N_RB_sc)
    ref_sss_baseband_IQ = s_p_l(ref_sss_re_array, l, N_DL_RB, N_RB_sc, N_DL_CP, delta_f)
    ref_sss_Uu_signal = downlink_modulate(ref_sss_baseband_IQ, t, f_0)
    ref_sss_received_IQ = downlink_downconvert(ref_sss_Uu_signal, t, f_0)
    
    ref_sss_baseband_IQ_fft = fft.fft(ref_sss_baseband_IQ[-1*N:], N)
    ref_sss_baseband_IQ_fft_conj = conjugate(ref_sss_baseband_IQ_fft)
    ref_sss_received_IQ_fft = fft.fft(ref_sss_received_IQ[-1*N:], N)
    ref_sss_received_IQ_fft_conj = conjugate(ref_sss_received_IQ_fft)
    
    # must we do coherent detection??
    channel_est = ref_sss_received_IQ_fft/ref_sss_baseband_IQ_fft
    est_ref_sss_received_IQ_fft_conj = conjugate(ref_sss_baseband_IQ_fft * channel_est)
    
    subframe_list = (0,5)
    N_ID_cell_list = (0,1,80,90,167)
    sss_received_IQ_dict = dict()
    for subframe in subframe_list:
        for N_ID_cell in N_ID_cell_list:
            sss_re_array = sss_symbol_array(subframe, N_ID_cell, N_DL_RB, N_RB_sc)
            sss_baseband_IQ = s_p_l(sss_re_array, l, N_DL_RB, N_RB_sc, N_DL_CP, delta_f)
            sss_Uu_signal = downlink_modulate(sss_baseband_IQ, t, f_0)
            sss_received_IQ = downlink_downconvert(sss_Uu_signal, t, f_0)
            sss_received_IQ_dict[(subframe,N_ID_cell)] = sss_received_IQ[-1*N:]

    corr_dict = dict()
    cs_list = arange(-1*(N/2), N/2, 1)
    max_dict = dict()

    #print cs_corr.shape, cs_list.shape
    legend_list = list()
    y_offsets = dict()

    for subframe in subframe_list:
        for N_ID_cell in N_ID_cell_list:
            corr_dict[(subframe,N_ID_cell)] = array( [0] *N )
            for i in arange(len(cs_list)):
                corr_dict[(subframe,N_ID_cell)][i] = abs(correlate(est_ref_sss_received_IQ_fft_conj, fft.fft(roll(sss_received_IQ_dict[(subframe,N_ID_cell)],cs_list[i]), N))[0])
            max_dict[(subframe,N_ID_cell)] = find_max(corr_dict[(subframe,N_ID_cell)])
    # normalize the correlation results
    overall_max_y = 0
    for k in max_dict.keys():
        x, y = max_dict[k]
        if y>overall_max_y:
            overall_max_y = y
    overall_max_y = float(overall_max_y)
    for k in max_dict.keys():
        x, y = max_dict[k]
        max_dict[k] = (x, y/overall_max_y)
    for k in corr_dict.keys():
        corr_dict[k] = corr_dict[k]/overall_max_y
    
    if to_draw:
        for subframe,N_ID_cell in max_dict.keys():
            y_offsets[(subframe,N_ID_cell)] = 60
        y_offsets[(ref_subframe,ref_N_ID_cell)] = -80
        for subframe,N_ID_cell in corr_dict.keys():
            plt.plot(cs_list, corr_dict[(subframe,N_ID_cell)], marker='+', linestyle='-')
            legend_list.append( 'subframe=%s, N_ID_cell=%s'%(subframe,N_ID_cell) )
            x, y = max_dict[(subframe,N_ID_cell)]
            plt.annotate('Max of subframe=%s, N_ID_cell=%s: %4.4s @ cs=%s'%(subframe,N_ID_cell,y,cs_list[x]), xy=(cs_list[x], y), arrowprops=dict(facecolor='black', shrink=0.15), textcoords='offset points', xytext=(-90, y_offsets[(subframe,N_ID_cell)]))
        plt.title('SSS received IQ correlation reference subframe=%s N_ID_cell=%s'%(ref_subframe,ref_N_ID_cell))
        plt.legend(legend_list)
        plt.xlabel("Cyclic Shift")
        plt.ylabel("Correlation (normalized to peak)")
        plt.show()
    
    return corr_dict

def SSS_received_IQ_spectrum_correlation():
    for ref_subframe in (0,):
        for ref_N_ID_cell in (0,):
            SSS_received_IQ_spectrum_correlation_ref(ref_subframe, ref_N_ID_cell, to_draw=True)
#@+node:michael.20120314211632.1431: *3* 09. SSS baseband detect
def SSS_baseband_detect(baseband_IQ_signal, local_t, to_draw=False):
    '''
    SSS_baseband_detect(baseband_IQ_signal, t): index
    return the index of the start of SSS in given baseband IQ signal sequence
    Note: for this function, the parameter t must be of the scale of second, and should not be decimated.
    '''
    #baseband_IQ_signal_conj = conjugate(baseband_IQ_signal)
    #baseband_IQ_signal_I = real(baseband_IQ_signal)
    #baseband_IQ_signal_Q = imag(baseband_IQ_signal)
    
    ref_subframe_list = (0,5)
    ref_N_ID_cell_list = (0,1)
    
    sss_ref_fft_conj_dict = dict()
    for ref_subframe in ref_subframe_list:
        for ref_N_ID_cell in ref_N_ID_cell_list:
            ref_sss_re_array = sss_symbol_array(ref_subframe, ref_N_ID_cell, N_DL_RB, N_RB_sc)
            ref_sss_baseband_IQ = s_p_l(ref_sss_re_array, l, N_DL_RB, N_RB_sc, N_DL_CP, delta_f)
            ref_sss_Uu_signal = downlink_modulate(ref_sss_baseband_IQ, t, f_0)
            ref_sss_received_IQ = downlink_downconvert(ref_sss_Uu_signal, t, f_0)
            ref_sss_baseband_IQ_fft = fft.fft(ref_sss_baseband_IQ[-1*N:], N)
            ref_sss_baseband_IQ_fft_conj = conjugate(ref_sss_baseband_IQ_fft)
            ref_sss_received_IQ_fft = fft.fft(ref_sss_received_IQ[-1*N:], N)
            ref_sss_received_IQ_fft_conj = conjugate(ref_sss_received_IQ_fft)
            # must we do coherent detection??
            channel_est = ref_sss_received_IQ_fft/ref_sss_baseband_IQ_fft
            est_ref_sss_received_IQ_fft_conj = conjugate(ref_sss_baseband_IQ_fft * channel_est)
            sss_ref_fft_conj_dict[(ref_subframe,ref_N_ID_cell)] = est_ref_sss_received_IQ_fft_conj
    
    tmp_t = arange(0, (N_CP_l+N)*T_s, T_s)
    
    legend_list = list()
    corr_dict = dict()
    offset_list = arange(0, len(local_t)-N+1, 1)
    max_dict = dict()
    
    for subframe in ref_subframe_list:
        for N_ID_cell in ref_N_ID_cell_list:
            corr_dict[(subframe,N_ID_cell)] = array( [0.0] * len(offset_list) )
            for offset in offset_list:
                corr_dict[(subframe,N_ID_cell)][offset] = abs(correlate(fft.fft(baseband_IQ_signal[offset:offset+N],N), sss_ref_fft_conj_dict[(subframe,N_ID_cell)])[0])
            max_dict[(subframe,N_ID_cell)] = find_max(corr_dict[(subframe,N_ID_cell)])
            #corr_list[N_ID][offset] = abs(correlate(baseband_IQ_signal_conj[offset:offset+N], pss_baseband_IQ_list[N_ID])[0])
            #corr_list[N_ID][offset] = abs(correlate(baseband_IQ_signal_conj[offset:offset+N], baseband_IQ_signal[144:144+N])[0])
    sframe, n_ID_cell, X, Y = -1, -1, -1, -1
    for subframe,N_ID_cell in max_dict.keys():
        x, y = max_dict[(subframe,N_ID_cell)]
        if y>Y:
            sframe, n_ID_cell, X, Y = subframe, N_ID_cell, x, y
    Y = float(Y)
    for subframe,N_ID_cell in max_dict.keys():
        corr_dict[(subframe,N_ID_cell)] = corr_dict[(subframe,N_ID_cell)]/abs(Y)
    #print Y
    Y = Y/abs(Y)
                
    if to_draw:
        for subframe,N_ID_cell in max_dict.keys():
            plt.plot(1000*local_t[:-1*(N-1)], corr_dict[(subframe,N_ID_cell)], marker='+', linestyle='-')
            legend_list.append( 'subframe=%s, N_ID_cell=%s'%(subframe,N_ID_cell) )
        plt.annotate('Highest peak with subframe=%s N_ID_cell=%s: %4.4s @start_index=%s'%(sframe,n_ID_cell,Y,X), xy=(1000*t[X], Y), arrowprops=dict(facecolor='black', shrink=0.15), textcoords='offset points', xytext=(-20, -30))
        plt.title('SSS baseband detect')
        plt.legend(legend_list)
        plt.xlabel("t (ms)")
        plt.ylabel("Correlation")
        max_t = 1000*local_t[:-1*(N-1)][-1]
        min_t = 1000*local_t[0]
        plt.axis([min_t-max_t*0.01, min(1000*local_t[-1], max_t*1.3), -0.1*Y, 1.1*Y])
        plt.show()
    
    return (sframe, n_ID_cell ,X, Y)

def test_SSS_detect_in_baseband_IQ():
    
    for subframe in (0,):
        for N_ID_cell in (0,):
            sss_re_array = sss_symbol_array(subframe, N_ID_cell, N_DL_RB, N_RB_sc)
            sss_baseband_IQ = s_p_l(sss_re_array, l, N_DL_RB, N_RB_sc, N_DL_CP, delta_f)
            sss_Uu_signal = downlink_modulate(sss_baseband_IQ, t, f_0)
            sss_received_IQ = downlink_downconvert(sss_Uu_signal, t, f_0)

        received_baseband_IQ = array( [0.0+0.0*1j] * len(sss_received_IQ) * 2 )
        received_baseband_IQ[:len(sss_received_IQ)] = sss_received_IQ
        long_t = arange(0, 2*(N_CP_l+N)*T_s, T_s)
        
        SSS_baseband_detect(received_baseband_IQ, long_t, to_draw=True)
#@+node:michael.20120312091134.1404: *3* 6.11.2.1 Sequence generation
def sss_x5(mask, n):
    x_init = 0b10000
    while n>4:
        tmp = 0
        for i in (0,1,2,3,4):
            if (mask>>i)&1 == 1: # feedback for this item is enabled
                tmp += int((x_init>>i)&1)
            tmp = tmp%2
        x_init = x_init>>1 ^ tmp*(2**4)
        n -= 1
    return int((x_init>>n)&1)


def sss_z_(i):
    if i >=0 and i <=30:
        return 1 - 2 * sss_x5(0b10111, i)

def sss_z_1(m, n):
    return sss_z_((n+(m%8))%31)

def sss_c_(i):
    if i>=0 and i<=30:
        return 1 - 2 * sss_x5(0b01001, i)

def sss_c_0(n, N_ID_2):
    return sss_c_((n+N_ID_2)%31)

def sss_c_1(n, N_ID_2):
    return sss_c_((n+N_ID_2+3)%31)

def sss_s_(i):
    return 1 - 2 * sss_x5(0b00101, i)

def sss_s(m, n):
    return sss_s_((n+m)%31)

def sss_d(n, subframe, N_ID_cell):
    N_ID_1 = N_ID_cell/3
    N_ID_2 = N_ID_cell%3
    q_ = N_ID_1/30
    q = (N_ID_1 + q_*(q_+1)/2)/30
    m_ = N_ID_1 + q*(q+1)/2
    m_0 = m_%31
    m_1 = (m_0 + m_/31 + 1)%31
    if n%2==0:
        n_ = n/2
        if subframe == 0:
            result = sss_s(m_0, n_) * sss_c_0(n_, N_ID_2)
        elif subframe == 5:
            result = sss_s(m_1, n_) * sss_c_0(n_, N_ID_2)
    else:
        n_ = (n-1)/2
        if subframe == 0:
            result = sss_s(m_1, n_) * sss_c_1(n_, N_ID_2) * sss_z_1(m_0, n_)
        elif subframe == 5:
            result = sss_s(m_0, n_) * sss_c_1(n_, N_ID_2) * sss_z_1(m_1, n_)
    return result

def sss_seq(subframe, N_ID_cell):
    sss = array([0] * 62)
    for i in range(62):
        sss[i] = sss_d(i, subframe, N_ID_cell)
    return sss
#@+node:Michael.20120314113327.1416: *3* 6.11.2.2 Mapping to REs
def sss_symbol_array(subframe, N_ID_cell, N_DL_RB, N_RB_sc):
    symbol_array = ndarray( shape=(N_DL_RB*N_RB_sc,), dtype=complex128 )
    for i in arange(len(symbol_array)):
        symbol_array[i] = 0
    for n in arange(0, 62):
        k = n-31+N_DL_RB*N_RB_sc/2
        symbol_array[k] = sss_d(n, subframe, N_ID_cell)
    return symbol_array
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
#@-others

test_enabling_bits = 0b111111111

# 01. SSS sequence generation
if test_enabling_bits & (1<<0):
    SSS_sequence_generation()

# 02. SSS baseband signal
if test_enabling_bits & (1<<1):
    sss_baseband_IQ()

# 03. SSS baseband IQ
if test_enabling_bits & (1<<2):
    sss_baseband_IQ_spectrum()

# 04. SSS baseband IQ correlation
if test_enabling_bits & (1<<3):
    sss_baseband_IQ_correlation()

# 05. SSS baseband IQ spectrum correlation
if test_enabling_bits & (1<<4):
    sss_baseband_IQ_spectrum_correlation()

# 06. SSS Uu signal
if test_enabling_bits & (1<<5):
    SSS_signal_Uu()

# 07. SSS received IQ
if test_enabling_bits & (1<<6):
    SSS_received_IQ()

# 08. SSS received IQ spectrum correlation
if test_enabling_bits & (1<<7):
    SSS_received_IQ_spectrum_correlation()

# 09. SSS baseband detect
if test_enabling_bits & (1<<8):
    test_SSS_detect_in_baseband_IQ()
#@-others
#@-leo
