#@+leo-ver=5-thin
#@+node:michael.20120305092148.1319: * @thin ./Simulation/PSS/PSS_gen_detect.py
#@+others
#@+node:michael.20120305092148.1318: ** source code
from scipy.signal import *
from numpy import *
import matplotlib.pyplot as plt

# time scale is in 1 s
T_s = 1.0/30720/1000 # in s

# configuration for PSS
l = 2
N_DL_RB = 110
N_RB_sc = 12
N_DL_CP = 0 # normal DL CP
N_ZC = 63
N_ID_2_tuple = (0,1,2)
delta_f = 15000
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
#@+node:michael.20120305092148.1322: *3* 01. PSS spectrum before OFDM generation
def PSS_spectrum_before_OFDM_generation(to_draw=True):
    
    subplot_pos_tuple = (221,222,223)
    
    for N_ID_2 in N_ID_2_tuple:
        
        plt.subplot(subplot_pos_tuple[N_ID_2])
        legend_list = list()
        pss_freq_symbol_array = pss_symbol_array(N_ID_2, N_DL_RB, N_RB_sc)
        plt.title('PSS spectrum for N_ID_2=%s'%N_ID_2)
        legend_list.append( 'Spectrum magnitude' )
        plt.plot(abs(pss_freq_symbol_array), linestyle='-')
        #legend_list.append( 'Spectrum phase' )
        #plt.plot(abs(pss_freq_symbol_array), linestyle='-')
        plt.xlabel('k (subcarrier index)')
        plt.ylabel('Spectrum magnitude')
        plt.legend(legend_list)
    plt.show()
        #plt.savefig('PSS_spectrum_before_OFDM_gen_for_N_ID_2=%s.png'%N_ID_2)
#@+node:michael.20120305092148.1320: *3* 02. PSS correlation in freq domain before OFDM gen
def PSS_corr_in_freq_domain_before_OFDM_gen(to_draw=True):
    
    zc_seq_d = dict()
    for N_ID_2 in N_ID_2_tuple:
        zc_seq_d[N_ID_2] = array([0]*N_ZC, dtype=complex128)
        for i in arange(N_ZC-1):
            if i<=30:
                zc_seq_d[N_ID_2][i] = pss_d(i, N_ID_2)
            else:
                zc_seq_d[N_ID_2][i] = pss_d(i, N_ID_2)
    
    legend_list = list()
    corr_dict = dict()
    max_dict = dict()
    cs_list = arange(-1*(N_ZC/2), N_ZC/2+1, 1)
    y_offsets = dict()
    y_offsets[(0,0)] = -20
    y_offsets[(0,1)] = 10
    y_offsets[(0,2)] = 20
    y_offsets[(1,1)] = -40
    y_offsets[(1,2)] = 20
    y_offsets[(2,2)] = -60
    
    overall_max_y = 0
    #print cs_corr.shape, cs_list.shape
    for p in N_ID_2_tuple:
        for q in N_ID_2_tuple:
            if p<=q:
                corr_dict[(p,q)] = array([0]*N_ZC)
                for i in arange(len(cs_list)):
                    corr_dict[(p,q)][i] = abs(correlate(zc_seq_d[p],roll(zc_seq_d[q],cs_list[i]))[0])
                max_dict[(p,q)] = find_max(corr_dict[(p,q)])
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
        for p,q in corr_dict.keys():
            plt.plot(cs_list, corr_dict[(p,q)], marker='+', linestyle='-')
            legend_list.append( 'N_ID_2 %s vs. %s'%(p,q) )
            x, y = max_dict[(p,q)]
            plt.annotate('Max of N_ID_2 %svs%s =%4.4s @ cs=%s'%(p,q,y,cs_list[x]), xy=(cs_list[x], y), arrowprops=dict(facecolor='black', shrink=0.15), textcoords='offset points', xytext=(60, y_offsets[(p,q)]))
        plt.title('PSS correlation in freq domain before OFDM generation')
        plt.legend(legend_list)
        plt.xlabel("Cyclic Shift")
        plt.ylabel("Correlation (normalized to peak)")
        plt.show()
    
    return corr_dict

    #plt.savefig('PSS_corr_in_freq_domain_before_OFDM_gen.png', dpi=300)
#@+node:michael.20120305092148.1321: *3* 03. PSS baseband IQ signal in time domain
def PSS_baseband_IQ_signal_in_time_domain():
    
    subplot_pos_tupe = (    (331,332,333),
                                    (334,335,336),
                                    (337,338,339)
                                )
    title_tuple = ('PSS baseband IQ OFDM signal magnitude','PSS baseband IQ OFDM signal real part','PSS baseband IQ OFDM signal imag part')
    y_label_tuple = ('IQ Magnitude', 'I part', 'Q part')
    func_tuple = (abs, real, imag)
    pss_baseband_symbol_list = [0]*3
    for N_ID_2 in N_ID_2_tuple:
        pss_freq_symbol_array = pss_symbol_array(N_ID_2, N_DL_RB, N_RB_sc)
        pss_baseband_symbol_list[N_ID_2] = ofdm_baseband_IQ_signal_generate(pss_freq_symbol_array, l, N_DL_RB, N_RB_sc, N_DL_CP, delta_f)
        
    for i in (0,1,2):
        for N_ID_2 in N_ID_2_tuple:
            plt.subplot(131+N_ID_2)
            if N_ID_2 == 1:
                plt.title(title_tuple[i])
            plt.plot(t*1000, func_tuple[i](pss_baseband_symbol_list[N_ID_2]))
            plt.xlabel('Time (ms)')
            plt.ylabel(y_label_tuple[i])
            plt.axis([-0.01, 0.075, 0, 15])
            plt.legend( ('N_ID_2=%s'%N_ID_2,) )
            
        plt.show()
#@+node:michael.20120310203114.1395: *3* 04. PSS baseband IQ spectrum
def PSS_baseband_IQ_spectrum(to_draw=True):
    
    subplot_pos_tuple = (221,222,223)
    
    for N_ID_2 in N_ID_2_tuple:
        
        plt.subplot(subplot_pos_tuple[N_ID_2])
        legend_list = list()
        pss_freq_symbol_array = pss_symbol_array(N_ID_2, N_DL_RB, N_RB_sc)
        #pss_freq_symbol_array_ext = array([0+0*1j]*N)
        #pss_freq_symbol_array_ext[N/2-31:N/2] = pss_freq_symbol_array[len(pss_freq_symbol_array)/2-31:len(pss_freq_symbol_array)/2]
        #pss_freq_symbol_array_ext[N/2:N/2+31] = pss_freq_symbol_array[len(pss_freq_symbol_array)/2:len(pss_freq_symbol_array)/2+31]
        #print pss_freq_symbol_array_ext[N/2]
        #pss_ifft = fft.ifft(pss_freq_symbol_array, N)
        #pss_fft = fft.fft(pss_ifft, N)
        pss_baseband_IQ = ofdm_baseband_IQ_signal_generate(pss_freq_symbol_array, l, N_DL_RB, N_RB_sc, N_DL_CP, delta_f)[-1*N:]
        pss_baseband_IQ_fft = fft.fft(pss_baseband_IQ, N)
        plt.title('PSS baseband IQ spectrum for N_ID_2=%s'%N_ID_2)
        legend_list.append( 'Spectrum magnitude' )
        plt.plot(abs(pss_baseband_IQ_fft), linestyle='-')
        #legend_list.append( 'Spectrum phase' )
        #plt.plot(abs(pss_freq_symbol_array), linestyle='-')
        plt.xlabel('n (FFT index)')
        plt.ylabel('Spectrum magnitude')
        plt.legend(legend_list)
    plt.show()
        #plt.savefig('PSS_spectrum_before_OFDM_gen_for_N_ID_2=%s.png'%N_ID_2)
#@+node:michael.20120312091134.1399: *3* 05. PSS baseband IQ spectrum correlation
def PSS_baseband_IQ_spectrum_correlation(to_draw=True):
    
    pss_iq_sig_list = [0]*3
    pss_baseband_IQ_FFT_list = [0] * 3
    corr_dict = dict()
    cs_list = arange(-1*(N/2), N/2, 1)
    max_dict = dict()
    
    for N_ID_2 in N_ID_2_tuple:
        
        pss_freq_symbol_array = pss_symbol_array(N_ID_2, N_DL_RB, N_RB_sc)
        pss_baseband_IQ = ofdm_baseband_IQ_signal_generate(pss_freq_symbol_array, l, N_DL_RB, N_RB_sc, N_DL_CP, delta_f)[-1*N:]
        pss_baseband_IQ_FFT_list[N_ID_2] = fft.fft(pss_baseband_IQ, N)

    #print cs_corr.shape, cs_list.shape
    legend_list = list()
    y_offsets = dict()
    y_offsets[(0,0)] = -20
    y_offsets[(0,1)] = 10
    y_offsets[(0,2)] = 20
    y_offsets[(1,1)] = -35
    y_offsets[(1,2)] = 40
    y_offsets[(2,2)] = -50
    for p in N_ID_2_tuple:
        for q in N_ID_2_tuple:
            if p<=q:
                corr_dict[(p,q)] = array( [0] *N )
                for i in arange(len(cs_list)):
                    corr_dict[(p,q)][i] = abs(correlate(pss_baseband_IQ_FFT_list[p], roll(pss_baseband_IQ_FFT_list[q],cs_list[i]))[0])
                max_dict[(p,q)] = find_max(corr_dict[(p,q)])
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
        for p,q in corr_dict.keys():
            plt.plot(cs_list, corr_dict[(p,q)], marker='+', linestyle='-')
            legend_list.append( 'N_ID_2 %s vs. %s'%(p,q) )
            x, y = max_dict[(p,q)]
            plt.annotate('Max of N_ID_2 %svs%s =%4.4s @ cs=%s'%(p,q,y,cs_list[x]), xy=(cs_list[x], y), arrowprops=dict(facecolor='black', shrink=0.15), textcoords='offset points', xytext=(60, y_offsets[(p,q)]))
        plt.title('PSS baseband IQ spectrum correlation after OFDM generation')
        plt.legend(legend_list)
        plt.xlabel("Cyclic Shift")
        plt.ylabel("Correlation (normalized to peak)")
        plt.show()
    
    return corr_dict
    #plt.savefig('PSS_Uu_signal_inner_products.png', figsize=(1280,800), dpi=200, pad_inches=2)

#@+node:michael.20120305092148.1323: *3* 06. PSS Uu signal
def PSS_signal_Uu():
    
    for N_ID_2 in N_ID_2_tuple:
        
        pss_freq_symbol_array = pss_symbol_array(N_ID_2, N_DL_RB, N_RB_sc)
        pss_baseband_symbol = ofdm_baseband_IQ_signal_generate(pss_freq_symbol_array, l, N_DL_RB, N_RB_sc, N_DL_CP, delta_f)
        pss_uu_sig = downlink_modulate(pss_baseband_symbol, t, f_0)
        
        #plt.plot(t, symbol_array)
        plt.cla()
        legend_list = list()
        plt.plot(t*1000, pss_uu_sig)
        legend_list.append( 'PSS signal @Uu for N_ID_2=%s'%(N_ID_2) )
        plt.title('PSS signal @Uu for N_ID_2=%s'%N_ID_2)
        plt.xlabel('Time (ms)')
        plt.ylabel('Signal level')
        plt.legend(legend_list)
        #plt.axis( [-0.01, 0.075, -0.1, 14] )
        plt.show()
        #plt.savefig('PSS_signal_Uu_for_N_ID_2=%s.png'%N_ID_2, dpi=300)
#@+node:michael.20120309091906.1387: *3* 07. PSS Uu signal downconversion
def PSS_signal_Uu_downconversion():
    
    
    for N_ID_2 in N_ID_2_tuple:
        
        pss_freq_symbol_array = pss_symbol_array(N_ID_2, N_DL_RB, N_RB_sc)
        pss_baseband_symbol = ofdm_baseband_IQ_signal_generate(pss_freq_symbol_array, l, N_DL_RB, N_RB_sc, N_DL_CP, delta_f)
        pss_uu_sig = downlink_modulate(pss_baseband_symbol, t, f_0)
        pss_uu_sig_down = downlink_downconvert(pss_uu_sig, t, f_0)
        
        #plt.plot(t, symbol_array)
        legend_list = list()
        plt.plot(t*1000, real(pss_uu_sig_down))
        legend_list.append( 'PSS Uu signal downconverted for N_ID_2=%s'%(N_ID_2) )
        plt.plot(t*1000, real(pss_baseband_symbol))
        legend_list.append( 'PSS baseband IQ for N_ID_2=%s'%(N_ID_2) )
        plt.title('PSS signal Uu downconverted for N_ID_2=%s'%N_ID_2)
        plt.xlabel('Time (ms)')
        plt.ylabel('IQ signal real part')
        plt.legend(legend_list)
        #plt.axis( [-0.01, 0.075, -0.1, 14] )
        plt.show()
        #plt.savefig('PSS_signal_Uu_for_N_ID_2=%s.png'%N_ID_2, dpi=300)
#@+node:michael.20120312091134.1401: *3* 08. PSS Uu signal downconverted correlation
def PSS_Uu_signal_downconverted_correlation():
    for correlation_type in ('IQ', 'I+Q', 'I', 'Q'):
        PSS_Uu_signal_downconverted_correlation_IQ(correlation_type)

def PSS_Uu_signal_downconverted_correlation_IQ(correlation_type='I+Q', to_draw=True):
    
    pss_Uu_signal_downconverted_IQ_list = [0]*3
    pss_Uu_signal_downconverted_I_list = [0]*3
    pss_Uu_signal_downconverted_Q_list = [0]*3
    pss_baseband_IQ_list = [0]*3
    pss_baseband_I_list = [0]*3
    pss_baseband_Q_list = [0]*3
    
    for N_ID_2 in N_ID_2_tuple:
        pss_freq_symbol_array = pss_symbol_array(N_ID_2, N_DL_RB, N_RB_sc)
        pss_baseband_IQ = ofdm_baseband_IQ_signal_generate(pss_freq_symbol_array, l, N_DL_RB, N_RB_sc, N_DL_CP, delta_f)
        pss_baseband_IQ_list[N_ID_2] = pss_baseband_IQ[-1*N:]
        pss_baseband_I_list[N_ID_2] = real(pss_baseband_IQ)[-1*N:]
        pss_baseband_Q_list[N_ID_2] = imag(pss_baseband_IQ)[-1*N:]
        pss_uu_sig = downlink_modulate(pss_baseband_IQ, t, f_0)
        pss_received_IQ = downlink_downconvert(pss_uu_sig, t, f_0)[-1*N:]
        pss_Uu_signal_downconverted_IQ_list[N_ID_2] = pss_received_IQ
        pss_Uu_signal_downconverted_I_list[N_ID_2] = real(pss_received_IQ)
        pss_Uu_signal_downconverted_Q_list[N_ID_2] = imag(pss_received_IQ)
    
    legend_list = list()
    corr_dict = dict()
    max_dict = dict()
    min_dict = dict()
    cs_list = arange(-1*(N/2), N/2, 1)
    #print cs_corr.shape, cs_list.shape
    legend_list = list()
    y_offsets = dict()
    y_offsets[(0,0)] = -20
    y_offsets[(0,1)] = 10
    y_offsets[(0,2)] = 20
    y_offsets[(1,1)] = -40
    y_offsets[(1,2)] = -80
    y_offsets[(2,2)] = -60
    for p in N_ID_2_tuple:
        for q in N_ID_2_tuple:
            if p<=q:
                corr_dict[(p,q)] = array([0]*N)
                for i in arange(len(cs_list)):
                    if correlation_type=='I+Q':
                        corr_dict[(p,q)][i] = abs(correlate(pss_baseband_I_list[p], roll(pss_Uu_signal_downconverted_I_list[q],cs_list[i]))[0]) + abs(correlate(pss_baseband_Q_list[p], roll(pss_Uu_signal_downconverted_Q_list[q],cs_list[i]))[0])
                    elif correlation_type=='I':
                        corr_dict[(p,q)][i] = abs(correlate(pss_baseband_I_list[p], roll(pss_Uu_signal_downconverted_I_list[q],cs_list[i]))[0])
                    elif correlation_type=='Q':
                        corr_dict[(p,q)][i] = correlate(pss_baseband_Q_list[p], roll(pss_Uu_signal_downconverted_Q_list[q],cs_list[i]))[0]
                    elif correlation_type=='IQ':
                        corr_dict[(p,q)][i] = abs(correlate(pss_baseband_IQ_list[p], roll(conjugate(pss_Uu_signal_downconverted_IQ_list[q]),cs_list[i]))[0])
                max_dict[(p,q)] = find_max(corr_dict[(p,q)])
                min_dict[(p,q)] = find_min(corr_dict[(p,q)])
                
    overall_max_y = 0
    for k in max_dict.keys():
        x, y = max_dict[k]
        if y>overall_max_y:
            overall_max_y = y
    for k in min_dict.keys():
        x, y = min_dict[k]
        if abs(y)>overall_max_y:
            overall_max_y = abs(y)
    overall_max_y = float(overall_max_y)
    for k in max_dict.keys():
        x, y = max_dict[k]
        max_dict[k] = (x, y/overall_max_y)
    for k in min_dict.keys():
        x, y = min_dict[k]
        min_dict[k] = (x, y/overall_max_y)
    for k in corr_dict.keys():
        corr_dict[k] = corr_dict[k]/overall_max_y
    if to_draw:
        for p,q in corr_dict.keys():
            plt.plot(cs_list, corr_dict[(p,q)], marker='+', linestyle='-')
            legend_list.append( 'N_ID_2 %s vs. %s'%(p,q) )
            x, y = max_dict[(p,q)]
            plt.annotate('Max of N_ID_2 %svs%s =%4.4s @ cs=%s'%(p,q,y,cs_list[x]), xy=(cs_list[x], y), arrowprops=dict(facecolor='black', shrink=0.15), textcoords='offset points', xytext=(60, y_offsets[(p,q)]))
            if correlation_type=='Q':
                x, y = min_dict[(p,q)]
                plt.annotate('Min of N_ID_2 %svs%s =%4.4s @ cs=%s'%(p,q,y,cs_list[x]), xy=(cs_list[x], y), arrowprops=dict(facecolor='black', shrink=0.15), textcoords='offset points', xytext=(30, -1*y_offsets[(p,q)]))
                #print x,y
        plt.title('PSS Uu signal downconverted correlation, type: %s'%correlation_type)
        plt.legend(legend_list)
        plt.xlabel("Cyclic Shift")
        plt.ylabel("Correlation (normalized to peak)")
        plt.show()
    
    return corr_dict
#@+node:michael.20120310203114.1394: *3* 09. PSS received IQ spectrum
def PSS_received_IQ_spectrum(to_draw=True):
    
    subplot_pos_tuple = (221,222,223)
    
    for N_ID_2 in N_ID_2_tuple:
        
        plt.subplot(subplot_pos_tuple[N_ID_2])
        legend_list = list()
        pss_freq_symbol_array = pss_symbol_array(N_ID_2, N_DL_RB, N_RB_sc)
        pss_baseband_symbol = ofdm_baseband_IQ_signal_generate(pss_freq_symbol_array, l, N_DL_RB, N_RB_sc, N_DL_CP, delta_f)
        pss_uu_sig = downlink_modulate(pss_baseband_symbol, t, f_0)
        pss_uu_sig_down = downlink_downconvert(pss_uu_sig, t, f_0)[-1*N:]
        pss_uu_sig_down_fft = fft.fft(pss_uu_sig_down, N)
        plt.title('Received PSS IQ spectrum for N_ID_2=%s'%N_ID_2)
        legend_list.append( 'Spectrum magnitude' )
        plt.plot(abs(pss_uu_sig_down_fft), linestyle='-')
        #legend_list.append( 'Spectrum phase' )
        #plt.plot(abs(pss_freq_symbol_array), linestyle='-')
        plt.xlabel('n (FFT index)')
        plt.ylabel('Spectrum magnitude')
        plt.legend(legend_list)
    plt.show()
        #plt.savefig('PSS_spectrum_before_OFDM_gen_for_N_ID_2=%s.png'%N_ID_2)
#@+node:michael.20120305092148.1326: *3* 10. PSS Uu signal downconverted decimated to 1/16 correlation
def PSS_Uu_signal_downconverted_decimated_1_16_correlation(to_draw=True):
    
    pss_Uu_signal_downconverted_list = [0]*3
    
    for N_ID_2 in N_ID_2_tuple:
        pss_freq_symbol_array = pss_symbol_array(N_ID_2, N_DL_RB, N_RB_sc)
        pss_baseband_IQ = ofdm_baseband_IQ_signal_generate(pss_freq_symbol_array, l, N_DL_RB, N_RB_sc, N_DL_CP, delta_f)
        pss_uu_sig = downlink_modulate(pss_baseband_IQ, t, f_0)
        pss_Uu_signal_downconverted_list[N_ID_2] = downlink_downconvert(pss_uu_sig, t, f_0)[-1*N:]
    
    N_16 = N/16
    pss_Uu_signal_downconverted_dec_list = [0]*3
    for N_ID_2 in N_ID_2_tuple:
        pss_Uu_signal_downconverted_dec_list[N_ID_2] = array([0.0+0.0*1j]*N_16)
        for i in arange(N_16):
            pss_Uu_signal_downconverted_dec_list[N_ID_2][i] = pss_Uu_signal_downconverted_list[N_ID_2][16*i]
    legend_list = list()
    corr_dict = dict()
    max_dict = dict()
    cs_list = arange(-1*(N_16/2), N_16/2, 1)
    #print cs_corr.shape, cs_list.shape
    legend_list = list()
    y_offsets = dict()
    y_offsets[(0,0)] = -20
    y_offsets[(0,1)] = 10
    y_offsets[(0,2)] = 20
    y_offsets[(1,1)] = -10
    y_offsets[(1,2)] = 10
    y_offsets[(2,2)] = -60
    for p in N_ID_2_tuple:
        for q in N_ID_2_tuple:
            if p<=q:
                corr_dict[(p,q)] = array([0]*N_16)
                for i in arange(len(cs_list)):
                    corr_dict[(p,q)][i] = abs(correlate(pss_Uu_signal_downconverted_dec_list[p], roll(pss_Uu_signal_downconverted_dec_list[q],cs_list[i]))[0])
                max_dict[(p,q)] = find_max(corr_dict[(p,q)])
                
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
        for p,q in corr_dict.keys():
            plt.plot(cs_list, corr_dict[(p,q)], marker='+', linestyle='-')
            legend_list.append( 'N_ID_2 %s vs. %s'%(p,q) )
            x, y = max_dict[(p,q)]
            plt.annotate('Max of N_ID_2 %svs%s =%4.4s @ cs=%s'%(p,q,y,cs_list[x]), xy=(cs_list[x], y), arrowprops=dict(facecolor='black', shrink=0.15), textcoords='offset points', xytext=(60, y_offsets[(p,q)]))
        plt.title('PSS Uu signal downconverted correlation')
        plt.legend(legend_list)
        plt.xlabel("Cyclic Shift")
        plt.ylabel("Correlation (normalized to peak)")
        plt.show()
    
    return corr_dict
#@+node:michael.20120312091134.1402: *3* 11. PSS detect
def PSS_baseband_detect(baseband_IQ_signal, local_t, to_draw=False):
    '''
    PSS_baseband_detect(baseband_IQ_signal, t): index
    return the index of the start of PSS in given baseband IQ signal sequence
    Note: for this function, the parameter t must be of the scale of second, and should not be decimated.
    '''
    baseband_IQ_signal_conj = conjugate(baseband_IQ_signal)
    baseband_IQ_signal_I = real(baseband_IQ_signal)
    baseband_IQ_signal_Q = imag(baseband_IQ_signal)
    
    pss_baseband_IQ_list = [0] * 3
    pss_baseband_I_list = [0]*3
    pss_baseband_Q_list = [0]*3
    
    tmp_t = arange(0, (N_CP_l+N)*T_s, T_s)
    #N_ID_2_tuple = (0,1)
    for N_ID_2 in N_ID_2_tuple:
        
        pss_seq = pss_symbol_array(N_ID_2, N_DL_RB, N_RB_sc)
        pss_baseband_IQ = ofdm_baseband_IQ_signal_generate(pss_seq, l, N_DL_RB, N_RB_sc, N_DL_CP, delta_f)
        pss_baseband_IQ_list[N_ID_2] = pss_baseband_IQ[-1*N:]
        pss_baseband_I_list[N_ID_2] = real(pss_baseband_IQ[-1*N:])
        pss_baseband_Q_list[N_ID_2] = imag(pss_baseband_IQ[-1*N:])
    
    legend_list = list()
    corr_list = [0]*3
    offset_list = arange(0, len(local_t)-N+1, 1)
    max_list = [0]*3
    for i in N_ID_2_tuple:
        corr_list[i] = array( [0.0] * len(offset_list) )
    
    legend_list = list()
    y_offsets = (-80, -80, -80)
    
    for offset in offset_list:
        for N_ID in N_ID_2_tuple:
            #corr_list[N_ID][offset] = correlate(baseband_IQ_signal_Q[offset:offset+N], pss_baseband_Q_list[N_ID])[0]
            corr_list[N_ID][offset] = abs(correlate(baseband_IQ_signal_conj[offset:offset+N], pss_baseband_IQ_list[N_ID])[0])
            #corr_list[N_ID][offset] = abs(correlate(baseband_IQ_signal_conj[offset:offset+N], baseband_IQ_signal[144:144+N])[0])
    n_ID_2, X, Y = -1, 0, 0
    for N_ID_2 in N_ID_2_tuple:
        x, y = find_max(corr_list[N_ID_2])
        if y>abs(Y):
            n_ID_2, X, Y = N_ID_2, x, y
        x, y = find_min(corr_list[N_ID_2])
        if abs(y)>abs(Y):
            n_ID_2, X, Y = N_ID_2, x, y
    Y = float(Y)
    for N_ID_2 in N_ID_2_tuple:
        corr_list[N_ID_2] = corr_list[N_ID_2]/abs(Y)
    Y = Y/abs(Y)
    if n_ID_2==1 and Y<0:
        n_ID_2 = 2
    
    if to_draw:
        for N_ID_2 in N_ID_2_tuple:
            plt.plot(1000*local_t[:-1*(N-1)], corr_list[N_ID_2], marker='+', linestyle='-')
            legend_list.append( 'N_ID_2=%s'%N_ID_2 )
        plt.annotate('Highest peak with N_ID_2=%s: %4.4s @start_index=%s'%(n_ID_2,Y,X), xy=(1000*t[X], Y), arrowprops=dict(facecolor='black', shrink=0.15), textcoords='offset points', xytext=(-20, y_offsets[N_ID_2]))
        plt.title('PSS baseband detect')
        plt.legend(legend_list)
        plt.xlabel("t (ms)")
        plt.ylabel("Correlation")
        max_t = 1000*local_t[:-1*(N-1)][-1]
        min_t = 1000*local_t[0]
        plt.axis([min_t-max_t*0.01, min(1000*local_t[-1], max_t*1.3), -0.1*Y, 1.1*Y])
        plt.show()
    
    return (n_ID_2,X,Y)

def test_PSS_detect_in_baseband_IQ():
    
    for n_ID_2 in N_ID_2_tuple:

        pss_sequence = pss_symbol_array(n_ID_2, N_DL_RB, N_RB_sc)
        pss_baseband_IQ = ofdm_baseband_IQ_signal_generate(pss_sequence, l, N_DL_RB, N_RB_sc, N_DL_CP, delta_f)
        pss_Uu_sig = downlink_modulate(pss_baseband_IQ, t, f_0)
        pss_received_IQ = downlink_downconvert(pss_Uu_sig, t, f_0)
        received_baseband_IQ = array( [0.0+0.0*1j] * (len(pss_received_IQ) * 2) )
        received_baseband_IQ[:len(pss_received_IQ)] = pss_received_IQ
        long_t = arange(0, 2*(N_CP_l+N)*T_s, T_s)
        
        PSS_baseband_detect(received_baseband_IQ, long_t, to_draw=True)
        #PSS_baseband_detect(pss_received_IQ, t, to_draw=True)
#@+node:Michael.20120316092234.1457: *3* 12. Channel estimation using PSS
def channel_estimation_using_PSS(baseband_IQ_array, N_ID_2, l, N_DL_RB, N_RB_sc, delta_f, to_draw=False):
    '''
    Note: len(baseband_IQ_array)==N must be True.
    '''
    pss_freq_symbol_array = pss_symbol_array(N_ID_2, N_DL_RB, N_RB_sc)
    #pss_baseband_symbol = ofdm_baseband_IQ_signal_generate(pss_freq_symbol_array, l, N_DL_RB, N_RB_sc, N_DL_CP, delta_f)
    pss_seq_send = get_pss_seq_from_RE_symbol_array(pss_freq_symbol_array, N_RB_sc)
    
    pss_seq_received = get_pss_seq_from_RE_symbol_array(baseband_IQ_array, N_RB_sc)
    pss_channel_estimation = pss_seq_received / pss_seq_send
    
    #pss_uu_sig = downlink_modulate(pss_baseband_symbol, t, f_0)
    #pss_uu_sig_down = downlink_downconvert(pss_uu_sig, t, f_0)[-1*N:]
    #pss_received_RE_IQ_array = ofdm_baseband_IQ_to_RE_IQ_array(pss_uu_sig_down, N_DL_RB, N_RB_sc, delta_f=15000)
    if to_draw:
        subplot_pos_tuple = (121,122)
        #plt.subplot(subplot_pos_tuple[N_ID_2])
        legend_list = list()
        
        plt.subplot(121)
        plt.title('Channel estimation magnitude')
        plt.xlabel('n (subcarrier index)')
        plt.ylabel('Channel est. magnitude')
        plt.plot(abs(pss_channel_estimation))
        #plt.show()
        
        plt.subplot(122)
        plt.title('Channel estimation phase')
        plt.xlabel('n (subcarrier index)')
        plt.ylabel('Channel est. phase')
        plt.plot(angle(pss_channel_estimation))
        
        plt.show()
    
    return pss_channel_estimation
    
def test_channel_estimation_using_PSS():
    N_ID_2 = 0
    pss_freq_symbol_array = pss_symbol_array(N_ID_2, N_DL_RB, N_RB_sc)
    pss_baseband_symbol = ofdm_baseband_IQ_signal_generate(pss_freq_symbol_array, l, N_DL_RB, N_RB_sc, N_DL_CP, delta_f)
    #pss_seq_send = get_pss_seq_from_RE_symbol_array(pss_freq_symbol_array, N_RB_sc)
    
    #pss_seq_received = get_pss_seq_from_RE_symbol_array(pss_received_RE_IQ_array, N_RB_sc)
    #pss_channel_estimation = pss_seq_received / pss_seq_send
    
    pss_uu_sig = downlink_modulate(pss_baseband_symbol, t, f_0)
    pss_received_IQ_array = downlink_downconvert(pss_uu_sig, t, f_0)[-1*N:]
    #pss_received_RE_IQ_array = ofdm_baseband_IQ_to_RE_IQ_array(pss_uu_sig_down, N_DL_RB, N_RB_sc, delta_f=15000)
    channel_estimation_using_PSS(pss_received_IQ_array, N_ID_2, l, N_DL_RB, N_RB_sc, delta_f, to_draw=True)
#@+node:michael.20120305092148.1303: *3* 6.11.1.1 seq gen
def pss_d(n, N_ID_2):
    u = (25, 29, 34)[N_ID_2]
    d_n = 0
    if n>=0 and n<=30:
        d_n = exp(-1j*pi*u*n*(n+1)/63)
    elif n>=31 and n<=61:
        d_n = exp(-1j*pi*u*(n+1)*(n+2)/63)
    return d_n
#@+node:michael.20120305092148.1305: *3* 6.11.1.2 mapping to REs
def pss_symbol_array(N_ID_2, N_DL_RB, N_RB_sc):
    symbol_array = ndarray( shape=(N_DL_RB*N_RB_sc,), dtype=complexfloating )
    for i in arange(len(symbol_array)):
        symbol_array[i] = 0
    for n in arange(0, 62):
        k = n-31+N_DL_RB*N_RB_sc/2
        symbol_array[k] = pss_d(n, N_ID_2)
    return symbol_array

def pss_k_range(N_DL_RB, N_RB_sc):
    
    start_index = 0-31+N_DL_RB*N_RB_sc/2
    end_index = 61-31+N_DL_RB*N_RB_sc/2
    
    return (start_index, end_index)

def get_pss_seq_from_RE_symbol_array(re_symbol_array, N_RB_sc):
    '''
    Note: len(re_symbol_array)==N_DL_RB*N_RB_sc must be True!
    '''
    pss_seq_received = array([0.0 + 0.0 * 1j] * (6 * N_RB_sc))
    pss_start_k, pss_end_k = pss_k_range(N_DL_RB, N_RB_sc)
    tmp_index = pss_start_k
    for i in arange(5, 6*N_RB_sc):
        pss_seq_received[i] = re_symbol_array[tmp_index]
        tmp_index += 1
    return pss_seq_received
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
#@+node:michael.20120305092148.1330: *3* test
def test(to_draw=True):
    
    subplot_pos_tuple = (221,222,223)
    
    for N_ID_2 in N_ID_2_tuple:
        
        plt.subplot(subplot_pos_tuple[N_ID_2])
        legend_list = list()
        all_subcarrier_array = array([0] * (N_DL_RB * N_RB_sc+1))
        pss_freq_symbol_array = pss_symbol_array(N_ID_2, N_DL_RB, N_RB_sc)
        all_subcarrier_array[:N_DL_RB*N_RB_sc/2]
        pss_baseband_IQ = ofdm_baseband_IQ_signal_generate(pss_freq_symbol_array, l, N_DL_RB, N_RB_sc, N_DL_CP, delta_f)[-1*N:]
        pss_baseband_IQ_fft = fft.fft(pss_baseband_IQ, N)
        plt.title('PSS baseband IQ spectrum for N_ID_2=%s'%N_ID_2)
        legend_list.append( 'Spectrum magnitude' )
        plt.plot(abs(pss_baseband_IQ_fft), linestyle='-')
        #legend_list.append( 'Spectrum phase' )
        #plt.plot(abs(pss_freq_symbol_array), linestyle='-')
        plt.xlabel('n (FFT index)')
        plt.ylabel('Spectrum magnitude')
        plt.legend(legend_list)
    plt.show()
        #plt.savefig('PSS_spectrum_before_OFDM_gen_for_N_ID_2=%s.png'%N_ID_2)
#@+node:michael.20120312091134.1398: *3* Trash Can
#@+node:michael.20120305092148.1327: *4* 08. PSS received IQ FFT
def PSS_received_IQ_FFT():
    
    func_tuple = (abs, real, imag)
    legend_tuple = ('Absolute', 'Real part', 'Image part')
    
    
    for N_ID_2 in N_ID_2_tuple:
        
        pss_freq_symbol_array = pss_symbol_array(N_ID_2, N_DL_RB, N_RB_sc)
        pss_baseband_symbol = ofdm_baseband_IQ_signal_generate(pss_freq_symbol_array, l, N_DL_RB, N_RB_sc, N_DL_CP, delta_f)
        pss_uu_sig = downlink_modulate(pss_baseband_symbol, t, f_0)
        pss_uu_sig_down = downlink_downconvert(pss_uu_sig, t, f_0)
        pss_uu_sig_down_no_cp = pss_uu_sig_down[-1*N:]
        pss_uu_fft = fft.fft( pss_uu_sig_down_no_cp, N )
        
        for i in (0,1,2):
            plt.subplot(131+i)
            if i==1:
                plt.title('PSS received IQ FFT for N_ID_2=%s'%N_ID_2)
            legend_list = list()
            plt.plot(func_tuple[i](pss_uu_fft))
            legend_list.append( legend_tuple[i] )
            plt.xlabel('N')
            #plt.ylabel('Spectrum absolute value')
            plt.legend(legend_list)
            #plt.axis( [-0.01, 0.075, -0.1, 14] )
        plt.show()
            #plt.savefig('PSS_signal_Uu_for_N_ID_2=%s.png'%N_ID_2, dpi=300)
#@+node:michael.20120312091134.1400: *4* 09. PSS received IQ time-domain correlation
def PSS_received_IQ_time_domain_correlation(to_draw=True):
    
    pss_received_IQ_list = [0]*3
    pss_baseband_IQ_list = [0] * 3
    
    for N_ID_2 in N_ID_2_tuple:
        
        pss_freq_symbol_array = pss_symbol_array(N_ID_2, N_DL_RB, N_RB_sc)
        pss_baseband_IQ = ofdm_baseband_IQ_signal_generate(pss_freq_symbol_array, l, N_DL_RB, N_RB_sc, N_DL_CP, delta_f)
        pss_baseband_IQ_list[N_ID_2] = pss_baseband_IQ[-1*N:]
        pss_uu_sig = downlink_modulate(pss_baseband_IQ, t, f_0)
        pss_received_IQ_list[N_ID_2] = downlink_downconvert(pss_uu_sig, t, f_0)[-1*N:]
    
    cs_corr = array([0]*N)
    cs_list = arange(-1*(N/2), N/2, 1)
    #print cs_corr.shape, cs_list.shape
    legend_list = list()
    y_offsets = dict()
    y_offsets[(0,0)] = -10
    y_offsets[(0,1)] = 10
    y_offsets[(0,2)] = 20
    y_offsets[(1,1)] = -10
    y_offsets[(1,2)] = 10
    y_offsets[(2,2)] = -30
    for p in N_ID_2_tuple:
        for q in N_ID_2_tuple:
            if p<=q:
                for i in arange(len(cs_list)):
                    cs_corr[i] = abs(correlate(pss_baseband_IQ_list[p], roll(pss_received_IQ_list[q],cs_list[i]))[0])
                plt.plot(cs_list, cs_corr, marker='+', linestyle='-')
                legend_list.append( 'N_ID_2 %s vs. %s'%(p,q) )
                # add anotation to the max
                x, y = find_max(cs_corr)
                plt.annotate('Max of N_ID_2 %svs%s =%s @ cs=%s'%(p,q,y,cs_list[x]), xy=(cs_list[x], y), arrowprops=dict(facecolor='black', shrink=0.15), textcoords='offset points', xytext=(40, y_offsets[(p,q)]))
                
    plt.title('PSS received IQ correlation')
    plt.xlabel('Cyclic Shift')
    plt.ylabel('Correlation value')
    plt.legend(legend_list)
    plt.show()
    #plt.savefig('PSS_Uu_signal_inner_products.png', figsize=(1280,800), dpi=200, pad_inches=2)
#@+node:michael.20120305092148.1329: *4* 10. PSS Uu signal FFT N/16-point correlation
def PSS_Uu_signal_FFT_N_16_correlation():
    
    N_16 = N/16
    
    pss_uu_signal_fft_list = [0]*3
    t = arange(0, (N_CP_l+N)*T_s, T_s)
    for N_ID_2 in N_ID_2_tuple:
        pss_freq_symbol_array = pss_symbol_array(N_ID_2, N_DL_RB, N_RB_sc)
        pss_baseband_symbol = ofdm_baseband_IQ_signal_generate(pss_freq_symbol_array, l, N_DL_RB, N_RB_sc, N_DL_CP, delta_f)
        pss_uu_sig = downlink_modulate(pss_baseband_symbol, t, f_0)
        pss_uu_sig_no_cp = pss_uu_sig[-1*N:]
        pss_uu_sig_no_cp_dec = array( [0]*N_16 )
        for i in range(N_16):
            pss_uu_sig_no_cp_dec[i] = pss_uu_sig_no_cp[i*16]
        pss_uu_signal_fft_list[N_ID_2] = fft.fft( pss_uu_sig_no_cp_dec, N_16 )
    
    
    legend_list = list()
    cs_corr = array([0]*N_16)
    cs_list = arange(-1*(N_16/2), N_16/2, 1)
    y_offsets = dict()
    y_offsets[(0,0)] = 0
    y_offsets[(0,1)] = 10
    y_offsets[(0,2)] = 20
    y_offsets[(1,1)] = 0
    y_offsets[(1,2)] = 20
    y_offsets[(2,2)] = -10
    #print cs_corr.shape, cs_list.shape
    for p in N_ID_2_tuple:
        for q in N_ID_2_tuple:
            if p<=q:
                for i in arange(len(cs_list)):
                    cs_corr[i] = abs(correlate(pss_uu_signal_fft_list[p],roll(pss_uu_signal_fft_list[q],cs_list[i]))[0])
                plt.plot(cs_list, cs_corr, marker='+', linestyle='-')
                legend_list.append( 'N_ID_2 %s vs. %s'%(p,q) )
                # add anotation to the max
                x, y = find_max(cs_corr)
                plt.annotate('Max of N_ID_2 %svs%s =%s @ cs=%s'%(p,q,y,cs_list[x]), xy=(cs_list[x], y), arrowprops=dict(facecolor='black', shrink=0.15), textcoords='offset points', xytext=(60, y_offsets[(p,q)]))
    
    plt.title('PSS Uu signal decimated to 1/16 FFT correlation')
    plt.legend(legend_list)
    plt.xlabel("Cyclic Shift")
    plt.ylabel("Correlation")
    plt.show()
#@+node:michael.20120305092148.1328: *4* 09. PSS received IQ FFT N-point correlation
def PSS_received_IQ_FFT_N_correlation():
    
    pss_baseband_IQ_fft_list = [0] * 3
    pss_received_IQ_fft_list = [0]*3
    pss_received_IQ_list = [0] * 3
    for N_ID_2 in N_ID_2_tuple:
        pss_freq_symbol_array = pss_symbol_array(N_ID_2, N_DL_RB, N_RB_sc)
        pss_baseband_IQ = ofdm_baseband_IQ_signal_generate(pss_freq_symbol_array, l, N_DL_RB, N_RB_sc, N_DL_CP, delta_f)
        pss_baseband_IQ_fft_list[N_ID_2] = fft.fft(pss_baseband_IQ, N)
        pss_uu_sig = downlink_modulate(pss_baseband_IQ, t, f_0)
        pss_received_IQ_list[N_ID_2] = downlink_downconvert(pss_uu_sig, t, f_0)[-1*N:]
        #pss_received_IQ_fft_list[N_ID_2] = fft.fft( pss_received_IQ, N )
    
    
    legend_list = list()
    cs_corr = array([0]*N)
    cs_list = arange(-1*(N/2), N/2, 1)
    y_offsets = dict()
    y_offsets[(0,0)] = 0
    y_offsets[(0,1)] = 10
    y_offsets[(0,2)] = 20
    y_offsets[(1,1)] = -40
    y_offsets[(1,2)] = 20
    y_offsets[(2,2)] = -50
    #print cs_corr.shape, cs_list.shape
    for p in N_ID_2_tuple:
        for q in N_ID_2_tuple:
            if p<=q:
                for i in arange(len(cs_list)):
                    cs_corr[i] = abs(correlate(pss_baseband_IQ_fft_list[p],fft.fft(roll(pss_received_IQ_list[q],cs_list[i])))[0])
                plt.plot(cs_list, cs_corr, marker='+', linestyle='-')
                legend_list.append( 'N_ID_2 %s vs. %s'%(p,q) )
                # add anotation to the max
                x, y = find_max(cs_corr)
                plt.annotate('Max of N_ID_2 %svs%s =%s @ cs=%s'%(p,q,y,cs_list[x]), xy=(cs_list[x], y), arrowprops=dict(facecolor='black', shrink=0.15), textcoords='offset points', xytext=(60, y_offsets[(p,q)]))
    
    plt.title('PSS Uu signal FFT correlation')
    plt.legend(legend_list)
    plt.xlabel("Cyclic Shift")
    plt.ylabel("Correlation")
    plt.show()
#@+node:michael.20120309091906.1380: *4* 12. PSS Uu signal FFT after 10MHz lowpass N-point correlation
def PSS_Uu_signal_FFT_after_10MHzlowpass_N_correlation():
    
    cutoff_freq = 10*1000*1000
    nyq = 20*1000*1000
    lp_fir = firwin(80, cutoff_freq, window=('kaiser',8), nyq=20*1000*1000)
    
    pss_uu_signal_fft_list = [0]*3

    for N_ID_2 in N_ID_2_tuple:
        pss_freq_symbol_array = pss_symbol_array(N_ID_2, N_DL_RB, N_RB_sc)
        pss_baseband_symbol = ofdm_baseband_IQ_signal_generate(pss_freq_symbol_array, l, N_DL_RB, N_RB_sc, N_DL_CP, delta_f)
        pss_uu_sig = downlink_modulate(pss_baseband_symbol, t, f_0)
        pss_uu_sig_no_cp = pss_uu_sig[-1*N:]
        # apply an low pass FIR filter
        pss_uu_sig_no_cp_after_lp = convolve(pss_uu_sig_no_cp, lp_fir)
        pss_uu_signal_fft_list[N_ID_2] = fft.fft( pss_uu_sig_no_cp_after_lp, N )
    
    
    legend_list = list()
    
    cs_list = arange(-1*(N/2), N/2, 1)
    y_offsets = dict()
    y_offsets[(0,0)] = 0
    y_offsets[(0,1)] = 10
    y_offsets[(0,2)] = 20
    y_offsets[(1,1)] = -40
    y_offsets[(1,2)] = 20
    y_offsets[(2,2)] = -60
    #corr_dict = dict()
    #print cs_corr.shape, cs_list.shape
    for p in N_ID_2_tuple:
        for q in N_ID_2_tuple:
            if p<=q:
                cs_corr = array([0]*N)
                for i in arange(len(cs_list)):
                    cs_corr[i] = abs(correlate(pss_uu_signal_fft_list[p],roll(pss_uu_signal_fft_list[q],cs_list[i]))[0])
                #corr_dict[(p,q)] = cs_corr
                plt.plot(cs_list, cs_corr, marker='+', linestyle='-')
                legend_list.append( 'N_ID_2 %s vs. %s'%(p,q) )
                # add anotation to the max
                x, y = find_max(cs_corr)
                plt.annotate('Max of N_ID_2 %svs%s =%s @ cs=%s'%(p,q,y,cs_list[x]), xy=(cs_list[x], y), arrowprops=dict(facecolor='black', shrink=0.15), textcoords='offset points', xytext=(60, y_offsets[(p,q)]))
    
    plt.title('PSS Uu signal FFT after 10MHz low-pass filter correlation')
    plt.legend(legend_list)
    plt.xlabel("Cyclic Shift")
    plt.ylabel("Correlation")
    plt.show()
#@+node:michael.20120309091906.1384: *4* 15. PSS Uu signal FFT after 540KHz lowpass N-point correlation
def PSS_Uu_signal_FFT_after_540KHzlowpass_N_correlation():
    
    cutoff_freq = 540*1000
    nyq = 20*1000*1000
    lp_fir = firwin(80, cutoff_freq, window=('kaiser',8), nyq=nyq)
    
    
    pss_uu_signal_fft_list = [0]*3

    for N_ID_2 in N_ID_2_tuple:
        pss_freq_symbol_array = pss_symbol_array(N_ID_2, N_DL_RB, N_RB_sc)
        pss_baseband_symbol = ofdm_baseband_IQ_signal_generate(pss_freq_symbol_array, l, N_DL_RB, N_RB_sc, N_DL_CP, delta_f)
        pss_uu_sig = downlink_modulate(pss_baseband_symbol, t, f_0)
        pss_uu_sig_no_cp = pss_uu_sig[-1*N:]
        # apply an low pass FIR filter
        pss_uu_sig_no_cp_after_lp = convolve(pss_uu_sig_no_cp, lp_fir)
        pss_uu_signal_fft_list[N_ID_2] = fft.fft( pss_uu_sig_no_cp_after_lp, N )
    
    
    legend_list = list()
    
    cs_list = arange(-1*(N/2), N/2, 1)
    y_offsets = dict()
    y_offsets[(0,0)] = 0
    y_offsets[(0,1)] = -10
    y_offsets[(0,2)] = -20
    y_offsets[(1,1)] = -40
    y_offsets[(1,2)] = 20
    y_offsets[(2,2)] = -60
    #corr_dict = dict()
    #print cs_corr.shape, cs_list.shape
    for p in N_ID_2_tuple:
        for q in N_ID_2_tuple:
            if p<=q:
                cs_corr = array([0]*N)
                for i in arange(len(cs_list)):
                    cs_corr[i] = abs(correlate(pss_uu_signal_fft_list[p],roll(pss_uu_signal_fft_list[q],cs_list[i]))[0])
                #corr_dict[(p,q)] = cs_corr
                plt.plot(cs_list, cs_corr, marker='+', linestyle='-')
                legend_list.append( 'N_ID_2 %s vs. %s'%(p,q) )
                # add anotation to the max
                x, y = find_max(cs_corr)
                plt.annotate('Max of N_ID_2 %svs%s =%s @ cs=%s'%(p,q,y,cs_list[x]), xy=(cs_list[x], y), arrowprops=dict(facecolor='black', shrink=0.15), textcoords='offset points', xytext=(60, y_offsets[(p,q)]))
    
    plt.title('PSS Uu signal FFT after 540KHz low-pass filter correlation')
    plt.legend(legend_list)
    plt.xlabel("Cyclic Shift")
    plt.ylabel("Correlation")
    plt.show()
#@+node:michael.20120309091906.1382: *4* 14. PSS Uu signal FFT after 540KHz lowpass filter
def PSS_Uu_signal_FFT_after_540KHzlowpass_filter():
    
    func_tuple = (abs, real, imag)
    legend_tuple = ('Absolute', 'Real part', 'Image part')
    
    cutoff_freq = 540*1000
    nyq = 20*1000*1000
    lp_fir = firwin(80, cutoff_freq, window=('kaiser',8), nyq=nyq)
    
    
    for N_ID_2 in N_ID_2_tuple:
        
        pss_freq_symbol_array = pss_symbol_array(N_ID_2, N_DL_RB, N_RB_sc)
        pss_baseband_symbol = ofdm_baseband_IQ_signal_generate(pss_freq_symbol_array, l, N_DL_RB, N_RB_sc, N_DL_CP, delta_f)
        pss_uu_sig = downlink_modulate(pss_baseband_symbol, t, f_0)
        pss_uu_sig_no_cp = pss_uu_sig[-1*N:]
        
        # apply an low pass FIR filter
        pss_uu_sig_no_cp_after_lp = convolve(pss_uu_sig_no_cp, lp_fir)
        pss_uu_fft = fft.fft( pss_uu_sig_no_cp_after_lp, N )
        
        for i in (0,1,2):
            plt.subplot(131+i)
            if i==1:
                plt.title('PSS signal Uu FFT after 540KHZ low-pass filter for N_ID_2=%s'%N_ID_2)
            legend_list = list()
            plt.plot(func_tuple[i](pss_uu_fft))
            legend_list.append( legend_tuple[i] )
            plt.xlabel('N')
            #plt.ylabel('Spectrum absolute value')
            plt.legend(legend_list)
            #plt.axis( [-0.01, 0.075, -0.1, 14] )
        plt.show()
            #plt.savefig('PSS_signal_Uu_for_N_ID_2=%s.png'%N_ID_2, dpi=300)
#@+node:michael.20120309091906.1381: *4* 13. PSS Uu signal FFT after 10MHz lowpass N/16-point correlation
def PSS_Uu_signal_FFT_after_10MHzlowpass_N_16_correlation():
    
    cutoff_freq = 10*1000*1000
    nyq = 20*1000*1000
    lp_fir = firwin(80, cutoff_freq, window=('kaiser',8), nyq=20*1000*1000)
    
    N_16 = N/16
    
    pss_uu_signal_fft_list = [0]*3

    for N_ID_2 in N_ID_2_tuple:
        pss_freq_symbol_array = pss_symbol_array(N_ID_2, N_DL_RB, N_RB_sc)
        pss_baseband_symbol = ofdm_baseband_IQ_signal_generate(pss_freq_symbol_array, l, N_DL_RB, N_RB_sc, N_DL_CP, delta_f)
        pss_uu_sig = downlink_modulate(pss_baseband_symbol, t, f_0)
        pss_uu_sig_no_cp = pss_uu_sig[-1*N:]
        pss_uu_sig_no_cp_lp = convolve(pss_uu_sig_no_cp, lp_fir)
        pss_uu_sig_no_cp_dec = array( [0]*N_16 )
        for i in range(N_16):
            pss_uu_sig_no_cp_dec[i] = pss_uu_sig_no_cp_lp[i*16]
        pss_uu_signal_fft_list[N_ID_2] = fft.fft( pss_uu_sig_no_cp_dec, N_16 )
    
    
    legend_list = list()
    cs_corr = array([0]*N_16)
    cs_list = arange(-1*(N_16/2), N_16/2, 1)
    y_offsets = dict()
    y_offsets[(0,0)] = 0
    y_offsets[(0,1)] = 10
    y_offsets[(0,2)] = 20
    y_offsets[(1,1)] = 0
    y_offsets[(1,2)] = 20
    y_offsets[(2,2)] = -10
    #print cs_corr.shape, cs_list.shape
    for p in N_ID_2_tuple:
        for q in N_ID_2_tuple:
            if p<=q:
                for i in arange(len(cs_list)):
                    cs_corr[i] = abs(correlate(pss_uu_signal_fft_list[p],roll(pss_uu_signal_fft_list[q],cs_list[i]))[0])
                plt.plot(cs_list, cs_corr, marker='+', linestyle='-')
                legend_list.append( 'N_ID_2 %s vs. %s'%(p,q) )
                # add anotation to the max
                x, y = find_max(cs_corr)
                plt.annotate('Max of N_ID_2 %svs%s =%s @ cs=%s'%(p,q,y,cs_list[x]), xy=(cs_list[x], y), arrowprops=dict(facecolor='black', shrink=0.15), textcoords='offset points', xytext=(60, y_offsets[(p,q)]))
    
    plt.title('PSS Uu signal after 10MHz low-pass filter then decimated to 1/16 FFT correlation')
    plt.legend(legend_list)
    plt.xlabel("Cyclic Shift")
    plt.ylabel("Correlation")
    plt.show()
#@+node:michael.20120309091906.1386: *4* 16. PSS Uu signal FFT after 540KHz lowpass N/16-point correlation
def PSS_Uu_signal_FFT_after_540KHzlowpass_N_16_correlation():
    
    cutoff_freq = 540*1000
    nyq = 20*1000*1000
    lp_fir = firwin(80, cutoff_freq, window=('kaiser',8), nyq=20*1000*1000)
    
    N_16 = N/16
    
    pss_uu_signal_fft_list = [0]*3

    for N_ID_2 in N_ID_2_tuple:
        pss_freq_symbol_array = pss_symbol_array(N_ID_2, N_DL_RB, N_RB_sc)
        pss_baseband_symbol = ofdm_baseband_IQ_signal_generate(pss_freq_symbol_array, l, N_DL_RB, N_RB_sc, N_DL_CP, delta_f)
        pss_uu_sig = downlink_modulate(pss_baseband_symbol, t, f_0)
        pss_uu_sig_no_cp = pss_uu_sig[-1*N:]
        pss_uu_sig_no_cp_lp = convolve(pss_uu_sig_no_cp, lp_fir)
        pss_uu_sig_no_cp_dec = array( [0]*N_16 )
        for i in range(N_16):
            pss_uu_sig_no_cp_dec[i] = pss_uu_sig_no_cp_lp[i*16]
        pss_uu_signal_fft_list[N_ID_2] = fft.fft( pss_uu_sig_no_cp_dec, N_16 )
    
    
    legend_list = list()
    cs_corr = array([0]*N_16)
    cs_list = arange(-1*(N_16/2), N_16/2, 1)
    y_offsets = dict()
    y_offsets[(0,0)] = 0
    y_offsets[(0,1)] = 10
    y_offsets[(0,2)] = 20
    y_offsets[(1,1)] = 0
    y_offsets[(1,2)] = 20
    y_offsets[(2,2)] = -10
    #print cs_corr.shape, cs_list.shape
    for p in N_ID_2_tuple:
        for q in N_ID_2_tuple:
            if p<=q:
                for i in arange(len(cs_list)):
                    cs_corr[i] = abs(correlate(pss_uu_signal_fft_list[p],roll(pss_uu_signal_fft_list[q],cs_list[i]))[0])
                plt.plot(cs_list, cs_corr, marker='+', linestyle='-')
                legend_list.append( 'N_ID_2 %s vs. %s'%(p,q) )
                # add anotation to the max
                x, y = find_max(cs_corr)
                plt.annotate('Max of N_ID_2 %svs%s =%s @ cs=%s'%(p,q,y,cs_list[x]), xy=(cs_list[x], y), arrowprops=dict(facecolor='black', shrink=0.15), textcoords='offset points', xytext=(60, y_offsets[(p,q)]))
    
    plt.title('PSS Uu signal after 540KHz low-pass filter then decimated to 1/16 FFT correlation')
    plt.legend(legend_list)
    plt.xlabel("Cyclic Shift")
    plt.ylabel("Correlation")
    plt.show()
#@+node:michael.20120309091906.1379: *4* 11. PSS Uu signal FFT after 10MHz lowpass filter
def PSS_Uu_signal_FFT_after_10MHzlowpass_filter():
    
    func_tuple = (abs, real, imag)
    legend_tuple = ('Absolute', 'Real part', 'Image part')
    
    cutoff_freq = 10*1000*1000
    nyq = 20*1000*1000
    lp_fir = firwin(80, cutoff_freq, window=('kaiser',8), nyq=20*1000*1000)
    
    
    for N_ID_2 in N_ID_2_tuple:
        
        pss_freq_symbol_array = pss_symbol_array(N_ID_2, N_DL_RB, N_RB_sc)
        pss_baseband_symbol = ofdm_baseband_IQ_signal_generate(pss_freq_symbol_array, l, N_DL_RB, N_RB_sc, N_DL_CP, delta_f)
        pss_uu_sig = downlink_modulate(pss_baseband_symbol, t, f_0)
        pss_uu_sig_no_cp = pss_uu_sig[-1*N:]
        # apply an low pass FIR filter
        pss_uu_sig_no_cp_after_lp = convolve(pss_uu_sig_no_cp, lp_fir)
        pss_uu_fft = fft.fft( pss_uu_sig_no_cp_after_lp, N )
        
        for i in (0,1,2):
            plt.subplot(131+i)
            if i==1:
                plt.title('PSS signal Uu FFT after 10MHZ low-pass filter for N_ID_2=%s'%N_ID_2)
            legend_list = list()
            plt.plot(func_tuple[i](pss_uu_fft))
            legend_list.append( legend_tuple[i] )
            plt.xlabel('N')
            #plt.ylabel('Spectrum absolute value')
            plt.legend(legend_list)
            #plt.axis( [-0.01, 0.075, -0.1, 14] )
        plt.show()
            #plt.savefig('PSS_signal_Uu_for_N_ID_2=%s.png'%N_ID_2, dpi=300)
#@+node:michael.20120305092148.1325: *4* 09. PSS received IQ decimated to 1/16
def PSS_received_IQ_decimated_16():
    
    for N_ID_2 in N_ID_2_tuple:
        
        pss_freq_symbol_array = pss_symbol_array(N_ID_2, N_DL_RB, N_RB_sc)
        pss_baseband_symbol = ofdm_baseband_IQ_signal_generate(pss_freq_symbol_array, l, N_DL_RB, N_RB_sc, N_DL_CP, delta_f)
        pss_uu_sig = downlink_modulate(pss_baseband_symbol, t, f_0)
        length = len(pss_uu_sig)/16
        pss_uu_sig_dec = array([0]*length)
        for i in arange(len(pss_uu_sig_dec)):
            pss_uu_sig_dec[i] = pss_uu_sig[i*16]
        t_ = arange(0, length*16*T_s, 16*T_s)
        #print len(t), len(pss_uu_sig_dec)
        
        #plt.plot(t, symbol_array)
        plt.cla()
        legend_list = list()
        plt.plot(t_, pss_uu_sig_dec)
        legend_list.append( 'N_ID_2=%s'%(N_ID_2) )
        plt.title('PSS downconverted Uu signal decimated to 1/16 for N_ID_2=%s'%N_ID_2)
        plt.xlabel('Time (ms)')
        plt.ylabel('Signal level')
        plt.legend(legend_list)
        plt.axis( [-0.01, 0.075, -11, 7] )
        plt.show()
        #plt.savefig('PSS_signal_Uu_decimated64_for_N_ID_2=%s.png'%N_ID_2, dpi=300)
#@+node:michael.20120305092148.1324: *4* 09. PSS received IQ spectrum correlation
def PSS_received_IQ_convolve(to_draw=True):
    
    pss_iq_sig_list = [0]*3
    pss_baseband_IQ_list = [0] * 3
    
    for N_ID_2 in N_ID_2_tuple:
        
        pss_freq_symbol_array = pss_symbol_array(N_ID_2, N_DL_RB, N_RB_sc)
        pss_baseband_IQ_list[N_ID_2] = ofdm_baseband_IQ_signal_generate(pss_freq_symbol_array, l, N_DL_RB, N_RB_sc, N_DL_CP, delta_f)
        pss_uu_sig = downlink_modulate(pss_baseband_IQ_list[N_ID_2], t, f_0)
        pss_iq_sig_list[N_ID_2] = downlink_downconvert(pss_uu_sig, t, f_0)[-1*N:]
    
    sig_sample_num = N
    cs_corr = array([0]*sig_sample_num)
    cs_list = arange(-1*(sig_sample_num/2), sig_sample_num/2, 1)
    #print cs_corr.shape, cs_list.shape
    plt.cla()
    legend_list = list()
    y_offsets = dict()
    y_offsets[(0,0)] = -10
    y_offsets[(0,1)] = 10
    y_offsets[(0,2)] = 20
    y_offsets[(1,1)] = -10
    y_offsets[(1,2)] = 10
    y_offsets[(2,2)] = -30
    for p in N_ID_2_tuple:
        for q in N_ID_2_tuple:
            if p<=q:
                for i in arange(len(cs_list)):
                    cs_corr[i] = abs(sum(convolve(pss_baseband_IQ_list[p], roll(pss_iq_sig_list[q],cs_list[i]))[0]))
                plt.plot(cs_list, cs_corr, marker='+', linestyle='-')
                legend_list.append( 'N_ID_2 %s vs. %s'%(p,q) )
                # add anotation to the max
                x, y = find_max(cs_corr)
                plt.annotate('Max of N_ID_2 %svs%s =%s @ cs=%s'%(p,q,y,cs_list[x]), xy=(cs_list[x], y), arrowprops=dict(facecolor='black', shrink=0.15), textcoords='offset points', xytext=(40, y_offsets[(p,q)]))
                
    plt.title('PSS received IQ correlation')
    plt.xlabel('Cyclic Shift')
    plt.ylabel('Correlation value')
    plt.legend(legend_list)
    plt.show()
    #plt.savefig('PSS_Uu_signal_inner_products.png', figsize=(1280,800), dpi=200, pad_inches=2)
#@+node:michael.20120305092148.1331: *3* z_fft_of_ZC
from numpy import *
import matplotlib.pyplot as plt



#@+others
#@+node:michael.20120305092148.1316: *4* Zadoff-Chu seq
def ZC( n, N_ZC, q ):
    '''
    give the n-th element of 0 cyclic shift Z-C sequence with root index q and length N_ZC.
    '''
    return exp(-1j*2*pi*q*n*(n+1)/2/N_ZC)
#@-others


def z_fft_of_ZC():
    
    l = 0
    N_DL_RB = 110
    N_RB_sc = 12
    N_DL_CP = 0 # normal DL CP
    N_ZC = 61
    
    zc_seq = array( [0]*(N_ZC), dtype=complex128 )
    #print len(zc_seq)
    for i in arange(N_ZC):
        zc_seq[i] = ZC(i, N_ZC, 7)
    
    plt.subplot(141)
    plt.plot(abs(zc_seq))
    
    #plt.subplot(142)
    #plt.plot(abs(fft.fft(zc_seq)))
    
    plt.subplot(142)
    n_seq = arange(N_ZC)
    x = 0
    for k in arange(N_ZC):
        x += zc_seq[k]*exp(1j*2*pi*k*n_seq/N_ZC)
    x = x/N_ZC
    plt.plot(abs(x))
    
    #plt.subplot(143)
    #plt.plot(abs(fft.ifft(fft.fft(zc_seq))))
    
    plt.subplot(144)
    plt.plot(abs(fft.ifft(zc_seq)))
    
    plt.show()

#@-others

test_enabling_bits = 0b111111111111

# 01. PSS spectrum before OFDM generation
if test_enabling_bits & (1<<0):
    PSS_spectrum_before_OFDM_generation()

# 02. PSS correlation in freq domain before OFDM gen
if test_enabling_bits & (1<<1):
    PSS_corr_in_freq_domain_before_OFDM_gen()

# 03. PSS baseband IQ signal in time domain
if test_enabling_bits & (1<<2):
    PSS_baseband_IQ_signal_in_time_domain()

# 04. PSS baseband IQ spectrum
if test_enabling_bits & (1<<3):
    PSS_baseband_IQ_spectrum()

# 05. PSS baseband IQ spectrum correlation
if test_enabling_bits & (1<<4):
    PSS_baseband_IQ_spectrum_correlation()

# 06. PSS Uu signal 
if test_enabling_bits & (1<<5):
    PSS_signal_Uu()

# 07. PSS signal Uu downconversion
if test_enabling_bits & (1<<6):
    PSS_signal_Uu_downconversion()

# 08. PSS Uu signal downconverted correlation
if test_enabling_bits & (1<<7):
    PSS_Uu_signal_downconverted_correlation()

# 09. PSS received IQ spectrum
if test_enabling_bits & (1<<8):
    PSS_received_IQ_spectrum()

# 10. PSS Uu signal downconverted decimated to 1/16 correlation
if test_enabling_bits & (1<<9):
    PSS_Uu_signal_downconverted_decimated_1_16_correlation()

# 11. PSS detect
if test_enabling_bits & (1<<10):
    test_PSS_detect_in_baseband_IQ()
    
# 12. Channel estimation using PSS
if test_enabling_bits & (1<<11):
    test_channel_estimation_using_PSS()


#@-others
#@-leo
