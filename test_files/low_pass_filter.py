#@+leo-ver=5-thin
#@+node:michael.20120309091906.1378: * @thin ./test_files/low_pass_filter.py
#@+others
#@+node:michael.20120309091906.1377: ** souce code
from numpy import *
import matplotlib.pyplot as plt
from scipy.signal import *


def low_pass_filter():
    
    T_s = 0.01/307200
    t = arange(0, 2048*T_s, T_s)
    delta_f = 15*1000
    
    f_0_start = 10*1000*1000 + delta_f * 12 * 50
    f_0_end = f_0_start + delta_f * 12
    
    f_1_start = 10*1000*1000 - delta_f * 0
    f_1_end = f_1_start - delta_f * 1
    
    sig_01 = 0
    for f in arange(f_0_start, f_0_end):
        sig_01 += cos(2*pi*f*t) - sin(2*pi*f*t)
    
    sig_02 = 0
    for f in arange(f_1_start, f_1_end, -1):
        sig_02 += cos(2*pi*f*t) - sin(2*pi*f*t)
    #sig_01 = exp(-1j*2*pi*f_0*t)
    #sig_02 = exp(-1j*2*pi*f_1*t)
    #plt.plot(t, sig_01)
    #plt.plot(t, sig_02)
    fft_01 = fft.fft(sig_01, 4096)
    fft_02 = fft.fft(sig_02, 4096)
    
    legend_list = list()
    plt.plot(abs(fft_01))
    legend_list.append('FFT of sig_01 %6.6sMHz to %6.6sMHz'%(f_0_start/1000./1000, f_0_end/1000./1000))
    plt.plot(abs(fft_02))
    legend_list.append('FFT of sig_02 %6.6sMHz to %6.6sMHz'%(f_1_start/1000./1000, f_1_end/1000./1000))
    plt.legend(legend_list)
    
    
    plt.axis( [-10, 2050, -10, 10 + max(abs(fft_01))] )
    #plt.axis( [-10, 2050, -10, 10 + max(max(abs(fft_01)), max(abs(fft_02)))] )
    plt.show()

def fir_low_pass(nyq, cutoff):
    
    b = firwin(80, cutoff, window=('kaiser', 8), nyq=nyq)
    h, w = freqz(b, worN=80)
    fig = plt.figure()
    plt.title('Digital filter frequency response')
    ax1 = fig.add_subplot(111)
    plt.semilogy(h, np.abs(w), 'b')
    plt.ylabel('Amplitude (dB)', color='b')
    plt.xlabel('Frequency (rad/sample)')
    plt.grid()
    plt.legend()
    ax2 = ax1.twinx()
    angles = unwrap(angle(w))
    plt.plot(h, angles, 'g')
    plt.ylabel('Angle (radians)', color='g')
    plt.show()
    

fir_low_pass(20*1000*1000, 10*1000*1000)
#@-others
#@-leo
