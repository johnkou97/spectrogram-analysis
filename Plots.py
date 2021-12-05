import h5py
import os
import numpy as np
import math
import scipy
import pywt
from scipy import signal
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
from scipy.interpolate import CubicSpline as spline
from scipy.fftpack import fft, fftshift ,ifft,rfft,fftfreq,rfftfreq
c=2.9979e10
G=6.67408e-8
Msun=1.989e33
Length = G*Msun/c**2
Time = Length/c
Frequency=1/Time

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

BAM=['002','003','004','009','010','022','035','036','046','048','053','057','058','059','061','065','070','080','089','090','091','092','093','098','107','121','122','123','124','126','128']

#functions for q=1
def f20(M,R6):
    return 8.943+4.059*M-1.332*R6-.358*(M**2)-.182*R6*M+.048*(R6**2)

def fspir(M,R8):
    return 6.264+1.929*M-.645*R8+.881*(M**2)-.311*R8*M+.03*(R8**2)

def fpeak(M,R6):
    return 13.822-0.576*M-1.375*R6+.479*(M**2)-.073*R6*M+.044*(R6**2)

#functions for all cases
def f20_a(M,R6):
    return 9.586+4.09*M-1.427*R6+.048*(M**2)-.261*R6*M+.055*(R6**2)

def fspir_a(M,R8):
    return 5.846+1.75*M-.555*R8+1.002*(M**2)-.316*R8*M+.026*(R8**2)

def fpeak_a(M,R8):
    return 10.942-.369*M-.987*R8+1.095*(M**2)-.201*R8*M+.036*(R8**2)

#functions for fft
def fre_do(x,y,mass):
    fd=fft(y)
    N=len(y)
    if (N % 2) == 1:
        N=N+1
    T=x[1]-x[0]
    xf = np.linspace(0.0, 1.0/(2.0*T), int(N/2))/mass
    fq=fftfreq(len(y))
    mask=fq>=0
    fd=2.0*(fd/N)
    fd=fd[mask]
    fd=abs(fd)
    return xf,fd

def analyze(rh,mass):



    rhM=rh[:,1]
    time=rh[:,0]

    peaks,prop=scipy.signal.find_peaks(abs(rhM))
    ampls=rhM[peaks]
    merg=np.amax(abs(ampls))
    merg=np.where(abs(ampls)==merg)
    merg=int(merg[0])
    t0=peaks[merg]

    ampl=rhM[t0:]
    tim=time[t0:]

    #ampl=rhM
    #tim=time

    tuk=signal.tukey(len(ampl),0.03)
    dat=ampl*tuk

    fq,fd=fre_do(tim,dat,mass)

    mx=np.where(fd==np.amax(fd))[0][0]
    freq=fq[mx]
    amp=fd[mx]
    #fig=plt.figure()
    #plt.plot((fq*Frequency),fd)
    #plt.xlim(0,5500)
    return fq,fd,tim,dat#,fig

if os.path.exists('results'):
    pass
else:
    os.mkdir('results')

#metadata
SLy=1
H4=2
MS1=3
MPA1=4
ALF2=5
MS1b=6
ms1b=6
ENG=7

q2 = np.zeros(31)
mas2=np.zeros(31)
eos=np.zeros(31)


i=0
for m in range(0,2):
    for k in range(0,10):
        for j in range(0,10):
            name = 'metadata/BAM:0%s%s%s.txt' %(m,k,j)
            try:
                f=open(name)
                lines=f.readlines()
                exec(lines[10])
                q2[i] = id_mass_ratio
                exec(lines[8])
                mas2[i]=id_mass
                if i==0:
                    eos[i]=8 #by hand because 2H can not be given a value
                if i>0:
                    exec(lines[15])
                    eos[i]=id_eos

                i=i+1
            except OSError:
                pass

#isolate q=1 cases
index=[]
j=0
for i in range(len(q2)):
    if q2[i]==1:
        index.insert(j,i)
        j=j+1

#calculate M_chirp
Mc=pow(q2/pow(1+q2,2),3/5)*mas2

#calculate R_1.6 and R_1.8
m_r1=np.load('tid_def/SLy.npy')
m_r2=np.load('tid_def/H4.npy')
m_r3=np.load('tid_def/MS1.npy')
m_r4=np.load('tid_def/MPA1.npy')
m_r5=np.load('tid_def/ALF2.npy')
m_r6=np.load('tid_def/MS1b.npy')
m_r7=np.load('tid_def/ENG.npy')
m_r8=np.load('tid_def/2H.npy')

mx=np.amax(m_r1[0])
idx=np.where(m_r1[0]==mx)
idx=idx[0][0]
#cs1=spline(m_r1[0][1:idx],k_l1[0][1:idx])
cs11=spline(m_r1[0][1:idx],m_r1[1][1:idx])

mx=np.amax(m_r2[0])
idx=np.where(m_r2[0]==mx)
idx=idx[0][0]
#cs2=spline(m_r2[0][1:idx],k_l2[0][1:idx])
cs21=spline(m_r2[0][1:idx],m_r2[1][1:idx])

mx=np.amax(m_r3[0])
idx=np.where(m_r3[0]==mx)
idx=idx[0][0]
#cs3=spline(m_r3[0][1:idx],k_l3[0][1:idx])
cs31=spline(m_r3[0][1:idx],m_r3[1][1:idx])

mx=np.amax(m_r4[0])
idx=np.where(m_r4[0]==mx)
idx=idx[0][0]
cs41=spline(m_r4[0][1:idx],m_r4[1][1:idx])

mx=np.amax(m_r5[0])
idx=np.where(m_r5[0]==mx)
idx=idx[0][0]
cs51=spline(m_r5[0][1:idx],m_r5[1][1:idx])

mx=np.amax(m_r6[0])
idx=np.where(m_r6[0]==mx)
idx=idx[0][0]
cs61=spline(m_r6[0][1:idx],m_r6[1][1:idx])

mx=np.amax(m_r7[0])
idx=np.where(m_r7[0]==mx)
idx=idx[0][0]
cs71=spline(m_r7[0][1:idx],m_r7[1][1:idx])

mx=np.amax(m_r8[0])
idx=np.where(m_r8[0]==mx)
idx=idx[0][0]
cs81=spline(m_r8[0][1:idx],m_r8[1][1:idx])

r68=np.zeros((len(BAM),2))

for i in range(len(BAM)):
    if eos[i]==1:
        r68[i,0]=cs11(1.6)*Length/1.0e5
        r68[i,1]=cs11(1.8)*Length/1.0e5

    elif eos[i]==2:
        r68[i,0]=cs21(1.6)*Length/1.0e5
        r68[i,1]=cs21(1.8)*Length/1.0e5
    elif eos[i]==3:
        r68[i,0]=cs31(1.6)*Length/1.0e5
        r68[i,1]=cs31(1.8)*Length/1.0e5

    elif eos[i]==4:
        r68[i,0]=cs41(1.6)*Length/1.0e5
        r68[i,1]=cs41(1.8)*Length/1.0e5

    elif eos[i]==5:
        r68[i,0]=cs51(1.6)*Length/1.0e5
        r68[i,1]=cs51(1.8)*Length/1.0e5

    elif eos[i]==6:
        r68[i,0]=cs61(1.6)*Length/1.0e5
        r68[i,1]=cs61(1.8)*Length/1.0e5

    elif eos[i]==7:
        r68[i,0]=cs71(1.6)*Length/1.0e5
        r68[i,1]=cs71(1.8)*Length/1.0e5

    elif eos[i]==8:
        r68[i,0]=cs81(1.6)*Length/1.0e5
        r68[i,1]=cs81(1.8)*Length/1.0e5


#typ=np.load('pdf/typ.npy')
#q=1 analysis
f_2=np.zeros(len(BAM))
f_s=np.zeros(len(BAM))
f_p=np.zeros(len(BAM))
f_0=np.zeros(len(BAM))
for i in index:
    f_2[i]=f20(Mc[i],r68[i,0])
    f_p[i]=fpeak(Mc[i],r68[i,0])
    f_s[i]=fspir(Mc[i],r68[i,1])
f_0=2*f_p-f_2


#for i in [0]:
#    name = 'data/BAM:0'+BAM[i]+'.h5'
#    file=h5py.File(name,'r')
#    dat = list(file["/rh_22"])
#    rh = np.array(file["/rh_22/%s" %dat[-1]])
#    freq,amp,time,dat,fig1=analyze(rh,mas2[i])
#    timems=np.zeros(len(time))
#    for j in range(len(timems)):
#        timems[j]=time[j]*Time*1000
#    ax=plt.subplot()
#    ax.axvline(x=(f_p[i]*Mc[i])*1000,color='r',label='peak')
#    ax.axvspan((f_p[i]*Mc[i])*1000-196, (f_p[i]*Mc[i])*1000+196, alpha=0.3, color='grey')
#    ax.axvline(x=(f_2[i]*Mc[i])*1000,color='g',label='2-0')
#    ax.axvspan((f_2[i]*Mc[i])*1000-229, (f_2[i]*Mc[i])*1000+229, alpha=0.3, color='yellow')
#    ax.axvline((f_s[i]*Mc[i])*1000,color='orange',label='spiral')
#    ax.axvspan((f_s[i]*Mc[i])*1000-286, (f_s[i]*Mc[i])*1000+286, alpha=0.3, color='cyan')
#    ax.axvline((f_0[i]*Mc[i])*1000,linestyle="--",color='grey',label='2+0')
#    plt.xlabel('frequency (Hz)')
#    ax.xaxis.set_major_locator(MultipleLocator(500))
#    ax.xaxis.set_minor_locator(MultipleLocator(100))
#    plt.xlabel('Time(ms)')
#    plt.ylabel('Frequency(Hz)')
#    plt.legend()
#    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
#    if typ[i]==1:
#        plt.text(0.01,0.9,'Type I',transform=ax.transAxes,fontsize=12, verticalalignment='top', bbox=props)
#    elif typ[i]==2:
#        plt.text(0.01,0.9,'Type II',transform=ax.transAxes,fontsize=12, verticalalignment='top', bbox=props)
#    else:
#        plt.text(0.01,0.9,'Type III',transform=ax.transAxes,fontsize=12, verticalalignment='top', bbox=props)
#    #plt.savefig('results/q1/linear/BAM:0'+BAM[i]+'.jpg')
#    plt.savefig('paper/q1/linear/BAM:0'+BAM[i]+'.pdf',format='pdf')

#    plt.close()

#    fc = 2.5
#    dt=(time[1]-time[0])*Time*1000*mas2[i]
#    band = 2.5
#    wavelet = 'cmor'+str(band)+'-'+str(fc)
#    widths = fc/np.linspace(fc-2.0, fc+2.0, 400)/dt
#    cwtmatr, freqs = pywt.cwt(dat, widths, wavelet, dt)
#    power = abs(cwtmatr)
#
#    fig, ax = plt.subplots(figsize=(10, 6))
#    ax.pcolormesh((timems-timems[0])*mas2[i], freqs, power,norm=colors.LogNorm(0.01,np.amax(power)),cmap='jet')
#    plt.xlabel('Time(ms)')
#    plt.ylabel('Frequency(Hz)')
#    ax.axhspan((f_p[i]*Mc[i])-.196, (f_p[i]*Mc[i])+.196, alpha=0.25, color='black')
#    plt.text(0.2,f_p[i]*Mc[i],r'$f_{\rm{peak}}$')
#    ax.axhspan((f_s[i]*Mc[i])-.286, (f_s[i]*Mc[i])+.286, alpha=0.25, color='black')
#    plt.text(0.2,f_s[i]*Mc[i],r'$f_{\rm{spiral}}$')
#    ax.axhspan((f_2[i]*Mc[i])-.229, (f_2[i]*Mc[i])+.229, alpha=0.25, color='black')
#    plt.text(0.2,f_2[i]*Mc[i],r'$f_{2-0}$')
#    plt.axhline(y=f_0[i]*Mc[i],linestyle=':',color='black')
#    plt.text(0.2,f_0[i]*Mc[i]+0.1,r'$f_{2+0}$')
#    #plt.plot([0,timems[-1]-timems[0]],[f_p[i]*Mc[i],f_p[i]*Mc[i]],'--',color='grey',lw=2)
#    #plt.text(0.2,f_p[i]*Mc[i]-.1,r'$f_{\rm{peak}}$')
#    #plt.plot([0,timems[-1]-timems[0]],[f_0[i]*Mc[i],f_0[i]*Mc[i]],':',color='grey',lw=2)
#    #plt.text(0.2,f_0[i]*Mc[i]-.1,r'$f_{2+0}$')
#    #plt.plot([0,timems[-1]-timems[0]],[f_s[i]*Mc[i],f_s[i]*Mc[i]],'--',color='grey',lw=2)
#    #plt.text(0.2,f_s[i]*Mc[i]-.1,r'$f_{\rm{spiral}}$')
#    #if f_2[i]*Mc[i]>1.2:
#    #    plt.plot([0,timems[-1]-timems[0]],[f_2[i]*Mc[i],f_2[i]*Mc[i]],'--',color='grey',lw=2)
#    #    plt.text(0.2,f_2[i]*Mc[i]-.1,r'$f_{2-0}$')
#    ax.xaxis.set_major_locator(MultipleLocator(1))
#    ax.xaxis.set_minor_locator(MultipleLocator(.2))
#    ax.yaxis.set_major_locator(MultipleLocator(.5))
#    ax.yaxis.set_minor_locator(MultipleLocator(.1))
#    plt.savefig('paper/spec1/BAM:0'+BAM[i]+'.png',format='png',dpi=1200)
#    #plt.savefig('paper/spec1/BAM:0'+BAM[i]+'.svg',format='svg')
#    plt.close()




#analysis for all cases
f_2_a=np.zeros(len(BAM))
f_s_a=np.zeros(len(BAM))
f_p_a=np.zeros(len(BAM))
f_0_a=np.zeros(len(BAM))
for i in range(len(BAM)):
    f_2_a[i]=f20_a(Mc[i],r68[i,0])
    f_p_a[i]=fpeak_a(Mc[i],r68[i,1])
    f_s_a[i]=fspir_a(Mc[i],r68[i,1])
f_0_a=2*f_p_a-f_2_a




for i in range(len(BAM)):
    name = 'data/BAM:0'+BAM[i]+'.h5'
    file=h5py.File(name,'r')
    dat = list(file["/rh_22"])
    rh = np.array(file["/rh_22/%s" %dat[-1]])
    freq,amp,time,dat=analyze(rh,mas2[i])
    timems=np.zeros(len(time))
    for j in range(len(timems)):
        timems[j]=time[j]*Time*1000
    #ax=plt.subplot()
    #ax.axvline(x=(f_p_a[i]*Mc[i])*1000,color='r',label='peak')
    #ax.axvspan((f_p_a[i]*Mc[i])*1000-196, (f_p_a[i]*Mc[i])*1000+196, alpha=0.3, color='grey')
    #ax.axvline(x=(f_2_a[i]*Mc[i])*1000,color='g',label='2-0')
    #ax.axvspan((f_2_a[i]*Mc[i])*1000-229, (f_2_a[i]*Mc[i])*1000+229, alpha=0.3, color='yellow')
    #ax.axvline((f_s_a[i]*Mc[i])*1000,color='orange',label='spiral')
    #ax.axvspan((f_s_a[i]*Mc[i])*1000-286, (f_s_a[i]*Mc[i])*1000+286, alpha=0.3, color='cyan')
    #ax.axvline((f_0_a[i]*Mc[i])*1000,linestyle="--",color='grey',label='2+0')
    #plt.xlabel('frequency (Hz)')
    #ax.xaxis.set_major_locator(MultipleLocator(500))
    #ax.xaxis.set_minor_locator(MultipleLocator(100))
    #plt.legend()
    #props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    #if typ[i]==1:
    #    plt.text(0.01,0.9,'Type I',transform=ax.transAxes,fontsize=12, verticalalignment='top', bbox=props)
    #elif typ[i]==2:
    #    plt.text(0.01,0.9,'Type II',transform=ax.transAxes,fontsize=12, verticalalignment='top', bbox=props)
    #else:
    #    plt.text(0.01,0.9,'Type III',transform=ax.transAxes,fontsize=12, verticalalignment='top', bbox=props)
    #plt.savefig('results/BAM:0'+BAM[i]+'.pdf',format='pdf')
    #plt.close()

    #fig=plt.figure()
    #plt.subplot(212)
    #plt.plot((freq*Frequency),amp)
    #plt.xlim(0,5000)
    #plt.xlabel('Frequency (Hz)')
    #plt.legend(['Postmerger only'])
    #plt.subplot(222)
    #plt.plot(time*Time*1000-time[0]*Time*1000,dat)
    #plt.title('Postmerger')
    #plt.subplot(221)
    #plt.plot(rh[:,0]*Time*1000,rh[:,1])
    #plt.title('Time Domain')
    #plt.savefig('results/3fig/linear/BAM:0'+BAM[i]+'.jpg')
    #plt.close()


    fc = 3
    dt=(time[1]-time[0])*Time*1000*mas2[i]
    band = 2.5
    wavelet = 'cmor'+str(band)+'-'+str(fc)
    widths = fc/np.linspace(fc-2.0, fc+2.0, 400)/dt
    cwtmatr, freqs = pywt.cwt(dat, widths, wavelet, dt)
    power = abs(cwtmatr)
    st=find_nearest(freqs,(f_p_a[i]*Mc[i]))
    maxfreq=np.zeros(power[1].size)
    position=st
    for j in range(len(maxfreq)):

        maxpos= np.where(power[:,j] == np.amax(power[position-30:position+30,j]))[0][0]
        maxfreq[j]=freqs[maxpos]
        position=maxpos

    #fig, ax = plt.subplots(figsize=(10, 6), dpi=150)
    #ax.pcolormesh(time*Time*1000-time[0]*Time*1000, freqs, power,cmap='jet')
    #plt.xlabel('Time(ms)')
    #plt.ylabel('Frequency(Hz)')
    #ax.xaxis.set_major_locator(MultipleLocator(1))
    #ax.xaxis.set_minor_locator(MultipleLocator(.2))
    #ax.yaxis.set_major_locator(MultipleLocator(.5))
    #ax.yaxis.set_minor_locator(MultipleLocator(.1))
    #plt.savefig('results/spec/linear/BAM:0'+BAM[i]+'.jpg')
    #plt.close()

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.pcolormesh(time*Time*mas2[i]*1000-time[0]*Time*mas2[i]*1000, freqs, power,norm=colors.LogNorm(0.005,np.amax(power)),cmap='jet')
    plt.xlabel('Time(ms)')
    plt.ylabel('Frequency(Hz)')
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.xaxis.set_minor_locator(MultipleLocator(.2))
    ax.yaxis.set_major_locator(MultipleLocator(.5))
    ax.yaxis.set_minor_locator(MultipleLocator(.1))
    ax.axhspan((f_p_a[i]*Mc[i])-.196, (f_p_a[i]*Mc[i])+.196, alpha=0.25, color='black')
    plt.text(0.2,f_p_a[i]*Mc[i],r'$f_{\rm{peak}}$')
    ax.axhspan((f_s_a[i]*Mc[i])-.286, (f_s_a[i]*Mc[i])+.286, alpha=0.25, color='black')
    plt.text(0.2,f_s_a[i]*Mc[i],r'$f_{\rm{spiral}}$')
    ax.axhspan((f_2_a[i]*Mc[i])-.229, (f_2_a[i]*Mc[i])+.229, alpha=0.25, color='black')
    plt.text(0.2,f_2_a[i]*Mc[i],r'$f_{2-0}$')
    plt.axhline(y=f_0_a[i]*Mc[i],linestyle=':',color='black')
    plt.text(0.2,f_0_a[i]*Mc[i]+0.1,r'$f_{2+0}$')
    plt.plot((timems-timems[0])*mas2[i],maxfreq,color='black')
    #plt.plot([0,timems[-1]-timems[0]],[f_p_a[i]*Mc[i],f_p_a[i]*Mc[i]],'--',color='grey',lw=2)
    #plt.text(0.2,f_p_a[i]*Mc[i]-.1,r'$f_{\rm{peak}}$')
    #plt.plot([0,timems[-1]-timems[0]],[f_0_a[i]*Mc[i],f_0_a[i]*Mc[i]],':',color='grey',lw=2)
    #plt.text(0.2,f_0_a[i]*Mc[i]-.1,r'$f_{2+0}$')
    #plt.plot([0,timems[-1]-timems[0]],[f_s_a[i]*Mc[i],f_s_a[i]*Mc[i]],'--',color='grey',lw=2)
    #plt.text(0.2,f_s_a[i]*Mc[i]-.1,r'$f_{\rm{spiral}}$')
    #if f_2_a[i]*Mc[i]>1.2:
    #    plt.plot([0,timems[-1]-timems[0]],[f_2_a[i]*Mc[i],f_2_a[i]*Mc[i]],'--',color='grey',lw=2)
    #    plt.text(0.2,f_2_a[i]*Mc[i]-.1,r'$f_{2-0}$')

    plt.savefig('results/BAM:0'+BAM[i]+'.png',format='png',dpi=1200)
    #plt.savefig('paper/spec/BAM:0'+BAM[i]+'.svg',format='svg')
    #plt.savefig('results/spec/log/BAM:0'+BAM[i]+'.jpg')
    plt.close()
