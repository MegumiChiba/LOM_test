import matplotlib.pyplot as plt
import numpy as np
import os, sys
import struct
import h5py 
from eventHist import *
from HDFwriterGen2DOMPMTinspection import HDFwriterGen2DOMPMTinspection

# allow ignoring hits with charge less than specified cut, e.g. 0.1 pC
# charge_cut = 0.2
charge_cut = None

# make plots of wuBase hit waveforms, ascii output format V1
# as recorded by take-data*.py script
#
# Usage: python plot-spe-data-cut.py basename infile1 [infile2] [...]
# 
# Output: out.pdf and hist.pdf

# assume data take in new high-bandwidth configuration, which was made
# with higher overall gain too
bandwidth="high"

# make a filename root from the first argument
basename='lom_wubase'#sys.argv[1]
infile=sys.argv[1]
  
# calibration constants to convert from raw adc values into uA & pC
ch1_uA_per_count = 0.2494
#ch1_pC_per_count_sample = 0.0512#0.004157  changed to Chiba value
ch1_pC_per_count_sample = 2.0/4096/1670.0/60e6*1e12 # pC/1 ADC sample calibrated by micro base  0.004873
#conversion_ch1 = 2.0/4096/1670.0/60e6*1e12 # pC/1 ADC sample calibrated by micro base  0.004873
conversion_ch1 = 2.0*1.85/4096/1670/60e6*1e12 # measured by Tsuji for v1.4c
#conversion_ch2 = conversion_ch1/0.0128 # measured by Chiba
conversion_ch2 = conversion_ch1/0.0115 # measured by Tsuji for v1.4c
yourname = "Chiba"
description = "Data collected with no injection."

#default number
MCUID    = 1234
PMTID = 9000
pwmfreq  = 12345
temperature = 25
PMTIDstr = "BB{0}".format(PMTID)
userdict = {"Lasersetting":-1, "position":np.array([0,0,0]), "B-field":np.array([0,0,0])}
wubaseID = 0000
runtype = 1
voltage10 = 80
dacvalue = 700
Nsample_max = 100

h_ch1all=eventHist(250.,350.,100,basename,"Ch1 Value (All)","Events")
h_ch1event=eventHist(250.,350.,100,basename,"Ch1 Value (Current event)","Events")
h_ch1ped=eventHist(250.,350.,100,basename,"Ch1 Pedestal","Events")
h_ch1max=eventHist(-20.,200.,110,basename,"Ch1 Maximum (subtracted)","Events")
h_ch1area=eventHist(-100.,1000.,110,basename,"Ch1 Pulse Sum (subtracted)","Events")
h_ch1imax=eventHist(-5.,50.,110,basename,"Ch1 Peak Current (uA)","Events")
h_ch1charge=eventHist(-0.2,2.0,110,basename,"Ch1 Charge (pC)","Events")
h_dt=eventHist(0.,0.1,100,basename,"Time between hits (sec)","Events")
h_log10dt=eventHist(-8.,2.,1000,basename,"Log10(dt/sec) between hits","Events")

# read mode fsm definitions
vlist = [76,78,80,82,84,86,88,90,92,94,96]

for i in vlist:
  print(i)
  lasttime=None
  filename    = infile + str(i) + '.dat'
  outfilename = infile + str(i)+ '.hd5'
  Timestamp = []
  Nsample = []
  hitData=[]  # will hold (ch1,ch2,discraw,discsync) for each saved hit
  maxDisplayWaveforms = 100
  waveformCount = 0
  
  try:
    open(filename,'rb')
      
  except:
    blank = i
    print(filename,'is not found')
    continue

  if os.path.getsize(filename)==0:
    print(filename,'is 0 size')
    continue

  with open(filename,'rb') as f:
    while True:
      try:
        hdr=f.read(19)
      except:
        print('failed to read data')
        break
      if not hdr: break
      try:
        n_samples=struct.unpack('<H',hdr[1:3])[0]
      except:
        print('failed to read data')
        break
      if n_samples>255: 
        print(n_samples>255)
        break
      if n_samples<=0: 
        print(n_samples<=0)
        break
      time=struct.unpack('<Q',hdr[5:11]+bytes([0,0]))[0]
      #print('time OK')
      try:
        ch1_data=f.read(2*n_samples)
      except:
        print('failed to read data')
        break
      if(len(ch1_data)<n_samples): 
        print((len(ch1_data)<n_samples))
        break
      try:
        ch1=struct.unpack('<%dH'%n_samples,ch1_data)
        ch2_data=f.read(2*n_samples)
      except:
        print('failed')
        break
      #ch2_data=f.read(2*n_samples)
      if(len(ch2_data)<n_samples): 
        print(len(ch2_data)<n_samples)
        break
      ch2=struct.unpack('<%dH'%n_samples,ch2_data)

      #print(len(Nsample),len(Timestamp))
      Ch1 = np.asarray(ch1).astype(np.uint16)
      if len(Ch1)>Nsample_max:
        print('Lengs',len(ch1),'is too long waveform')
        break
      while len(Ch1)<=Nsample_max:
        Ch1=np.append(Ch1,0)
        #print(len(Ch1),'/',Nsample_max)
      # Ch1 = np.asarray(ch1).astype(np.uint16)
      if waveformCount == 0:
          adcs_ch1 = Ch1
      elif waveformCount > 0:
          adcs_ch1 = np.vstack((adcs_ch1,Ch1))

      Ch2 = np.asarray(ch2).astype(np.uint16)
      while len(Ch2)!=Nsample_max:
        Ch2=np.append(Ch2,0)
      if waveformCount == 0:
          adcs_ch2 = Ch2
      elif waveformCount > 0:
          adcs_ch2 = np.vstack((adcs_ch2,Ch2))

      ymax=max(ch1)
      if ymax==0: raise ValueError   # empty waveform = not good
      #h_ch1event.clear()
      for ch1_value in ch1: 
        h_ch1event.increment(ch1_value)
      imax,nmax=h_ch1event.getMaximum()
      # compute sloppily interpolated max in (imax-.5,imax+.5), call that pedestal
      nleft=h_ch1event.getBinValue(imax-1)
      nright=h_ch1event.getBinValue(imax+1)
      if nmax+nleft+nright<=0: continue
      pedestal=h_ch1event.xmin+h_ch1event.dx*(imax+(nright-nleft)/(nmax+nleft+nright))
      # add up samples in window (-3,+6) around max of waveform
      tmax=ch1.index(ymax)
      tleft=tmax-3
      if tleft<0: tleft=0
      tright=tmax+6+1
      if tright>=len(ch1): tright=len(ch1)-1
      area=sum(ch1[tleft:tright])-(tright-tleft)*pedestal
      charge=area*ch1_pC_per_count_sample
      if charge_cut is not None:
        if charge<charge_cut: continue
      h_ch1area.increment(area)
      h_ch1charge.increment(charge)
      h_ch1ped.increment(pedestal)
      # subtract pedestal from peak, add to histogram
      h_ch1max.increment(ymax-pedestal)
      h_ch1imax.increment((ymax-pedestal)*ch1_uA_per_count)
      for ch1_value in ch1: 
        h_ch1all.increment(ch1_value)
      #chiba
      nsamples = len(ch1)
      Timestamp.append(time)
      Nsample.append(nsamples)
      waveformCount += 1
      if len(hitData)<maxDisplayWaveforms:
        hitData.append((ch1,ch2))
      # compute time separation between this hit and the previous one (assuming time stamps at 60Msps)
      if lasttime:
        if time<=lasttime:
          h_log10dt.increment(-100.) # use underflow bin for unphysical values
          h_dt.increment(-100.)
        else:
          h_log10dt.increment(math.log10((time-lasttime)/60e6))
          h_dt.increment((time-lasttime)/60e6)
      lasttime=time

    print(f"Total reported hits = {waveformCount}")
    print(f"Number of underflows in log(dt) histogram = {h_dt.underflow}")


    # a class to save HDF file
    hdf = HDFwriterGen2DOMPMTinspection(Nsample_max)

    # set the metadata
    hdf.fill_metadata(PMTID, PMTIDstr, wubaseID, MCUID, yourname, runtype, voltage10, \
      pwmfreq, dacvalue, temperature, conversion_ch1, conversion_ch2, description, userdict)

    # save the waveform and their basic values
    Nsample = np.array(Nsample).astype(np.uint16)
    try:
      adcs_ch1 = np.array(adcs_ch1).astype(np.uint16)
      
    except:
      blank = i
      print(filename,'does not have data')
      continue

    adcs_ch2 = np.array(adcs_ch2).astype(np.uint16)
    hdf.fill(Nsample, Timestamp, Timestamp, adcs_ch1, adcs_ch2)

    hdf.fit_v0()

    hdf.write("{0}".format(outfilename))

'''
# plot sample waveforms    
plt.clf()
fig,(ax1,ax2,ax3,ax4) = plt.subplots(4)
times=[i*1000./60. for i in range(10000)]
nmax=0
ch1_min=5000.
for (ch1,ch2) in hitData:
  n=len(ch1)-6
  if n>nmax: nmax=n
  ax1.plot(times[:n],ch1[6:],'-',linewidth=0.5)
  ax2.plot(times[:n],ch2[6:],'-',linewidth=0.5)
  ax3.plot(times[:n],ch1[6:],'-',linewidth=0.5)
  ax4.plot(times[:n],ch2[6:],'-',linewidth=0.5)
  ch1_min=min(ch1_min,min(ch1))
ax1.set(ylabel='Channel 1')
ax2.set(ylabel='Channel 2')
ax3.set(ylabel='Channel 1')
ax4.set(ylabel='Channel 2')
ax4.set(xlabel='Time (nsec)')
ax1.set_xlim(0,600)
ax1.set_ylim(250,500)
ax2.set_xlim(0,600)
ax2.set_ylim(3750,3850)
ax3.set_xlim(0,600)
ax3.set_ylim(0,4096)
ax3.set_yticks((0,1000,2000,3000,4000))
ax4.set_xlim(0,600)
ax4.set_ylim(0,4096)
ax4.set_yticks((0,1000,2000,3000,4000))
ax1.label_outer()
ax2.label_outer()
ax3.label_outer()
ax4.label_outer()
plt.savefig(basename+'_waveforms.pdf')
plt.show()

# plot histograms of raw adc values
plt.clf()
f=plt.figure(figsize=(10.,10.5))
f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=.3, hspace=.4)
hist_list=[h_ch1all,h_ch1ped,h_ch1max,h_ch1area]
for i in range(len(hist_list)):
    ax=f.add_subplot(3,2,i+1)
    hist_list[i].plot(ax,'b')
    hist_list[i].autoSetLimits(ax,scale='lin')
    ax.set_xlabel(hist_list[i].xlabel)
    ax.set_ylabel(hist_list[i].ylabel)
    ax.set_title(hist_list[i].title)
plt.savefig(basename+'_adc_hist.pdf')
plt.show()

# plot histograms in physical units
plt.clf()
f=plt.figure(figsize=(10.,10.5))
f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=.3, hspace=.4)
hist_list=[h_ch1imax,h_ch1charge]
for i in range(len(hist_list)):
    ax=f.add_subplot(3,2,i+1)
    hist_list[i].plot(ax,'b')
    hist_list[i].autoSetLimits(ax,scale='lin')
    ax.set_xlabel(hist_list[i].xlabel)
    ax.set_ylabel(hist_list[i].ylabel)
    ax.set_title(hist_list[i].title)
plt.savefig(basename+'_pulse_hist.pdf')
plt.show()

# plot more histograms
plt.clf()
f=plt.figure(figsize=(10.,10.5))
f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=.3, hspace=.4)
hist_list=[h_log10dt,h_log10dt,h_dt,h_dt]
scale_list=['lin','log','lin','log']  # which plots use lin or log scaling on y axis
for i in range(len(hist_list)):
    ax=f.add_subplot(3,2,i+1)
    hist_list[i].plot(ax,'b')
    hist_list[i].autoSetLimits(ax,scale=scale_list[i])
    ax.set_xlabel(hist_list[i].xlabel)
    ax.set_ylabel(hist_list[i].ylabel)
    ax.set_title(hist_list[i].title)
plt.savefig(basename+'_dt_hist.pdf')
plt.show()
'''
