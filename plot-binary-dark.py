import matplotlib.pyplot as plt
import numpy as np
import os, sys
import struct

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
basename=sys.argv[1]
ch = sys.argv[2]
  
# calibration constants to convert from raw adc values into uA & pC
ch1_uA_per_count = 0.2494
ch1_pC_per_count_sample = 0.004157
ch1_pC_per_count_sample = 2.0/4096/1670.0/60e6*1e12 # pC/1 ADC sample calibrated by micro base  0.004873
conversion_ch1 = 2.0/4096/1670.0/60e6*1e12 # pC/1 ADC sample calibrated by micro base  0.004873
conversion_ch2 = conversion_ch1/0.0128 # measured by Chiba
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
Nsample_max = 100#comment

h_ch1all=eventHist(250.,350.,100,basename,"Ch1 Value (All)","Events")
h_ch1event=eventHist(250.,350.,100,basename,"Ch1 Value (Current event)","Events")
h_ch1ped=eventHist(250.,350.,100,basename,"Ch1 Pedestal","Events")
h_ch1max=eventHist(-20.,200.,110,basename,"Ch1 Maximum (subtracted)","Events")
h_ch1area=eventHist(-100.,1000.,110,basename,"Ch1 Pulse Sum (subtracted)","Events")
h_ch1imax=eventHist(-5.,50.,110,basename,"Ch1 Peak Current (uA)","Events")
h_ch1charge=eventHist(-0.2,2.0,110,basename,"Ch1 Charge (pC)","Events")
h_dt=eventHist(0.,0.1,100,basename,"Time between hits (sec)","Events")
h_log10dt=eventHist(-8.,2.,1000,basename,"Log10(dt/sec) between hits","Events")

infile=sys.argv[3]

for i in range(1):
  #print(i)
  lasttime=None
  filename    = infile + '.dat'#str(i) + '.dat'
  outfilename = infile + '.hd5'#str(i)+ '.hd5'
  Timestamp = []
  Nsample = []
  hitData=[]  # will hold (ch1,ch2,discraw,discsync) for each saved hit
  DT = []
  timelist = []
  timeori = []
  maxDisplayWaveforms = 100
  waveformCount = 0

  try:
    open(filename,'rb')
    print('Reading ',filename)
  except:
    blank = i
    print(filename,'is not found')
    continue

  if os.path.getsize(filename)==0:
    print(filename,'is 0 size')
    continue

  with open(filename,'rb') as f:
    while True:
      hdr=f.read(19)
      if not hdr:
          print("not hdr")
          break

      n_samples=struct.unpack('<H',hdr[1:3])[0]
      #print(n_samples) #comment out
      if n_samples>255:
          print(n_samples,'n_samples>255 continue')
          continue
          #break
      if n_samples==0: 
        print('n_samples==0 continue')
        continue
        #break
      time=struct.unpack('<Q',hdr[5:11]+bytes([0,0]))[0]
      ch1_data=f.read(2*n_samples)
      #print('data',time,ch1_data) #comment out
      if(len(ch1_data)<n_samples): break
      ch1=struct.unpack('<%dH'%n_samples,ch1_data)
      ch2_data=f.read(2*n_samples)
      if(len(ch2_data)<n_samples): break
      ch2=struct.unpack('<%dH'%n_samples,ch2_data)

      Ch1 = np.asarray(ch1).astype(np.uint16)
      if len(Ch1)>Nsample_max:
        print('Lengs',len(ch1),'is too long waveform continue')
        continue
      while len(Ch1)<=Nsample_max:
        Ch1=np.append(Ch1,0)
        #print(Ch1) #comment out
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
      h_ch1event.clear()
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
          DT.append((time-lasttime)/60e6)
          timelist.append(time/60e6)
          timeori.append(time)
      lasttime=time
    # a class to save HDF file
    hdf = HDFwriterGen2DOMPMTinspection(Nsample_max)

    # set the metadata
    hdf.fill_metadata(PMTID, PMTIDstr, wubaseID, MCUID, yourname, runtype, voltage10, \
      pwmfreq, dacvalue, temperature, conversion_ch1, conversion_ch2, description, userdict)

    # save the waveform and their basic values
    Nsample = np.array(Nsample).astype(np.uint16)
    #print(len(adcs_ch1),type(adcs_ch1))
    try:
      adcs_ch1 = np.array(adcs_ch1).astype(np.uint16)
      
    except:
      blank = i
      print(filename,'does not have data')
      continue

    adcs_ch2 = np.array(adcs_ch2).astype(np.uint16)
    hdf.fill(Nsample, timeori, Timestamp, adcs_ch1, adcs_ch2)

    hdf.fit_v0()

    hdf.write("{0}".format(outfilename))

    print(f"Total reported hits = {waveformCount}")
    print(f"Number of underflows in log(dt) histogram = {h_dt.underflow}")

    counts=0
    sumt = 0
    deadtime=1
    for dt in DT:
        if dt<=0.1:
            counts+=1
            sumt+=dt
        else:
            deadtime += 1

    darkrate=counts/sumt
    print('Darkrate in {:.2g} s is {:.4g} Hz'.format(sumt,darkrate))

    plt.clf()
    f=plt.figure(figsize=(10.,10.5))
    f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=.3, hspace=.4)
    hist_list=[h_log10dt,h_log10dt,h_dt,h_dt]
    hist_name=['h_log10dt','h_log10dt','h_dt','h_dt']
    scale_list=['lin','log','lin','log']  # which plots use lin or log scaling on y axis
    for i in range(len(hist_list)):
        ax=f.add_subplot(3,2,i+1)
        hist_list[i].plot(ax,'b')
        hist_list[i].autoSetLimits(ax,scale=scale_list[i])
        ax.set_xlabel(hist_list[i].xlabel)
        ax.set_ylabel(hist_list[i].ylabel)
        ax.set_title(hist_name[i])
    ax.set_title('ch{} Darkrate in {:.2g} s: {:.3g} Hz'.format(ch,sumt,darkrate))
    #plt.show()
plt.savefig(basename+'_dt_hist.pdf')



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
plt.show()
plt.savefig(basename+'_waveforms.pdf')

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
plt.show()
plt.savefig(basename+'_adc_hist.pdf')

# plot histograms in physical units
plt.clf()
f=plt.figure(figsize=(10.,10.5))
f.subplots_adjust(left=None, bottom=Non
='lin')
    ax.set_xlabel(hist_list[i].xlabel)
    ax.set_ylabel(hist_list[i].ylabel)
    ax.set_title(hist_list[i].title)
plt.show()
plt.savefig(basename+'_pulse_hist.pdf')
#''
# plot more histograms
    plt.clf()
    f=plt.figure(figsize=(10.,10.5))
    f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=.3, hspace=.4)
    hist_list=[h_log10dt,h_log10dt,h_dt,h_dt]
    hist_name=['h_log10dt','h_log10dt','h_dt','h_dt']
    scale_list=['lin','log','lin','log']  # which plots use lin or log scaling on y axis
    for i in range(len(hist_list)):
        ax=f.add_subplot(3,2,i+1)
        hist_list[i].plot(ax,'b')
        hist_list[i].autoSetLimits(ax,scale=scale_list[i])
        ax.set_xlabel(hist_list[i].xlabel)
        ax.set_ylabel(hist_list[i].ylabel)
        ax.set_title(hist_name[i])
    ax.set_title('ch{} Darkrate in {:.2g} s: {:.3g} Hz'.format(ch,sumt,darkrate))
    #plt.show()
plt.savefig(basename+'_dt_hist.pdf')

plt.plot(DT)
plt.xlabel('count')
plt.ylabel('dt')
plt.show()
'''
