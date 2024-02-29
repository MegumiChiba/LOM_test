
import h5py 
import numpy as np
import datetime
import time
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit  
import warnings
warnings.filterwarnings("ignore")
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages     
import sys 

file    = sys.argv[1]
pdfname = sys.argv[2]
ch = sys.argv[3]

pdfname = f'{pdfname}_dark.pdf'
#pdf = PdfPages(f"data/lom_test/BB{PMTID}/{pdfname}")
pdf = PdfPages(f"{pdfname}")
#pdf = PdfPages(f"data/BB{PMTID}/{pdfname}")

timelist = np.array([])
sum_time = np.array([])
ratelist = np.array([])
countlist = np.array([])
blank = np.array([])

def Gauss1(x,Aspe,Vspe,Mspe):
	# Mspe=mean
	gaussian = Aspe *np.exp(-(x-Mspe)**2/(2*Vspe**2))
	return gaussian
def GainFit(Vct,Again):
    k = 0.75
    n = 10
    gainfit = Again *Vct**(k*n)
    return gainfit

for loop in range(10):
    error=0
    filename    = file + str(loop) + '.hd5'
    try:
        f = h5py.File(filename, mode="r")
    except:
        blank = np.append(blank,loop)
        print(filename,'is not found')
        continue
    print('Open',filename)

    nsample          = f["data"]["nsample"][()]
    FPGAtime         = f["data"]["FPGAtime"][()]
    FPGAtcword       = f["data"]["FPGAtcword"][()]
    charge_ch1       = f["data"]["charge_ch1"][()]        
    peak_ch1         = f["data"]["peak_ch1"][()]          
    time_ch1         = f["data"]["time_ch1"][()]          
    charge_fit_ch1   = f["data"]["charge_fit_ch1"][()]    
    peak_fit_ch1     = f["data"]["peak_fit_ch1"][()]      
    time_fit_ch1     = f["data"]["time_fit_ch1"][()]      
    pedestal_ch1     = f["data"]["pedestal_ch1"][()]     
    ADC_ch1          = f["data"]["ADC_ch1"][()]               
    charge_ch2       = f["data"]["charge_ch2"][()]        
    peak_ch2         = f["data"]["peak_ch2"][()]          
    time_ch2         = f["data"]["time_ch2"][()]          
    charge_fit_ch2   = f["data"]["charge_fit_ch2"][()]    
    peak_fit_ch2     = f["data"]["peak_fit_ch2"][()]      
    time_fit_ch2     = f["data"]["time_fit_ch2"][()]      
    pedestal_ch2     = f["data"]["pedestal_ch2"][()]     
    ADC_ch2          = f["data"]["ADC_ch2"][()]               
    Nwfm             = f["metadata"]["Nwfm"][()]
    Nsample_max      = f["metadata"]["Nsample_max"][()]          
   # voltage          = f["metadata"]["voltage"][()]       
    DACvalue         = f["metadata"]["DACvalue"][()]       
    conversion_ch1   = f["metadata"]["conversion_ch1"][()]
    conversion_ch2   = f["metadata"]["conversion_ch2"][()]
    date             = f["metadata"]["date"][()]          
    #PMTID            = f["metadata"]["PMTID"][()]   
    PMTIDstr         = f["metadata"]["PMTIDstr"][()]
    wubaseID         = f["metadata"]["wubaseID"][()]      
    runtype          = f["metadata"]["runtype"][()]       
    #OSCfreq          = f["metadata"]["OSCfreq"][()]       
    creatorname      = f["metadata"]["creatorname"][()]#.decode("utf-8")       
    description      = f["metadata"]["description"][()]#.decode("utf-8")   
    fitversion       = f["metadata"]["fitversion"][()]   
    temperature      = f["metadata"]["temperature"][()]

    print(Nwfm,'data')
    if Nwfm<500:
        continue
#    timelist = np.append(timelist,loop)

    bins_set = 100
    range_set = (0,200)#v1.4b
    range_set = (0,120)#v1.4c

    hist, bins = np.histogram(peak_ch1, bins = bins_set,range = range_set)
    x_list = []
    for i in range(len(hist)):
        x_list.append((bins[i] + bins[i+1]) / 2)

    x_data = np.array(x_list)
    y_data = hist
    limit = int(bins_set*0.15)
    limit2 = int(bins_set*0.85)
    mean = np.mean(peak_ch1[limit:])
    mean = 90 #v1.4b
    mean = 30 #v1.4c
    try:
        popt, pcov = curve_fit(Gauss1,x_data[limit:limit2],y_data[limit:limit2],p0=[1000,20,mean])
        ADCspe = popt[2]
        #print(popt)
        print('SPE  ',popt[2],'counts')
        #print('0.2PE',popt[2]*0.2,'counts')
        print('0.25PE',popt[2]*0.25,'counts')
    except RuntimeError:
        popt=[1,1,1]
        pcov= [[1, 1, 1],[1, 1, 1],[1, 1, 1]]
        print( 'Optimal parameters not found')
 #   '''
    fig, ax = plt.subplots()
    plt.hist(peak_ch1,bins = bins_set,log = True,range = range_set,histtype="step")
    plt.plot(x_data[limit:limit2],Gauss1(x_data[limit:limit2], popt[0],popt[1],popt[2]),'-',label='1PE: {:.3g} counts'.format(popt[2]))   #result line
    plt.vlines(popt[2]*0.2,0.8,100000,label='0.2PE: {:.3g} counts'.format(popt[2]*0.2))
    plt.vlines(popt[2]*0.25,0.8,100000,label='0.25PE: {:.3g} counts'.format(popt[2]*0.25))
    plt.ylim(0.8,np.max(hist)*2)
    plt.title("Peak hight distribution")
    plt.legend()
    pdf.savefig(fig)
    plt.show()
#'''
    nwf = Nwfm-1
    #nwf = 50000
    threshold = 0.25
    deltat = np.array([])
    peak_dark = np.array([])
    #peak_dark=np.append(peak_dark,1)
    timestamp = FPGAtime[0]   # *6e7[s]
    print(timestamp)
    print(nwf,'calcuration')

    for i in range(nwf-10):
        if ADC_ch1[i][35]>=int(pedestal_ch1[i])-2:
            sigindex = np.where(ADC_ch1[i][:]>=ADCspe*threshold+pedestal_ch1[i]) 
            sigindex = np.asarray(sigindex)
            
            for j in range(len(sigindex[0,:])):
                #print(i)
                if int(ADC_ch1[i][sigindex[0,j]])-int(ADC_ch1[i][sigindex[0,j]-1])>0:
                    if int(ADC_ch1[i][sigindex[0,j]])-int(ADC_ch1[i][sigindex[0,j]+1])>0:
                        try:
                            delta = ((FPGAtime[i]+sigindex[0,j])-timestamp)/6*1e-7 #[s]
                        except:
                            error=1
                            print('Error of size')
                        if delta<0.1:
                            try:
                                deltat = np.append(deltat,(((FPGAtime[i]+sigindex[0,j])-timestamp)/6*1e-7))
                            except:
                                error=1
                                print('Error of size')
                                break
                            peak_dark = np.append(peak_dark,ADC_ch1[i][sigindex[0,j]]-pedestal_ch1[i])
                        timestamp = FPGAtime[i]+sigindex[0,j]
        elif peak_ch1[i] >= ADCspe*threshold:
            #print(np.shape(FPGAtime),np.shape(time_fit_ch1),np.shape(timestamp))
            try:
                delta = ((FPGAtime[i]+time_fit_ch1[i])-timestamp)/6*1e-7 #[s]
            except:
                error=1
                print('Error of data size', len(FPGAtime),'i=',i)
                break
            if delta<0.1:
                try:
                    deltat = np.append(deltat,(((FPGAtime[i]+time_fit_ch1[i])-timestamp)/6*1e-7))
                except:
                    error=1
                    print('Error of data size',len(FPGAtime),'i=',i)
                    break
                peak_dark = np.append(peak_dark,peak_ch1[i])
            #print((((FPGAtime[i]+time_fit_ch1[i])-timestamp)/6*1e-7))
            try:
                timestamp = FPGAtime[i]+time_fit_ch1[i] # *6e7[s]
            except:
                error=1
                print('Error of data size')
                break
    if error == 0:
        deltat=deltat[np.where(deltat>0)]
        print(len(deltat[np.where(deltat<1)]),'wfs in','{:.2f}'.format(sum(deltat[np.where(deltat<1)])),'s')
        #print((sum(deltat[np.where(deltat<1)])/len(deltat[np.where(deltat<1)]))**(-1),'Hz')
        #darkcount = (sum(deltat[np.where(deltat<1)])/len(deltat[np.where(deltat<1)]))**(-1)
        darkcount = len(deltat)/np.sum(deltat)
        sum_dt = np.sum(deltat)
        print('Darkrate:{:.4g}Hz'.format(darkcount))
    
        #timelist = np.append(timelist,int(description))
        ratelist = np.append(ratelist,1/popt[0])
        if darkcount>0:
            countlist = np.append(countlist,darkcount)
            timelist = np.append(timelist,loop)
            sum_time = np.append(sum_time,sum_dt)

meanrate = np.mean(countlist)
print('mean of darkrate: {:.3g} Hz'.format(np.mean(countlist)))
print('sigma of countrate: {:.3g} Hz'.format(np.std(countlist)))
# print('mean of darkrate: {:.5g} Hz'.format(np.mean(countlist[np.where(countlist>=meanrate-10)])))
# print('sigma of countrate: {:.5g} Hz'.format(np.std(countlist[np.where(countlist>=meanrate-10)])))
#'''
fig, ax = plt.subplots()
#plt.annotate('mean of countrate: {:.5g} Hz'.format(np.mean(countlist[np.where(countlist>=meanrate-10)])),xy=(0.5,2e3), fontsize=14)
plt.scatter(timelist,countlist,label='sigma: {:.3g} Hz'.format(np.std(countlist)))
plt.ylim(np.mean(countlist[np.where(abs(countlist-meanrate)<=10)])-20,np.mean(countlist[np.where(abs(countlist-meanrate)<=10)])+20)
plt.title(f'ch{ch}  '+'Darkrate in {:.2g}s: {:.3g} Hz'.format(np.sum(sum_time),np.mean(countlist)))
plt.title(f'ch{ch}  '+'Darkrate in {:.2g}s: {:.3g} Hz'.format(np.sum(sum_time),np.mean(countlist[np.where(abs(countlist-meanrate)<=10)])))
plt.xlabel('Time')
plt.ylabel('Dark Rate[Hz]')
plt.legend()
fig.tight_layout()
pdf.savefig(fig)
plt.show()
#'''
pdf.close()

