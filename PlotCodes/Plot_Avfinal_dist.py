#Plots the Av magnitude due to the balmer decerment
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import sys, os, string
import pandas as pd
from astropy.io import fits
import collections

#Folder to save the figures
figout = '/Users/blorenz/COSMOS/Reports/2018/Images/'

#The location with the file for all of our data
fluxdatapath = '/Users/blorenz/COSMOS/COSMOSData/lineflux.txt'

#File for the MAD of the difference in flux of duplicates in each line
maddatapath = '/Users/blorenz/COSMOS/COSMOSData/linemad.txt'

#File for the Avs of every object
avdatapath = '/Users/blorenz/COSMOS/COSMOSData/balmer_avs.txt'

#Read in the mad of the lines 
mad_df = ascii.read(maddatapath).to_pandas()

#Read in the vs 
Av_df = ascii.read(avdatapath).to_pandas()

#Read the datafile:
fluxdata = ascii.read(fluxdatapath).to_pandas()

#Fontsizes for plotting
axisfont = 18
ticksize = 16
titlefont = 24
legendfont = 16
textfont = 16

#Division function
def divz(X,Y):
        return X/np.where(Y,Y,Y+1)*np.not_equal(Y,0)

Avs = []
Averrs = []
for i in range(0,len(Av_df)):
    Avs.append(Av_df.iloc[i][Av_df.iloc[i]['useAv']])
    if (Av_df.iloc[i]['useAv']=='AvMedian'):
        Averrs.append(0)
    else: Averrs.append(Av_df.iloc[i][Av_df.iloc[i]['useAv']+'err'])
Av_df['Av'] = Avs
Av_df['eAv'] = Averrs

fluxdata = pd.merge(fluxdata,Av_df,on='fluxfile')

fig,ax2 = plt.subplots(1,figsize=(10,7))

HaHbidx = (fluxdata['useAv'] == 'AvHaHb')
HaHgidx = (fluxdata['useAv'] == 'AvHaHg')
Ha_avgidx = (fluxdata['useAv'] == 'AvHa_avg')
Medianidx = (fluxdata['useAv'] == 'AvMedian')

#colors = ['red','blue','green','orange','black']
colors = ['blue','blue','blue','blue','blue']


ax2.errorbar(fluxdata[HaHbidx]['LMASS'],fluxdata[HaHbidx]['Av'],yerr=fluxdata[HaHbidx]['eAv'],color=colors[0],marker='o',ms=4,lw=0.5,label='Ha/Hb',ls='None')
#ax2.plot(fluxdata[HbHgidx]['LMASS'],fluxdata[HbHgidx]['Av'],color=colors[1],marker='o',ms=4,lw=0.5,ls='None',label='Hb/Hg')
ax2.errorbar(fluxdata[HaHgidx]['LMASS'],fluxdata[HaHgidx]['Av'],yerr=fluxdata[HaHgidx]['eAv'],color=colors[2],marker='o',ms=4,lw=0.5,ls='None', label='Ha/Hg')
ax2.errorbar(fluxdata[Ha_avgidx]['LMASS'],fluxdata[Ha_avgidx]['Av'],yerr=fluxdata[Ha_avgidx]['eAv'],color=colors[3],marker='o',ms=4,lw=0.5,ls='None',label='Average')
ax2.errorbar(fluxdata[Medianidx]['LMASS'],fluxdata[Medianidx]['Av'],yerr=fluxdata[Medianidx]['eAv'],color=colors[4],marker='o',ms=4,lw=0.5,ls='None',label='Median Av')
#ax2.legend()

#ax2.set_title('Av vs Mass',fontsize = titlefont)
ax2.set_ylabel('Av (mag)',fontsize = axisfont)
ax2.set_xlabel('log(Stellar Mass) (Msun)',fontsize = axisfont)
'''

ax2.errorbar(np.log10(fluxdata[HaHbidx]['6563_fix_flux']),fluxdata[HaHbidx]['Av'],yerr=fluxdata[HaHbidx]['eAv'],color=colors[0],marker='o',ms=4,lw=0.5,label='Ha/Hb',ls='None')
#ax2.plot(fluxdata[HbHgidx]['LMASS'],fluxdata[HbHgidx]['Av'],color=colors[1],marker='o',ms=4,lw=0.5,ls='None',label='Hb/Hg')
ax2.errorbar(np.log10(fluxdata[HaHgidx]['6563_fix_flux']),fluxdata[HaHgidx]['Av'],yerr=fluxdata[HaHgidx]['eAv'],color=colors[2],marker='o',ms=4,lw=0.5,ls='None', label='Ha/Hg')
ax2.errorbar(np.log10(fluxdata[Ha_avgidx]['6563_fix_flux']),fluxdata[Ha_avgidx]['Av'],yerr=fluxdata[Ha_avgidx]['eAv'],color=colors[3],marker='o',ms=4,lw=0.5,ls='None',label='Average')
#ax2.errorbar(np.log10(fluxdata[Medianidx]['6563_fix_flux']),fluxdata[Medianidx]['Av'],yerr=fluxdata[Medianidx]['eAv'],color=colors[4],marker='o',ms=4,lw=0.5,ls='None',label='Median Av')
#ax2.legend()

#ax2.set_title('Av vs Mass',fontsize = titlefont)
ax2.set_ylabel('Av (mag)',fontsize = axisfont)
ax2.set_xlabel('log(H$\\alpha$ Flux) (Msun)',fontsize = axisfont)
ax2.set_xlim(0,5)
'''


ax2.tick_params(labelsize = ticksize)
    

plt.tight_layout()
fig.savefig(figout+'AvFinal_dist.pdf')
plt.close()
