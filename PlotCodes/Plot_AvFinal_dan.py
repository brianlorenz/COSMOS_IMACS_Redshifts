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
fluxdatapath = '/Users/blorenz/COSMOS/COSMOSData/lineflux_red.txt'

#File for the Avs of every object
avdatapath = '/Users/blorenz/COSMOS/COSMOSData/balmer_avs.txt'
#Read in the vs 
Av_df = ascii.read(avdatapath).to_pandas()
Av_df = Av_df.drop('LMASS',axis=1)

#Read the datafile:
fluxdata = ascii.read(fluxdatapath).to_pandas()


fluxdata = pd.merge(fluxdata,Av_df,on='fluxfile')

#The location to store the scale and its stddev of each line
qualdatapath = '/Users/blorenz/COSMOS/COSMOSData/dataqual.txt'
#Read in the scale of the lines 
dataqual = ascii.read(qualdatapath).to_pandas()
d = {'True': True, 'False': False}

#Fontsizes for plotting
axisfont = 18
ticksize = 16
titlefont = 24
legendfont = 16
textfont = 16
ticks = 8

#Division function
def divz(X,Y):
        return X/np.where(Y,Y,Y+1)*np.not_equal(Y,0)

lines=['6563_fix','4861']

#Filter the data
goodlines = [dataqual[line+'_good'].map(d) for line in lines]
#Needs to be good in all lines to be good
allgood = np.logical_and.reduce(goodlines)
#Needs to be bad in any line to be bad
badlines = [dataqual[line+'_bad'].map(d) for line in lines]
baddata = np.logical_or.reduce(badlines)
lowlines = [dataqual[line+'_low'].map(d) for line in lines]
#Needs to be low in any line to be low, and also not bad in a line
somelow = np.logical_and(np.logical_or.reduce(lowlines),np.logical_not(baddata))


#Find the duplicates
dupids = [item for item, count in collections.Counter(fluxdata[allgood]['OBJID']).items() if count > 1]
dupobjsarr = []
for obj in dupids:
    dupobjsarr.append(fluxdata[fluxdata.OBJID == obj].OBJID)

dup=1


colorsmask=['blue','blue','green','green','orange','orange','black']
masks=['b6','b7','e6','e7','i6','j7','f6']
labels = ['b6-b7','e6-e7','i6-j7']



if dup:
    fig,dupax = plt.subplots(1,1,figsize=(8.5,7),sharex=True,sharey=True)
    axarr = [1,2]
else: fig,axarr = plt.subplots(1,3,figsize=(24,7),sharex=True,sharey=True)


#colors = ['red','blue','green','orange','black']
#colors = ['blue','blue','blue','blue','blue']


#ax2.errorbar(fluxdata[HaHbidx]['LMASS'],fluxdata[allgood]['av'],yerr=fluxdata[allgood]['dav1'],color=colors[0],marker='o',ms=4,lw=0.5,ls='None')

#Plot the data with error bars
#Counter
c = 0
for ax in axarr:
    if ax==2:
        break
    if c==0:
        col = 'good'
        filt = allgood
        color='blue'
    elif c==1:
        col = 'low'
        filt = somelow
        color='yellow'
    else:
        col = 'bad'
        filt = baddata
        color='red'
    #p = fluxdata
    if (c==0 and dup):
        ax = dupax
        for j in range(0,3):
            ax.plot(-100,-100,color=colorsmask[2*j],ms=4,marker='o',label=labels[j])
        ax.legend()
    if dup: ax.errorbar(fluxdata[filt]['LMASS'],fluxdata[filt]['av'],yerr=fluxdata[filt]['dav1'],color='grey',marker='o',ms=4,lw=0.5,ls='None',zorder=3)
    else: ax.errorbar(fluxdata[filt]['LMASS'],fluxdata[filt]['av'],yerr=fluxdata[filt]['dav1'],color=color,marker='o',ms=4,lw=0.5,ls='None')
    #Titles, axes, legends
    if c==0:
        ax.set_ylabel('Av (mag)',fontsize = axisfont)
    ax.set_xlabel('log(Stellar Mass) (Msun)',fontsize = axisfont)
    ax.tick_params(labelsize = ticksize, size=ticks)
    #Draw lines between duplicats:
    if (c==0 and dup):
            for i in range(0,len(dupobjsarr)):
                pdup = fluxdata.iloc[dupobjsarr[i].index]                
                ax.plot(pdup['LMASS'],pdup['av'],color='red',marker='.',ms=0,lw=0.5,ls='-',zorder=15)
                ax.plot(pdup['LMASS'].iloc[0],pdup['av'].iloc[0],color=colorsmask[masks.index(pdup.iloc[0]['Mask'])],ls='-',ms=4,marker='o',lw=1,zorder=10,label=None)
                ax.plot(pdup['LMASS'].iloc[1],pdup['av'].iloc[1],color=colorsmask[masks.index(pdup.iloc[1]['Mask'])],ls='-',ms=4,marker='o',lw=1,zorder=11,label=None)
    c = c+1
    ax.set_xlim(8.95,10.05)
    ax.set_ylim(-0.5,4)
fig.tight_layout()
if dup: fig.savefig(figout + 'av_mass_dup.pdf')
else: fig.savefig(figout + 'av_mass.pdf')
plt.close(fig)



dup = 0


HaHbok = (fluxdata['AvHaHbok'] > 0)
allgoodavs = np.logical_and(allgood,HaHbok)

fig,ax = plt.subplots(1,1,figsize=(8.5,7),sharex=True,sharey=True)


#Plot the data with error bars
#Counter
c = 0
if 1==1:
    if c==0:
        col = 'good'
        filt = allgoodavs
        color='blue'
    if (c==0 and dup):
        for j in range(0,3):
            ax.plot(-100,-100,color=colorsmask[2*j],ms=4,marker='o',label=labels[j])
        ax.legend()
    if dup: ax.errorbar(fluxdata[filt]['AvHaHb'],fluxdata[filt]['av'],xerr=fluxdata[filt]['AvHaHberr'],yerr=fluxdata[filt]['dav1'],color='grey',marker='o',ms=4,lw=0.5,ls='None',zorder=3)
    else: ax.errorbar(fluxdata[filt]['AvHaHb'],fluxdata[filt]['av'],xerr=fluxdata[filt]['AvHaHberr'],yerr=fluxdata[filt]['dav1'],color=color,marker='o',ms=4,lw=0.5,ls='None')
    #Titles, axes, legends
    if c==0:
        ax.set_xlabel('Av H$\\alpha$/H$\\beta$ (mag)',fontsize = axisfont)
    ax.set_ylabel('Av Fit (mag)',fontsize = axisfont)
    ax.tick_params(labelsize = ticksize, size=ticks)
    #Draw lines between duplicats:
    if (c==0 and dup):
            for i in range(0,len(dupobjsarr)):
                pdup = fluxdata.iloc[dupobjsarr[i].index]                
                ax.plot(pdup['AvHaHb'],pdup['av'],color='red',marker='.',ms=0,lw=0.5,ls='-',zorder=15)
                ax.plot(pdup['AvHaHb'].iloc[0],pdup['av'].iloc[0],color=colorsmask[masks.index(pdup.iloc[0]['Mask'])],ls='-',ms=4,marker='o',lw=1,zorder=10,label=None)
                ax.plot(pdup['AvHaHb'].iloc[1],pdup['av'].iloc[1],color=colorsmask[masks.index(pdup.iloc[1]['Mask'])],ls='-',ms=4,marker='o',lw=1,zorder=11,label=None)
    c = c+1
    #Plot a unity line
    ax.plot((-10,1000),(-10,1000),color='black',ls='--')
    ax.set_xlim(-0.5,4)
    ax.set_ylim(-0.5,4)
fig.tight_layout()
if dup: fig.savefig(figout + 'avfit_avHaHb_dup.pdf')
else: fig.savefig(figout + 'avfit_avHaHb_mass.pdf')
plt.close(fig)
