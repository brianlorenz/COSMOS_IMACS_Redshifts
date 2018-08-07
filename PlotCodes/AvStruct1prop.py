#Creates a BPT diagram for all objects, and a second figure that shows objects for which single lines are low
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import sys, os, string
import pandas as pd
from astropy.io import fits
import collections
from astropy.cosmology import WMAP9 as cosmo
from astropy.stats import biweight_midvariance

#Folder to save the figures
figout = '/Users/blorenz/COSMOS/Reports/2018/Images/'

#The location with the file for all of our data
fluxdatapath = '/Users/blorenz/COSMOS/COSMOSData/lineflux_red.txt'

#Location of the equivalent width data
ewdata = '/Users/blorenz/COSMOS/COSMOSData/lineew.txt'
#Read in the ew of the lines 
ew_df = ascii.read(ewdata).to_pandas()

#The location to store the scale and its stddev of each line
qualdatapath = '/Users/blorenz/COSMOS/COSMOSData/dataqual.txt'
#Read in the scale of the lines 
dataqual = ascii.read(qualdatapath).to_pandas()
d = {'True': True, 'False': False}

#File with the error array
errdatapath = '/Users/blorenz/COSMOS/COSMOSData/errs.txt'
#Read in the scale of the lines 
err_df = ascii.read(errdatapath,data_start=1,header_start=0,format='csv').to_pandas()

#Read the datafile:
fluxdata = ascii.read(fluxdatapath).to_pandas()

#File with the structural properties
spropdatapath = '/Users/blorenz/COSMOS/COSMOSData/struct_prop.txt'
#Read in the scale of the lines 
sprop_df = ascii.read(spropdatapath).to_pandas()
sprop_df = sprop_df.rename(columns={'id':'OBJID'})
fluxdata = pd.merge(fluxdata,sprop_df)

#Fontsizes for plotting
axisfont = 24
ticksize = 18
ticks = 8
titlefont = 24
legendfont = 16
textfont = 16

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



fig,axarr = plt.subplots(2,3,figsize=(24,7),sharex=True,sharey=True)
axarr = np.reshape(axarr,6)

#Gets rid of objects with bad ellipticities
filtar = fluxdata['ar']>0


#Plot the data with error bars
#Counter
c = 0
plotdata = 'ar'
ylabels = 'Axis Ratio'
savename = 'AxisRatio'

#fluxdata['n']=np.log10(fluxdata['n'])
#fluxdata['SMD']=np.log10(divz(fluxdata['LMASS'],(4*np.pi*fluxdata['re_kpc']**2)))

mr1 = (fluxdata[allgood]['ar']<0.25)
mr2 = np.logical_and(fluxdata[allgood]['ar']>=0.25,fluxdata[allgood]['ar']<0.5)
mr3 = np.logical_and(fluxdata[allgood]['ar']>=0.5,fluxdata[allgood]['ar']<0.75)
mr4 = (fluxdata[allgood]['ar']>=0.75)
med1 = np.median(fluxdata[allgood][mr1].av)
med2 = np.median(fluxdata[allgood][mr2].av)
med3 = np.median(fluxdata[allgood][mr3].av)
med4 = np.median(fluxdata[allgood][mr4].av)
med751 = np.percentile(fluxdata[allgood][mr1].av,75)
med752 = np.percentile(fluxdata[allgood][mr2].av,75)
med753 = np.percentile(fluxdata[allgood][mr3].av,75)
med754 = np.percentile(fluxdata[allgood][mr4].av,75)
emed1 = np.sqrt(biweight_midvariance(fluxdata[allgood][mr1].av))
emed2 = np.sqrt(biweight_midvariance(fluxdata[allgood][mr2].av))
emed3 = np.sqrt(biweight_midvariance(fluxdata[allgood][mr3].av))
emed4 = np.sqrt(biweight_midvariance(fluxdata[allgood][mr4].av))
ms=12
lwbw=2
        
for ax in axarr:
    if c in [0,1,2]:
        massfilt = fluxdata['LMASS']<9.5
    else:
        massfilt = fluxdata['LMASS']>=9.5
    if c in [0,3]:
        col = 'good'
        filt = allgood
        color='blue'
    elif c in [1,4]:
        col = 'low'
        filt = somelow
        color='orange'
    else:
        col = 'bad'
        filt = baddata
        color='red'
    #ax.errorbar(fluxdata[filt][filtar]['av'],fluxdata[filt][filtar]['ar'],xerr=fluxdata[filt][filtar]['dav1'],color=color,marker='o',ms=4,lw=0.5,ls='None')
    #ax2.errorbar(fluxdata[filt][filtar]['av'],fluxdata[filt][filtar]['re_kpc'],xerr=fluxdata[filt][filtar]['dav1'],color=color,marker='o',ms=4,lw=0.5,ls='None')
    #Titles, axes, legends
    acount = 0
    filttype = (fluxdata[plotdata[acount]]>-98.9)
    if c==0:
        ax.set_ylabel(ylabels[acount],fontsize = axisfont)
    ax.set_xlabel('Av (mag)',fontsize = axisfont)
    ax.tick_params(labelsize = ticksize, size=ticks)
    ax.errorbar(fluxdata[filt][massfilt][filttype]['av'],fluxdata[filt][massfilt][filttype][plotdata[acount]],xerr=fluxdata[filt][massfilt][filttype]['dav1'],color=color,marker='o',ms=4,lw=0.5,ls='None')
    if c==0:
        ax3.errorbar(med1,0.125,xerr=emed1,color='black',marker='o',ms=ms,lw=lwbw,ls='None')
        ax3.errorbar(med2,0.375,xerr=emed2,color='black',marker='o',ms=ms,lw=lwbw,ls='None')
        ax3.errorbar(med3,0.625,xerr=emed3,color='black',marker='o',ms=ms,lw=lwbw,ls='None')
        ax3.errorbar(med4,0.875,xerr=emed4,color='black',marker='o',ms=ms,lw=lwbw,ls='None')
        ax3.errorbar(med751,0.125,xerr=emed1,color='red',marker='o',ms=ms,lw=lwbw,ls='None')
        ax3.errorbar(med752,0.375,xerr=emed2,color='red',marker='o',ms=ms,lw=lwbw,ls='None')
        ax3.errorbar(med753,0.625,xerr=emed3,color='red',marker='o',ms=ms,lw=lwbw,ls='None')
        ax3.errorbar(med754,0.875,xerr=emed4,color='red',marker='o',ms=ms,lw=lwbw,ls='None')
    ax.set_xlim(-0.5,5)
    ax.set_ylim(-2.2,1)
    c = c+1


fig.tight_layout()
fig.savefig(figout + 'Av_'+savename+'_mass.pdf')
plt.close(fig)

