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


combinemass = 1

if not combinemass:
    fig,axarr = plt.subplots(2,3,figsize=(24,15),sharex=True,sharey=True)
    axarr = np.reshape(axarr,6)
else:
    fig,axarr = plt.subplots(2,2,figsize=(16,15),sharex=True,sharey=True)
    axarr = np.reshape(axarr,4)
#Gets rid of objects with bad ellipticities
filtar = fluxdata['ar']>0


#Plot the data with error bars
#Counter
c = 0
plotdata = 'ar'
ylabel = 'Axis Ratio'
savename = 'AxisRatio'

#fluxdata['n']=np.log10(fluxdata['n'])
#fluxdata['SMD']=np.log10(divz(fluxdata['LMASS'],(4*np.pi*fluxdata['re_kpc']**2)))



ms=12
lwbw=2

notbad = np.logical_not(baddata)

for ax in axarr:
    if c in [0,1,2]:
        massfilt = fluxdata['LMASS']<9.5
    else:
        massfilt = fluxdata['LMASS']>=9.5
    if c in [0,3]:
        col = 'good'
        filt = allgood
        if combinemass: filt = notbad
        color='blue'
    elif c in [1,4]:
        col = 'low'
        filt = somelow
        if combinemass: filt = (fluxdata.OBJID < 0)
        color='orange'
    else:
        col = 'bad'
        filt = baddata
        color='red'
    #ax.errorbar(fluxdata[filt][filtar]['av'],fluxdata[filt][filtar]['ar'],xerr=fluxdata[filt][filtar]['dav1'],color=color,marker='o',ms=4,lw=0.5,ls='None')
    #ax2.errorbar(fluxdata[filt][filtar]['av'],fluxdata[filt][filtar]['re_kpc'],xerr=fluxdata[filt][filtar]['dav1'],color=color,marker='o',ms=4,lw=0.5,ls='None')
    #Titles, axes, legends
    acount = 0
    filttype = (fluxdata[plotdata]>-98.9)
    if c==0:
        ax.set_ylabel(ylabel+', LMASS < 9.5',fontsize = axisfont)
    if c==3:
        ax.set_ylabel(ylabel+', LMASS >= 9.5',fontsize = axisfont)
    ax.set_xlabel('Av (mag)',fontsize = axisfont)
    ax.tick_params(labelsize = ticksize, size=ticks)
    filters = np.logical_and(filt,massfilt)
    filters = np.logical_and(filters,filttype)
    ax.errorbar(fluxdata[filters]['av'],fluxdata[filters][plotdata],xerr=fluxdata[filters]['dav1'],color=color,marker='o',ms=4,lw=0.5,ls='None')
    if c in [0,3]:
        loc1 = np.sin(22.5/180*np.pi)
        loc2 = np.sin(45.0/180*np.pi)
        loc3 = np.sin(67.5/180*np.pi)
        mr1 = (fluxdata[filters]['ar']<loc1)
        mr2 = np.logical_and(fluxdata[filters]['ar']>=loc1,fluxdata[filters]['ar']<loc2)
        mr3 = np.logical_and(fluxdata[filters]['ar']>=loc2,fluxdata[filters]['ar']<loc3)
        mr4 = (fluxdata[filters]['ar']>=loc3)
        med1 = np.median(fluxdata[filters][mr1].av)
        med2 = np.median(fluxdata[filters][mr2].av)
        med3 = np.median(fluxdata[filters][mr3].av)
        med4 = np.median(fluxdata[filters][mr4].av)
        med751 = np.percentile(fluxdata[filters][mr1].av,75)
        med752 = np.percentile(fluxdata[filters][mr2].av,75)
        med753 = np.percentile(fluxdata[filters][mr3].av,75)
        med754 = np.percentile(fluxdata[filters][mr4].av,75)
        emed1 = np.sqrt(biweight_midvariance(fluxdata[filters][mr1].av))
        emed2 = np.sqrt(biweight_midvariance(fluxdata[filters][mr2].av))
        emed3 = np.sqrt(biweight_midvariance(fluxdata[filters][mr3].av))
        emed4 = np.sqrt(biweight_midvariance(fluxdata[filters][mr4].av))
        ax.errorbar(med1,loc1/2,xerr=emed1,color='black',marker='o',ms=ms,lw=lwbw,ls='None')
        ax.errorbar(med2,(loc1+loc2)/2,xerr=emed2,color='black',marker='o',ms=ms,lw=lwbw,ls='None')
        ax.errorbar(med3,(loc2+loc3)/2,xerr=emed3,color='black',marker='o',ms=ms,lw=lwbw,ls='None')
        ax.errorbar(med4,(1+loc3)/2,xerr=emed4,color='black',marker='o',ms=ms,lw=lwbw,ls='None')
        ax.errorbar(med751,loc1/2,xerr=emed1,color='red',marker='o',ms=ms,lw=lwbw,ls='None')
        ax.errorbar(med752,(loc1+loc2)/2,xerr=emed2,color='red',marker='o',ms=ms,lw=lwbw,ls='None')
        ax.errorbar(med753,(loc2+loc3)/2,xerr=emed3,color='red',marker='o',ms=ms,lw=lwbw,ls='None')
        ax.errorbar(med754,(1+loc3)/2,xerr=emed4,color='red',marker='o',ms=ms,lw=lwbw,ls='None')
        ax.plot((-100,100),(loc1,loc1),color='black',ls='--')
        ax.plot((-100,100),(loc2,loc2),color='black',ls='--')
        ax.plot((-100,100),(loc3,loc3),color='black',ls='--')
    ax.set_xlim(-0.5,5)
    ax.set_ylim(0,1)
    c = c+1
    if (combinemass and (c in [1,4])): c = c+1


fig.tight_layout()
if combinemass: fig.savefig(figout + 'Av_'+savename+'_combmass.pdf')
else:fig.savefig(figout + 'Av_'+savename+'_mass.pdf')
plt.close(fig)

