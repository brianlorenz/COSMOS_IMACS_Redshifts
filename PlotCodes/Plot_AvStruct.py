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



fig,axarr = plt.subplots(1,3,figsize=(24,7),sharex=True,sharey=True)
fig2,axarr2 = plt.subplots(1,3,figsize=(24,7),sharex=True,sharey=True)
fig3,axarr3 = plt.subplots(1,3,figsize=(24,7),sharex=True,sharey=True)
fig4,axarr4 = plt.subplots(1,3,figsize=(24,7),sharex=True,sharey=True)
fig5,axarr5 = plt.subplots(1,3,figsize=(24,7),sharex=True,sharey=True)

#Gets rid of objects with bad ellipticities
filtar = fluxdata['ar']>0


#Plot the data with error bars
#Counter
c = 0
plotdata = ['n','re_kpc','ar','BoT','SMD']
ylabels = ['Sersic Index','Re (kpc)','Axis Ratio','Bulge to Total Ratio','log(Surface Mass Density (Msun/$kpc^2$))']
savename = ['Sersic','Re','AxisRatio','BulgeTotal','SMD']

#fluxdata['n']=np.log10(fluxdata['n'])
#fluxdata['BoT']=np.log10(fluxdata['BoT'])
fluxdata['SMD']=np.log10(divz(fluxdata['LMASS'],(4*np.pi*fluxdata['re_kpc']**2)))
uvjcolors = np.logical_and(fluxdata['UmV']>1.3,fluxdata['VmJ']<1.5)
uvjquies = np.logical_and(uvjcolors,fluxdata['UmV']>(fluxdata['VmJ']*0.88+0.69))
uvjsf = np.logical_not(uvjquies)

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
    ax2 = axarr2[c]
    ax3 = axarr3[c]
    ax4 = axarr4[c]
    ax5 = axarr5[c]
    if c==0:
        col = 'good'
        filt = allgood
        color='blue'
    elif c==1:
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
    for a in [ax,ax2,ax3,ax4,ax5]:
        filttype = (fluxdata[plotdata[acount]]>-98.9)
        if c==0:
            a.set_ylabel(ylabels[acount],fontsize = axisfont)
        a.set_xlabel('Av (mag)',fontsize = axisfont)
        a.tick_params(labelsize = ticksize, size=ticks)
        if acount == 4:
            filt5 = np.logical_and(filt,filttype)
            quifilt = np.logical_and(filt5,uvjquies)
            sffilt = np.logical_and(filt5,uvjsf)
            a.errorbar(fluxdata[quifilt]['av'],fluxdata[quifilt][plotdata[acount]],xerr=fluxdata[quifilt]['dav1'],color='red',marker='o',ms=4,lw=0.5,ls='None')
            a.errorbar(fluxdata[sffilt]['av'],fluxdata[sffilt][plotdata[acount]],xerr=fluxdata[sffilt]['dav1'],color='blue',marker='o',ms=4,lw=0.5,ls='None')
        else:
            a.errorbar(fluxdata[np.logical_and(filt,filttype)]['av'],fluxdata[np.logical_and(filt,filttype)][plotdata[acount]],xerr=fluxdata[np.logical_and(filt,filttype)]['dav1'],color=color,marker='o',ms=4,lw=0.5,ls='None')
        acount = acount + 1
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
    ax.set_ylim(0,5)
    ax2.set_xlim(-0.5,5)
    ax3.set_xlim(-0.5,5)
    ax4.set_xlim(-0.5,5)
    ax5.set_xlim(-0.5,5)
    ax5.set_ylim(-2.2,1)
    #Draw lines between duplicats:
    c = c+1
    #ax.set_xlim(8.95,10.05)
    #ax.set_ylim(-0.5,4)

fig.tight_layout()
fig2.tight_layout()
fig3.tight_layout()
fig4.tight_layout()
fig5.tight_layout()
fig.savefig(figout + 'Av_'+savename[0]+'.pdf')
fig2.savefig(figout + 'Av_'+savename[1]+'.pdf')
fig3.savefig(figout + 'Av_'+savename[2]+'.pdf')
fig4.savefig(figout + 'Av_'+savename[3]+'.pdf')
fig5.savefig(figout + 'Av_'+savename[4]+'.pdf')
plt.close(fig)
plt.close(fig2)
plt.close(fig3)
plt.close(fig4)
plt.close(fig5)
