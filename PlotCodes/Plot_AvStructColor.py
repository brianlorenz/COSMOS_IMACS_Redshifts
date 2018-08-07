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

#colorssplit = (fluxdata['BoT']<np.median(fluxdata['BoT']))
colordata = fluxdata['BoT']
colorname = 'Bulge/Total'
colorbins = 3
#colordata = np.log10(fluxdata['6563_fix_flux'])

for ax in axarr:
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
    filttype = (fluxdata[plotdata]>-98.9)
    if c==0:
        ax.set_ylabel(ylabel,fontsize = axisfont)
    ax.set_xlabel('Av (mag)',fontsize = axisfont)
    ax.tick_params(labelsize = ticksize, size=ticks)
    filters = np.logical_and(filt,filttype)
    ax.errorbar(fluxdata[filters]['av'],fluxdata[filters][plotdata],xerr=fluxdata[filters]['dav1'],color='grey',marker='o',ms=4,lw=0.5,ls='None',zorder=1,label='None')
    if colorbins:
        bins = [np.round(np.percentile(colordata,(100/colorbins)*i),3) for i in np.arange(1,colorbins)]
        print bins
        colors = ['blue','green','orange','purple']
        for i in range(0,len(bins)+1):
            if i==len(bins):
                filt2 =  np.logical_and(filters,colordata>bins[i-1])
                label = colorname+' > ' + str(bins[i-1])
            else:
                filt2 = np.logical_and(filters,colordata<bins[i])
                label = colorname+' < ' + str(bins[i])
                if i>0:
                    filt2 = np.logical_and(filt2,colordata>bins[i-1])
                    label = str(bins[i-1])+' < '+colorname+' < ' + str(bins[i])
                
            ax.scatter(fluxdata[filt2]['av'],fluxdata[filt2][plotdata],c=colors[i],marker='o',s=30,zorder=3,label=label)
            
        
        
    #colorbar = ax.scatter(fluxdata[filters]['av'],fluxdata[filters][plotdata],c=colordata[filters],marker='o',s=30,zorder=3)
    if c in [0]:
        ax.legend(loc=4,fontsize=12)
        mr1 = (fluxdata[allgood]['ar']<0.25)
        mr2 = np.logical_and(fluxdata[filters]['ar']>=0.25,fluxdata[filters]['ar']<0.5)
        mr3 = np.logical_and(fluxdata[filters]['ar']>=0.5,fluxdata[filters]['ar']<0.75)
        mr4 = (fluxdata[filters]['ar']>=0.75)
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
        ax.errorbar(med1,0.125,xerr=emed1,color='black',marker='o',ms=ms,lw=lwbw,ls='None')
        ax.errorbar(med2,0.375,xerr=emed2,color='black',marker='o',ms=ms,lw=lwbw,ls='None')
        ax.errorbar(med3,0.625,xerr=emed3,color='black',marker='o',ms=ms,lw=lwbw,ls='None')
        ax.errorbar(med4,0.875,xerr=emed4,color='black',marker='o',ms=ms,lw=lwbw,ls='None')
        ax.errorbar(med751,0.125,xerr=emed1,color='red',marker='o',ms=ms,lw=lwbw,ls='None')
        ax.errorbar(med752,0.375,xerr=emed2,color='red',marker='o',ms=ms,lw=lwbw,ls='None')
        ax.errorbar(med753,0.625,xerr=emed3,color='red',marker='o',ms=ms,lw=lwbw,ls='None')
        ax.errorbar(med754,0.875,xerr=emed4,color='red',marker='o',ms=ms,lw=lwbw,ls='None')
    ax.set_xlim(-0.5,5)
    ax.set_ylim(0,1)
    c = c+1

#cb = fig.colorbar(colorbar)
#cb.ax.tick_params(labelsize = ticksize)
    

fig.tight_layout()
fig.savefig(figout + 'Av_'+savename+'_color.pdf')
plt.close(fig)

