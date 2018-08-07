#Creates a UVJ diagram split by good, low, and bad measurements for specified lines
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

lines=['6583_fix','6563_fix','4861','5007']

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

fig,axarr = plt.subplots(1,3,figsize=(24,7),sharex=True,sharey=True)
#Gets rid of objects with bad ellipticities
filtar = fluxdata['ar']>0


#Plot the data with error bars
#Counter
c = 0
plotdata = 'ar'
ylabel = 'Axis Ratio'
savename = 'AxisRatio'


ms=12
lwbw=2

notbad = np.logical_not(baddata)

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
    ax.plot(fluxdata[filt]['VmJ'],fluxdata[filt]['UmV'],color=color,marker='o',ms=4,lw=0.5,ls='None')
    ax.plot((-100,0.69),(1.3,1.3),color='black')
    ax.plot((1.5,1.5),(2.01,100),color='black')
    xline = np.arange(0.69,1.5,0.001)
    yline = xline*0.88+0.69
    ax.plot(xline,yline,color='black')
    #Titles, axes, legends
     ax.set_ylabel('U-V',fontsize = axisfont)
    ax.set_xlabel('V-J',fontsize = axisfont)
    ax.tick_params(labelsize = ticksize, size=ticks)
    ax.set_xlim(-0.5,2.5)
    ax.set_ylim(0,2.5)
    c=c+1

fig.tight_layout()
fig.savefig(figout + 'UVJ_bptlines.pdf')
plt.close(fig)

