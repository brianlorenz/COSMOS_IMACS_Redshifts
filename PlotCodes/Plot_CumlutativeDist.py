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

print sys.argv
combine = int(sys.argv[1])

if not combine:
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
plotdata = fluxdata['av']
ylabel = ''
savename = 'CumuDist'

#fluxdata['n']=np.log10(fluxdata['n'])
#fluxdata['SMD']=np.log10(divz(fluxdata['LMASS'],(4*np.pi*fluxdata['re_kpc']**2)))

cutvar = fluxdata['LMASS']
cutvarname = 'LMASS'
cutoff = 9.5
ycutname = 'LMASS'


ms=12
lwbw=2

notbad = np.logical_not(baddata)


for ax in axarr:
    if c in [0,1,2]:
        cutfilt = cutvar<cutoff
    else:
        cutfilt = cutvar>=cutoff
    if c in [0,3]:
        col = 'good'
        filt = allgood
        if combine: filt = notbad
        color='blue'
    elif c in [1,4]:
        col = 'low'
        filt = somelow
        if combine: filt = (fluxdata.OBJID < 0)
        color='orange'
    else:
        col = 'bad'
        filt = baddata
        color='red'
    #Titles, axes, legends
    acount = 0
    filttype = (plotdata>-98.9)
    if c==0:
        ax.set_ylabel(ylabel + ' ' + ycutname + ' < ' + str(cutoff),fontsize = axisfont)
    if c==3:
        ax.set_ylabel(ylabel + ' ' + ycutname + ' >= ' + str(cutoff),fontsize = axisfont)
    
    ax.set_xlabel('Av (mag)',fontsize = axisfont)
    ax.tick_params(labelsize = ticksize, size=ticks)
    filters = np.logical_and(filt,cutfilt)
    filters = np.logical_and(filters,filttype)
    ydist = np.arange(len(plotdata[filters]))/float(len(plotdata[filters]))
    xdist = np.sort(plotdata[filters])
    ax.text(0.81,0.02,'N = ' + str(np.round(len(xdist),3)),fontsize = axisfont-2, transform=ax.transAxes)
    print len(xdist)
    ax.plot(xdist,ydist,color=color,lw=2)
    c = c+1
    if (combine and (c in [1,4])): c = c+1
    

fig.tight_layout()
if combine: fig.savefig(figout + 'Av_'+savename+'_'+cutvarname+'_combine.pdf')
else:fig.savefig(figout + 'Av_'+savename+'_'+cutvarname+'.pdf')
plt.close(fig)

