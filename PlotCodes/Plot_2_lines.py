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

#Merge our data with the UVISTA catalog
#The location of the muzzin et al data:
mdatapath = '/Users/blorenz/COSMOS/muzzin_data/UVISTA_final_colors_sfrs_v4.1.dat'

#Read in the mad of the lines 
mad_df = ascii.read(maddatapath).to_pandas()


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

#Plot of Haflux vs Hgflux of duplicates
l1 = sys.argv[1]
l2 = sys.argv[2]
sig = int(sys.argv[3])


minlim= 0.005
maxlim= 1000

lines = [l1,l2]

def findgoodfits(pd_df=fluxdata,lines=lines,sig=sig):
    goodidxs = [(pd_df[line+'_flag'] == 0) for line in lines]
    lowidxs = [(divz(pd_df[line+'_flux'],pd_df[line+'_scale']) < sig*mad_df[line+'_mad'][0]) for line in lines]
    goodidx = np.logical_and.reduce(goodidxs)
    lowidx = np.logical_or.reduce(lowidxs)
    badflux = pd_df[np.logical_not(goodidx)]
    lowflux = pd_df[np.logical_and(goodidx,lowidx)]
    goodflux = pd_df[np.logical_and(goodidx,np.logical_not(lowidx))]
    return goodflux,lowflux,badflux


goodflux,lowflux,badflux = findgoodfits()

yerr = mad_df[l2+'_mad'][0]
xerr = mad_df[l1+'_mad'][0]


xdata = divz(goodflux[l1+'_flux'],goodflux[l1+'_scale'])
ydata = divz(goodflux[l2+'_flux'],goodflux[l2+'_scale'])
xdatalow = divz(lowflux[l1+'_flux'],lowflux[l1+'_scale'])
ydatalow = divz(lowflux[l2+'_flux'],lowflux[l2+'_scale'])
   
lw=0.25
mark='.'
ms=6

fig,ax = plt.subplots(figsize=(8,7))
ax.errorbar(xdatalow,ydatalow,xerr,yerr,ls='None',lw=lw,ms=ms,color='grey',marker=mark)
ax.errorbar(xdata,ydata,xerr,yerr,ls='None',lw=lw,color='blue',ms=ms,marker=mark)

ax.plot((0,1000),(0,1000),color='black',ls='--')

#Titles, axes, legends
ax.set_title(l2+ ' Flux vs ' + l1 + ' Flux',fontsize = titlefont)
ax.legend(fontsize = legendfont)
ax.set_xlabel(l1  +' Flux ($10^{-17}$ erg/s/$cm^2$)',fontsize = axisfont)
ax.set_ylabel(l2 + ' Flux ($10^-{17}$ erg/s/$cm^2$)',fontsize = axisfont)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(minlim,maxlim)
ax.set_ylim(minlim,maxlim)
ax.tick_params(labelsize = ticksize)
plt.show()
fig.savefig(figout + 'Flux_' + l2 + '_' + l1 + '.pdf')
plt.close(fig)
