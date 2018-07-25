
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


#Get the strings of each line
lines = ['4861','4959','5007','6563']#'6548','6583']

#Plot of Hbflux vs Hbflux of duplicates

minlim= 0.005
maxlim= 100

#Filter to get the good and back liens
flagidx =(fluxdata['4861_flag']==0)
goodflux = fluxdata[flagidx]
badflux = fluxdata[np.logical_not(flagidx)]
dupids = [item for item, count in collections.Counter(goodflux['OBJID']).items() if count > 1]
dupobjsarr = []
hbflux1 = []
hbflux2 = []
for obj in dupids:
    dupobjsarr.append(fluxdata[fluxdata.OBJID == obj])
for i in range(0,len(dupobjsarr)):
    hbflux1.append(divz(dupobjsarr[i].iloc[0]['4861_flux'],dupobjsarr[i].iloc[0]['4861_scale']))
    hbflux2.append(divz(dupobjsarr[i].iloc[1]['4861_flux'],dupobjsarr[i].iloc[1]['4861_scale']))

dupidsb = [item for item, count in collections.Counter(badflux['OBJID']).items() if count > 1]
dupobjsarrb = []
hbflux1b = []
hbflux2b = []
for obj in dupidsb:
    dupobjsarrb.append(fluxdata[fluxdata.OBJID == obj])
for i in range(0,len(dupobjsarrb)):
    hbflux1b.append(divz(dupobjsarrb[i].iloc[0]['4861_flux'],dupobjsarrb[i].iloc[0]['4861_scale']))
    hbflux2b.append(divz(dupobjsarrb[i].iloc[1]['4861_flux'],dupobjsarrb[i].iloc[1]['4861_scale']))

lw=0.25
mark='.'

fig,ax = plt.subplots(figsize=(8,7))
ax.scatter(hbflux1,hbflux2,color='blue',marker=mark,label='Galaxy with good fit')
ax.scatter(hbflux1b,hbflux2b,color='red',marker=mark,label='Galaxy with bad fit')
ax.plot((0,1000),(0,1000),color='black',ls='--')

#Titles, axes, legends
ax.set_title('H$\\beta$ Flux in Duplicate Measurements',fontsize = titlefont)
ax.legend(fontsize = legendfont)
ax.set_xlabel('H$\\beta$ Flux ($10^{-17}$ erg/s/$cm^2$)',fontsize = axisfont)
ax.set_ylabel('H$\\beta$ Flux ($10^-{17}$ erg/s/$cm^2$)',fontsize = axisfont)
ax.text(0.02,0.79,'1.49*MAD:     ' + str(round(mad_df['4861_mad'],2)),fontsize = textfont, transform=ax.transAxes)      
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(minlim,maxlim)
ax.set_ylim(minlim,maxlim)
ax.tick_params(labelsize = ticksize)
plt.show()
fig.savefig(figout + 'HB_dup_flux.pdf')
plt.close(fig)
