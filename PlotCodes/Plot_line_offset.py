
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

#Line offsets against each other
#We want good Ha, HB, and OIII
flagidx = np.logical_and((fluxdata['6563_flag']==0), (fluxdata['6563_stddev']>1))
flagidx = np.logical_and(flagidx,(fluxdata['4861_flag']==0))
flagidx = np.logical_and(flagidx,(fluxdata['5007_flag']==0))
#flagidx = np.logical_and(flagidx,(fluxdata['4959_flag']==0))
#flagidx = np.logical_and(flagidx,(fluxdata['4959_stddev']>1))
flagidx = np.logical_and(flagidx,(fluxdata['4861_stddev']>1))
flagidx = np.logical_and(flagidx,(fluxdata['5007_stddev']>1))
goodflux = fluxdata[flagidx]
badflux = fluxdata[np.logical_not(flagidx)]

                                        
#ax.hist([h1,h2,h3,h4],color=['darkblue','mediumslateblue','dodgerblue','lightblue'],bins=bins,label=[lines[0],lines[1],lines[2],lines[3]])
fig,ax = plt.subplots(figsize=(13,7))



colors = goodflux['5007_mean']-5007*(1+goodflux['zcc'])
p1 = goodflux['6563_mean']-6562.80*(1+goodflux['zcc'])
p2 = goodflux['4861_mean']-4861.3*(1+goodflux['zcc'])
figure = ax.scatter(p1,p2,c=colors)

cb = fig.colorbar(figure)
ax.text(1.145,0.7,'O[III] 4959 Offset $(\AA)$',fontsize = axisfont, transform=ax.transAxes,rotation=90)      

#Titles, axes, legends
ax.set_ylabel('H$\\beta$ Offset $(\AA)$',fontsize = axisfont)
ax.set_title('Offset Comparison',fontsize = titlefont)
ax.set_xlabel('H$\\alpha$ Offset $(\AA)$',fontsize = axisfont)
ax.tick_params(labelsize = ticksize)
cb.ax.tick_params(labelsize = ticksize)
ax.set_xlim(-5,5)
ax.set_ylim(-5,5)
fig.tight_layout()
plt.show()
fig.savefig(figout + 'offset_by_gal.pdf')
plt.close(fig)
##Fit slope to this line and see if it is ratio of the wavelengths
