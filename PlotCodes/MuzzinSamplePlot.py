#Plots the sample in our mass and redshift range from Muzzin's catalog

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import sys, os, string
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy import units as u

#Folder to save the figures
figout = '/Users/blorenz/COSMOS/Reports/2018/Images/'

#The location with the file for all of our data
ourdatapath = '/Users/blorenz/COSMOS/COSMOSData/all_c_hasinger.txt'

#Merge our data with the UVISTA catalog
#The location of the muzzin et al data:
mdatapath = '/Users/blorenz/COSMOS/muzzin_data/UVISTA_final_colors_sfrs_v4.1.dat'

#Read in all of our data
ourdata = ascii.read(ourdatapath).to_pandas()

#Read in the muzzin data
mdata = ascii.read(mdatapath).to_pandas()
zidx = np.logical_and(mdata['z_peak'] > 0.3, mdata['z_peak'] < 0.4)
midx = np.logical_and(mdata['LMASS'] > 9, mdata['LMASS'] < 10)
idx = np.logical_and(zidx,midx)

total = mdata[idx]

#Fontsizes for plotting
axisfont = 18
ticksize = 16
titlefont = 24
legendfont = 16
textfont = 16
mark = 'o'
mark2 = 'o'
marks = 4.0
mark2s = 4.0
c1 = 'C0'
c2 = 'gold'
a2 = 1

fig,ax = plt.subplots(figsize=(8,7))
ax.plot(total.z_peak,total.LMASS,color=c2,markersize = mark2s,marker=mark2,ls='None',alpha=a2)
for i in range(0,len(ourdata)):
    plotobj = total[total.ID == ourdata.iloc[i]['OBJID']]
    ax.plot(plotobj.z_peak,plotobj.LMASS,color=c1,markersize = marks,marker=mark,ls='None')

h1, = ax.plot((0,1000),(0,1000),color=c1,marker=mark,ls='None',markersize = marks)
h2, = ax.plot((0,1000),(0,1000),color=c2,marker=mark2,ls='None',markersize = mark2s,alpha=a2)

#Titles, axes, legends
#ax.set_title('',fontsize = titlefont)
ax.legend(bbox_to_anchor=(0.6, 1.165),handles=[h1,h2],labels=['Sampled','Catalog'],fontsize = legendfont)
ax.set_xlabel('Photometric Redshift',fontsize = axisfont)
ax.set_ylabel('Log(Stellar Mass ($M_{sun}$)) ',fontsize = axisfont)
ax.set_xlim(0.295,0.405)
ax.set_ylim(8.95,10.05)
ax.tick_params(labelsize = ticksize)
plt.show()
fig.savefig(figout + 'muzzin_sample_plot.pdf')
