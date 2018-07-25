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

sig = 3
lines = ['4861','5007','4959','6583','6563_fix']
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

def scaleline(dataframe,line):
    newflux = divz(dataframe[line+'_flux'],dataframe[line+'_scale'])
    return newflux

Hagood = scaleline(goodflux,'6563_fix')
N2good = scaleline(goodflux,'6583')
O5007good = scaleline(goodflux,'5007')
O4959good = scaleline(goodflux,'4959')
O2good = scaleline(goodflux,'3727')
Hbgood = scaleline(goodflux,'4861')
R23numgood = O5007good+O4959good+O2good
goodR23 = divz(R23numgood,Hbgood)
goodHaNII = divz(N2good,Hagood)


Halow = scaleline(lowflux,'6563_fix')
N2low = scaleline(lowflux,'6583')
O5007low = scaleline(lowflux,'5007')
O4959low = scaleline(lowflux,'4959')
O2low = scaleline(lowflux,'3727')
Hblow = scaleline(lowflux,'4861')
R23numlow = O5007low+O4959low+O2low
lowR23 = divz(R23numlow,Hblow)
lowHaNII = divz(N2low,Halow)
                    
errHa = mad_df['6563_fix_mad'][0]
errN2 = mad_df['6583_fix_mad'][0]
errHb = mad_df['4861_mad'][0]
errO5007 = mad_df['5007_mad'][0]
errO4959 = mad_df['4959_mad'][0]
errO2 = mad_df['3727_mad'][0]
errR23num = errO5007+errO4959+errO2

eHaNIIgood = np.sqrt((divz(1,N2good) * errHa)**2 + (divz(-Hagood,(N2good**2)) * errN2)**2)
eR23good = np.sqrt((divz(1,R23numgood) * errHb)**2 + (divz(-Hbgood,(R23numgood**2)) * errR23num)**2)
eHaNIIlow = np.sqrt((divz(1,N2low) * errHa)**2 + (divz(-Halow,(N2low**2)) * errN2)**2)
eR23low = np.sqrt((divz(1,R23numlow) * errHb)**2 + (divz(-Hblow,(R23numlow**2)) * errR23num)**2)

#Not done yet, need errors

lw=0.125
mark='.'
ms=6

fig,ax = plt.subplots(figsize=(13,6))
ax.errorbar(np.log10(lowHaNII),np.log10(lowR23),eHaNIIlow,eR23low,ls='None',lw=lw,ms=ms,color='grey',marker=mark,label='Low in some line')
ax.errorbar(np.log10(goodHaNII),np.log10(goodR23),eHaNIIgood,eR23good,ls='None',lw=lw,color='blue',ms=ms,marker=mark,label='Above 3-sigma')

ax.plot((-1.1,-1.1),(-10,10),c='black',ls='--')
ax.plot((-1.3,-1.3),(-10,10),c='black',ls='--')


#Titles, axes, legends
ax.set_title('R23 vs NII/Ha',fontsize = titlefont)
ax.legend(fontsize = legendfont)
ax.set_xlabel('log10(NII/Ha)',fontsize = axisfont)
ax.set_ylabel('log10(R23)',fontsize = axisfont)
#ax.set_xscale("log")
#ax.set_yscale("log")
ax.set_xlim(-2.5,0)
ax.set_ylim(-0.2,1.2)
ax.tick_params(labelsize = ticksize)
plt.show()
fig.savefig(figout + 'R23.pdf')
plt.close(fig)

