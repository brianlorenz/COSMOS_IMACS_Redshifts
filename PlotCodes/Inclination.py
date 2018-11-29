#Examines 24um flux for our objects to investigate the location of dust
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




#Plot the data with error bars
#Counter
c = 0
ylabel = ''

ms=12
lwbw=2

x = np.arange(0.1,89.9,0.1)
y = divz(0.01,np.sin(x*np.pi/180))/0.01
x2 = np.sin(x*np.pi/180)


fig,ax = plt.subplots(figsize=(8,7))

c=0
color = 'blue'
ax.plot(x2,y,color='blue',marker='o',ms=4,lw=0.5,ls='None')
ax.set_xlabel('Axis Ratio',fontsize = axisfont)
ax.set_ylabel('Multiples of dust disk',fontsize = axisfont)
ax.set_ylim(0,20)
#ax.legend(loc=4,fontsize=axisfont-8)
ax.tick_params(labelsize = ticksize, size=ticks)

fig.tight_layout()
fig.savefig(figout + 'Inclination.pdf')
plt.close(fig)

