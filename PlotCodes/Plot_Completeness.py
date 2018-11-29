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

#Read all of the data
alldatapath = '/Users/blorenz/COSMOS/COSMOSData/all_c_hasinger.txt'
alldata = ascii.read(alldatapath).to_pandas()

#Read the muzzin sfrs
muzdatapath = '/Users/blorenz/COSMOS/muzzin_data/UVISTA_final_colors_sfrs_v4.1.dat'
muzdata = ascii.read(muzdatapath).to_pandas()

alldata = pd.merge(alldata,muzdata,left_on='id',right_on='ID')

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

lines=['6563_fix']

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



#Plot the data with error bars
#Counter
c = 0
ylabel = ''

fluxids = []
for i in range(0,len(fluxdata)):
    fluxids.append(fluxdata.fluxfile.iloc[i][4:])

for string in fluxids:
    alldata.drop(alldata[alldata.ImageName == 'cor_'+string].index[0],inplace=True)
nostar = (alldata.Star==0)
    


ms=6
lwbw=2

notbad = np.logical_not(baddata)

paper = 1

if not paper:
    fig = plt.figure(figsize = (15,7))
    ax = fig.add_axes((0.15,0.15,0.42,0.8))
else:
    fig,ax = plt.subplots(figsize=(8,7))

    
c=0
color = 'blue'
ax.plot(fluxdata[allgood]['LMASS'],fluxdata[allgood]['UmV'],color=color,marker='o',ms=ms,lw=0.5,ls='None',label='Significant H$\\alpha$ Detection')
ax.plot(fluxdata[somelow]['LMASS'],fluxdata[somelow]['UmV'],color=color,marker='v',mfc='None',ms=ms,lw=0.5,ls='None',label='Upper Limit on H$\\alpha$')
ax.plot(fluxdata[baddata]['LMASS'],fluxdata[baddata]['UmV'],color=color,marker='o',mfc='None',ms=ms,lw=0.5,ls='None',label='No H$\\alpha$ Measurement',mew=1.2)
ax.plot(alldata[nostar]['LMASS'],alldata[nostar]['UmV'],color='red',marker='o',mfc='None',ms=ms,lw=0.5,ls='None',label='No Redshift Measurement',mew=1.2)
if paper:
    ax.set_ylim(-0.2,2.2)
else:
    ax.set_ylim(0.1,2.2)
ax.set_xlim(8.95,10.05)
ax.set_xlabel('log(Stellar Mass) (M$_\odot)$',fontsize = axisfont)
ax.set_ylabel('U-V (mag)',fontsize = axisfont)
if not paper:
    ax.legend(fontsize=axisfont-2,bbox_to_anchor=(1.01, 0.5))
else:
    ax.legend(loc=4,fontsize=axisfont-8)
ax.tick_params(labelsize = ticksize, size=ticks)

fig.tight_layout()
fig.savefig(figout + 'Completeness.pdf')
plt.close(fig)

