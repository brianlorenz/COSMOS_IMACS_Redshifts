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


#File with the error array
errreddatapath = '/Users/blorenz/COSMOS/COSMOSData/errs_red.txt'
#Read in the scale of the lines 
err_dfred = ascii.read(errreddatapath,data_start=1,header_start=0,format='csv').to_pandas()

#Read the datafile:
fluxdata = ascii.read(fluxdatapath).to_pandas()

#File with the structural properties
spropdatapath = '/Users/blorenz/COSMOS/COSMOSData/struct_prop.txt'
#Read in the scale of the lines 
sprop_df = ascii.read(spropdatapath).to_pandas()
sprop_df = sprop_df.rename(columns={'id':'OBJID'})
fluxdata = pd.merge(fluxdata,sprop_df)

#The location with the file for the filter data
filtdatapath = '/Users/blorenz/COSMOS/COSMOSData/all_c_hasinger.txt'
#Read in the data
filtdata = ascii.read(filtdatapath).to_pandas()
cordata = filtdata[['id','Ks','eKs','Ks_tot','eKs_tot']]
cordata = cordata.rename(columns={'id':'OBJID'})

fluxdata = pd.merge(fluxdata,cordata,on='OBJID',how='inner')
fluxdata = fluxdata.drop_duplicates()
fluxdata = fluxdata.reset_index()

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

err = err_dfred['6563_fix_u']
errSFR = err*10**-17*4*np.pi*((cosmo.luminosity_distance(fluxdata['zcc'])*3.086*10**24)**2)*10**-41.27
Halum = fluxdata['6563_fix_flux_red']*10**-17*4*np.pi*((cosmo.luminosity_distance(fluxdata['zcc'])*3.086*10**24)**2)*(fluxdata.Ks_tot/fluxdata.Ks)
fluxdata['SFR'] = Halum*10**-41.27 #Kennicutt and Evans (2012)
fluxdata['SFR'] = fluxdata['SFR'].astype(float)


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


combinemass = 1


filtSFR = fluxdata['SFR']<10000000

ms=12
lwbw=2

notbad = np.logical_not(baddata)

ssfr = 1

fig,ax = plt.subplots(figsize=(8,7))
c=0

for w in range(0,2):
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
    filt = np.logical_and(filt,filtSFR)
    xdata = fluxdata[filt]['LMASS']
    ydata = np.log10(fluxdata[filt]['SFR'])
    #yerr = 1/np.log(10)*(errSFR[filt].astype(float)/fluxdata[filt]['SFR'])
    ax.set_xlabel('log(Stellar Mass) (M$_{sun})$',fontsize = axisfont)
    ax.set_ylabel('log(SFR) (M$_{sun})$/yr)',fontsize = axisfont)
    if ssfr:
        ydata = np.log10(divz(10**ydata,10**fluxdata[filt]['LMASS']))
        ax.set_ylabel('log(sSFR) (yr$^{-1}$)',fontsize = axisfont)
        if c==0:
            coeff = np.polyfit(xdata,ydata,1)
            x = np.arange(9,10,0.001)
            fit = coeff[0]*x+coeff[1]
            ax.plot(x,fit,c='red',lw=4)
            ax.text(0.665,0.93,'Slope: ' + str(np.round(coeff[0],2)),fontsize = axisfont, transform=ax.transAxes)
    ax.plot(xdata,ydata,color=color,marker='o',ms=4,lw=0.5,ls='None')
    ax.tick_params(labelsize = ticksize, size=ticks)
    ax.set_xlim(8.95,10.05)
    ax.set_ylim(-3,1.6)
    if ssfr:
        ax.set_ylim(-12,-8)
    c=c+1

fig.tight_layout()
if ssfr: fig.savefig(figout + 'sSFR_Mass.pdf')
else: fig.savefig(figout + 'SFR_Mass.pdf')
plt.close(fig)
