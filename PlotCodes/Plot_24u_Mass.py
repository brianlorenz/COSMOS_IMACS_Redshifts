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

#The location with the file for the filter data
filtdatapath = '/Users/blorenz/COSMOS/COSMOSData/all_c_hasinger.txt'
#Read in the data
filtdata = ascii.read(filtdatapath).to_pandas()
mips24data = filtdata[['id','mips24','emips24','Ks','eKs']]
mips24data = mips24data.rename(columns={'id':'OBJID'})

fluxdata = pd.merge(fluxdata,mips24data,on='OBJID',how='inner')
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



#Gets rid of objects with no mips measurement
filtmips = fluxdata['mips24']>0
filtks = fluxdata['Ks']>0
filtmuz = np.logical_and(filtmips,filtks)

emed24 = np.median(fluxdata[filtmips]['emips24'])
fluxdata.loc[np.logical_not(filtmips)]['mips24'] = emed24

#Finds the magnitudes of these bands:
fluxdata['mips24mag'] = -2.5*np.log10(fluxdata['mips24'])+25
fluxdata['Ksmag'] = -2.5*np.log10(fluxdata['Ks'])+25


#Plot the data with error bars
#Counter
c = 0
ylabel = ''



ms=12
lwbw=2

notbad = np.logical_not(baddata)

'''
fig,axarr = plt.subplots(1,3,figsize=(24,7),sharex=True,sharey=True)

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
    filt = np.logical_and(filt,filtmuz)
    #k24 = fluxdata[filt]['Ks'] - fluxdata[filt]['mips24']
    k24 = fluxdata[filt]['Ksmag'] - fluxdata[filt]['mips24mag']
    ax.plot(k24,fluxdata[filt]['UmV'],color=color,marker='o',ms=4,lw=0.5,ls='None')
    #ax.plot((-100,0.69),(1.3,1.3),color='black')
    #ax.plot((1.5,1.5),(2.01,100),color='black')
    #xline = np.arange(0.69,1.5,0.001)
    #yline = xline*0.88+0.69
    #ax.plot(xline,yline,color='black')
    #Titles, axes, legends
    ax.set_ylabel('U-V',fontsize = axisfont)
    ax.set_xlabel('Ks (mag) - [24]',fontsize = axisfont)
    ax.tick_params(labelsize = ticksize, size=ticks)
    ax.set_xlim(-3,3)
    ax.set_ylim(0,2.5)
    c=c+1

fig.tight_layout()
fig.savefig(figout + 'UV_Ks24.pdf')
plt.close(fig)
'''


cutvar = 'k24'
cutvarname = 'k24mag'
cutoff = 0.5
ycutname = 'Ks-[24]'


fluxdata['k24'] = fluxdata['Ksmag'] - fluxdata['mips24mag']
emipsmag = -2.5*(1/np.log(10))*divz(fluxdata['emips24'],fluxdata['mips24'])
ekmag = -2.5*(1/np.log(10))*divz(fluxdata['eKs'],fluxdata['Ks'])
fluxdata['ek24'] = np.sqrt(emipsmag**2+ekmag**2)


combinegoodlow = 1


fig,ax = plt.subplots(figsize=(8,7))
#axarr = np.reshape(axarr,6)


col = 'good'
filt=np.logical_not(baddata)
filt = np.logical_and(filt,filtmuz)
ax.errorbar(fluxdata[filt]['LMASS'],fluxdata[filt]['k24'],yerr=fluxdata[filt]['ek24'],color='black',marker='o',ms=4,lw=0.5,ls='None',label=None)
cutoff = np.median(fluxdata[filt]['av'])
binlow = fluxdata['av']<cutoff
binhigh = fluxdata['av']>=cutoff
filtlow = np.logical_and(filt,binlow)
filthigh = np.logical_and(filt,binhigh)
ax.scatter(fluxdata[filtlow]['LMASS'],fluxdata[filtlow]['k24'],c='red',marker='o',s=16,zorder=3,label='Av < '+ str(round(cutoff,2)))
ax.scatter(fluxdata[filthigh]['LMASS'],fluxdata[filthigh]['k24'],c='blue',marker='o',s=16,zorder=3,label='Av >= '+ str(round(cutoff,2)))
ax.legend()
ax.set_ylim(-3,3)
ax.set_xlim(8.95,10.05)
ax.set_xlabel('log(Stellar Mass) (M$_{sun}$)',fontsize = axisfont)
ax.set_ylabel('Ks-[24] (mag)',fontsize = axisfont)
#Titles, axes, legends
ax.tick_params(labelsize = ticksize, size=ticks)
 

fig.tight_layout()
fig.savefig(figout + 'Av_Ks24_LMASS.pdf')
plt.close(fig)

