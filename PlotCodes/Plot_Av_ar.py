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

#Location of the Battisti data
bdata = '/Users/blorenz/COSMOS/COSMOSData/Battisti_data.csv'
#Read in the ew of the lines 
bat_df = ascii.read(bdata).to_pandas()

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

notbad = np.logical_not(baddata)

#Plot the data with error bars
#Counter
c = 0
ylabel = ''

ms=12
lwbw=2

paper=1

x = np.arange(0.1,89.9,0.1)
y = divz(0.01,np.sin(x*np.pi/180))/0.01
x2 = np.sin(x*np.pi/180)

if not paper:
    fig = plt.figure(figsize = (19.5,7))
    ax = fig.add_axes((0.15,0.15,0.315,0.8))
else: fig,ax = plt.subplots(figsize = (8,7))
    
filt = notbad
filt = np.logical_and(filt,fluxdata['ar']>0.0)
filt = np.logical_and(filt,fluxdata['av']>0.0)
c=0
color = 'blue'
yerru = fluxdata['dav2']
yerrd = fluxdata['dav1']
ax.errorbar(fluxdata[filt]['ar'],fluxdata[filt]['av'],yerr=np.array([yerrd[filt],yerru[filt]]),color='cornflowerblue',marker='o',ms=4,lw=0.5,ls='None',label='Av error >0.17 dex')
errfilt = yerru<0.17
errfilt = np.logical_and(filt,errfilt)
pcirdark = ax.errorbar(fluxdata[errfilt]['ar'],fluxdata[errfilt]['av'],yerr=np.array([yerrd[errfilt],yerru[errfilt]]),color='indigo',marker='o',ms=4,lw=0.5,ls='None',label='Av error <0.17  dex')
ax.set_xlabel('Axis Ratio',fontsize = axisfont)
ax.set_ylabel('Av (mag)',fontsize = axisfont)
ax.set_ylim(0,4.5)
ax.set_xlim(0.1,1)
q1 = np.percentile(fluxdata[filt]['ar'],25)
q2 = np.percentile(fluxdata[filt]['ar'],50)
q3 = np.percentile(fluxdata[filt]['ar'],75)
mr1 = (fluxdata[filt]['ar']<q1)
mr2 = np.logical_and(fluxdata[filt]['ar']>=q1,fluxdata[filt]['ar']<q2)
mr3 = np.logical_and(fluxdata[filt]['ar']>=q2,fluxdata[filt]['ar']<q3)
mr4 = (fluxdata[filt]['ar']>=q3)
mrs = np.array([mr1,mr2,mr3,mr4])
meds = np.array([np.median(fluxdata['av'][filt][i]) for i in mrs])
emeds = 1.49*np.array([np.median(np.abs(fluxdata['av'][filt][i]-np.median(fluxdata['av'][filt][i]))) for i in mrs])
        
bins = np.array([(0.1+q1)/2,(q1+q2)/2,(q2+q3)/2,(q3+1)/2])
msbw = 12
lwbw = 3
mew=3
offset=0.01
bcolor='black'
ax.errorbar(bins,meds,yerr=emeds,marker='o',ms=msbw,lw=lwbw,ls='None',markerfacecolor='None', markeredgecolor=bcolor,mew=mew,ecolor=bcolor,label='Median in Axis Ratio bin')
bat_df = bat_df[bat_df['col2']<0.8]
bAv = (bat_df['col2']/0.276)
ax.plot(bat_df['col1'],bAv,color='red',label='Galex Sample (Battisti+ 2017)',lw=3)
handles, labels = ax.get_legend_handles_labels()
hand = [handles[2],handles[1],handles[3],handles[0]]
if not paper: ax.legend(handles=hand,fontsize=axisfont-4,bbox_to_anchor=(1.1, 0.5))
else: ax.legend(handles=hand,fontsize=axisfont-8)
ax.tick_params(labelsize = ticksize, size=ticks)

fig.tight_layout()
fig.savefig(figout + 'Av_ar.pdf')
plt.close(fig)

