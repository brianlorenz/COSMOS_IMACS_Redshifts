#Creates a BPT diagram for all objects, and a second figure that shows objects for which single lines are low
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

#Location of the equivalent width data
ewdata = '/Users/blorenz/COSMOS/COSMOSData/lineew.txt'
#Read in the ew of the lines 
ew_df = ascii.read(ewdata).to_pandas()

#The location to store the scale and its stddev of each line
qualdatapath = '/Users/blorenz/COSMOS/COSMOSData/dataqual.txt'
#Read in the scale of the lines 
dataqual = ascii.read(qualdatapath).to_pandas()

#File with the error array
errdatapath = '/Users/blorenz/COSMOS/COSMOSData/errs.txt'
#Read in the scale of the lines 
err_df = ascii.read(errdatapath,data_start=1,header_start=0,format='csv').to_pandas()


#Read the datafile:
fluxdata = ascii.read(fluxdatapath).to_pandas()

#Check if bpt correlates with stellar mass
#The location of the muzzin et al data:
mdatapath = '/Users/blorenz/COSMOS/muzzin_data/UVISTA_final_colors_sfrs_v4.1.dat'
#Read in the muzzin data
mdata = ascii.read(mdatapath).to_pandas()
mdata = mdata.rename(columns={'ID':'OBJID'})

fluxdata = pd.merge(fluxdata,mdata)


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


#BPT
#Strings of all of the lines needed for bpt
lines = ['4861','5007']
Hb = lines[0]
O3 = lines[1]



#fig2,axarr2 = plt.subplots(2,2,figsize=(15,12))
#ax1,ax2,ax3,ax4 = axarr2[0,0],axarr2[0,1],axarr2[1,0],axarr2[1,1]

#Takes the dataframe and the four lines to combine into the bpt
def getO3Hb(pd_df,err_df,N2,O3):
    errHb = err_df[Hb]
    errO3 = err_df[O3]
    #Divide by the scale to calibrate the flux
    calHb = divz(pd_df[Hb+'_flux'],pd_df[Hb+'_scale'])
    calO3 = divz(pd_df[O3+'_flux'],pd_df[O3+'_scale'])
    #Find the ratios
    Hbrat = np.log10(divz(calO3,calHb))
    #Find the errors
    eHbrat = (1/np.log(10))*divz(calHb,calO3)*np.sqrt((divz(1,calHb) * errO3)**2 + (divz(-calO3,(calHb**2)) * errHb)**2)
    return (Hbrat,eHbrat)

#Plotting parameters
ms = 3
lw=0.5
mark='o'

d = {'True': True, 'False': False}

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


plotframe = pd.DataFrame()
plotframe['fluxfile'] = fluxdata['fluxfile']

#Get the x and y axes for the good, low, and bad data
for i in range(0,3):
    if i==0:
        filt = allgood
        colname = 'good'
    elif i==1:
        filt = somelow
        colname = 'low'
    else:
        filt = baddata
        colname = 'bad'
    data = getO3Hb(fluxdata[filt],err_df[filt],Hb,O3)
    plotframe[colname+'y'] = data[0]
    plotframe[colname+'ey'] = data[1]

#Line that defines the agn cutoff
xline = np.arange(-3.0,0.49,0.001)
yline = 0.61/(xline-0.5)+1.3 #Groves and Kewley


#Make the figures
fig,axarr = plt.subplots(1,3,figsize=(24,7),sharex=True,sharey=True)

#Plot the data with error bars
c = 0
for ax in axarr:
    if c==0:
        col = 'good'
        color = 'blue'
        filt = allgood
    elif c==1:
        col = 'low'
        color = 'orange'
        filt = somelow
    else:
        col = 'bad'
        color = 'red'
        filt = baddata
    p = plotframe
    ssfr_gy = divz(fluxdata[filt]['SFR_tot'],10**fluxdata[filt]['LMASS'])*10**9
    ax.errorbar(ssfr_gy,p[filt][col+'y'],yerr=p[filt][col+'ey'],color=color,ls='None',ms=ms,marker=mark,lw=lw)
    #Titles, axes, legends
    if c==0:
        ax.set_ylabel('log10(O[III] 5007 / H$\\beta$)',fontsize = axisfont)
    ax.set_xlabel('sSFR ($Gyr^{-1}$)',fontsize = axisfont)
    ax.tick_params(labelsize = ticksize, size=ticks)
    ax.set_xlim(0.01,100)
    ax.set_xscale('log')
    ax.set_ylim(-1.25,1.5)
    #Plot the bpt line
    c = c+1
fig.tight_layout()
fig.savefig(figout + 'O3Hb_ssfr.pdf')
plt.close(fig)


