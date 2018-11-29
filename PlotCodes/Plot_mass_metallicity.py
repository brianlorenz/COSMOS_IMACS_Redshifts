#Creates a BPT diagram for all objects, and a second figure that shows objects for which single lines are low
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import sys, os, string
import pandas as pd
from astropy.io import fits
import collections
from cycler import cycler
from matplotlib import cm


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
lines = ['6563_fix','6583_fix','4861','5007']
Ha = lines[0]
N2 = lines[1]
Hb = lines[2]
O3 = lines[3]



#fig2,axarr2 = plt.subplots(2,2,figsize=(15,12))
#ax1,ax2,ax3,ax4 = axarr2[0,0],axarr2[0,1],axarr2[1,0],axarr2[1,1]

#Takes the dataframe and the four lines to combine into the bpt
def getbpt(pd_df,err_df,Ha,Hb,N2,O3):
    errHa = err_df[Ha]
    errN2 = err_df[N2]
    errHb = err_df[Hb]
    errO3 = err_df[O3]
    #Divide by the scale to calibrate the flux
    calHa = divz(pd_df[Ha+'_flux'],pd_df[Ha+'_scale'])
    calHb = divz(pd_df[Hb+'_flux'],pd_df[Hb+'_scale'])
    calN2 = divz(pd_df[N2+'_flux'],pd_df[N2+'_scale'])
    calO3 = divz(pd_df[O3+'_flux'],pd_df[O3+'_scale'])
    #Find the ratios
    Harat = np.log10(divz(calN2,calHa))
    Hbrat = np.log10(divz(calO3,calHb))
    #Find the errors
    eHarat = (1/np.log(10))*divz(calHa,calN2)*np.sqrt((divz(1,calHa) * errN2)**2 + (divz(-calN2,(calHa**2)) * errHa)**2)
    eHbrat = (1/np.log(10))*divz(calHb,calO3)*np.sqrt((divz(1,calHb) * errO3)**2 + (divz(-calO3,(calHb**2)) * errHb)**2)
    return (Harat,Hbrat,eHarat,eHbrat)

#Plotting parameters
ms = 6
lw=0.5
mark='o'

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
plotframe['OBJID'] = fluxdata['OBJID']
plotframe['Mask'] = fluxdata['Mask']

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
    bpts = getbpt(fluxdata[filt],err_df[filt],Ha,Hb,N2,O3)
    plotframe[colname+'x'] = bpts[0]
    plotframe[colname+'y'] = bpts[1]
    plotframe[colname+'ex'] = bpts[2]
    plotframe[colname+'ey'] = bpts[3]


#Find the duplicates
dupids = [item for item, count in collections.Counter(fluxdata[allgood]['OBJID']).items() if count > 1]
dupobjsarr = []
for obj in dupids:
    dupobjsarr.append(fluxdata[fluxdata.OBJID == obj].OBJID)

    
#Line that defines the agn cutoff
xline = np.arange(-3.0,0.469,0.001)
yline = 0.61/(xline-0.47)+1.19 #Kewley (2001)
xlineemp = np.arange(-3.0,0.049,0.001)
ylineemp = 0.61/(xlineemp-0.05)+1.3 #Kauffman (2003)





#If duplicates, set to 1, otherwise 0
dup = 1




#Make the figures
fig2,axarr2 = plt.subplots(1,1,figsize=(8,7),sharex=True,sharey=True)

HaNII = 1

#Plot the data with error bars
c = 0
for i in range(0,2):
    ax=axarr2
    if c==0:
        col = 'good'
        filt = allgood
    elif c==1:
        break
        col = 'low'
        filt = somelow
    else:
        col = 'bad'
        filt = baddata
    #colors = fluxdata[filt]['LMASS']
    p = plotframe
    if HaNII:
        ax.errorbar(fluxdata[filt]['LMASS'],p[filt][col+'x'],yerr=p[filt][col+'ex'],c='cornflowerblue',ls='None',ms=ms,marker=mark,lw=lw,zorder=1,label=None)
    else:
        ax.errorbar(fluxdata[filt]['LMASS'],p[filt][col+'y'],yerr=p[filt][col+'ey'],c='cornflowerblue',ls='None',ms=ms,marker=mark,lw=lw,zorder=1,label=None)
    #colordata = ax.scatter(p[filt][col+'x'],p[filt][col+'y'],c=colors,marker=mark,zorder=3,label=None)
    #Titles, axes, legends
    q1 = 9.25
    q2 = 9.5
    q3 = 9.73
    mr1 = (fluxdata[filt]['LMASS']<q1)
    mr2 = np.logical_and(fluxdata[filt]['LMASS']>=q1,fluxdata[filt]['LMASS']<q2)
    mr3 = np.logical_and(fluxdata[filt]['LMASS']>=q2,fluxdata[filt]['LMASS']<q3)
    mr4 = (fluxdata[filt]['LMASS']>=q3)
    mrs = np.array([mr1,mr2,mr3,mr4])
    if HaNII: meds = np.array([np.median(p[col+'x'][filt][i]) for i in mrs])
    else: meds = np.array([np.median(p[col+'y'][filt][i]) for i in mrs])
    emeds = 1.49*np.array([np.median(np.abs(fluxdata['LMASS'][filt][i]-np.median(fluxdata['LMASS'][filt][i]))) for i in mrs])
    bins = np.array([(9+q1)/2,(q1+q2)/2,(q2+q3)/2,(q3+10)/2])
    msbw = 12
    lwbw = 3
    mew=3
    offset=0.01
    bcolor='black'
    #ax.errorbar(bins,meds,yerr=emeds,marker='o',ms=msbw,lw=lwbw,ls='None',markerfacecolor='None', markeredgecolor=bcolor,mew=mew,ecolor=bcolor,label='Median in mass bin')
    if c==0:
        if HaNII: ax.set_ylabel('log(N[II] 6583 / H$\\alpha$)',fontsize = axisfont)
        else: ax.set_ylabel('log(O[III] 5007 / H$\\beta$)',fontsize = axisfont)
    ax.set_xlabel('log(Stellar Mass) (M$_\odot)$',fontsize = axisfont)
    ax.tick_params(labelsize = ticksize, size=ticks)
    ax.set_xlim(8.95,10.05)
    if HaNII: ax.set_ylim(-1.5,0.2)
    else: ax.set_ylim(-1,0.8)
    ax.plot(0,0,ls='-',color='black',alpha=0.5,label='Template')
    ax.plot(0,0,ls='-',color='cornflowerblue',label='Observed Spectrum')
    if c==0: ax.legend(fontsize=legendfont)
    c = c+1
fig2.tight_layout()
if HaNII: fig2.savefig(figout + 'Mass_N2Ha.pdf')
else: fig2.savefig(figout + 'Mass_O3Hb.pdf')
plt.close(fig2)
    

