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
def getbpt(pd_df,err_df,Ha,Hb,N2,O3,p=0):
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
    Harat = divz(calN2,calHa)
    Hbrat = divz(calO3,calHb)
    #Find the errors
    eHarat = divz(calHa,calN2)*np.sqrt((divz(1,calHa) * errN2)**2 + (divz(-calN2,(calHa**2)) * errHa)**2)
    eHbrat = divz(calHb,calO3)*np.sqrt((divz(1,calHb) * errO3)**2 + (divz(-calO3,(calHb**2)) * errHb)**2)
    #Get their logs
    lHarat = np.log10(Harat)
    lHbrat = np.log10(Hbrat)
    #Get the log errors
    elHaratu = np.log10(Harat+eHarat)-lHarat
    ldiff = Harat-eHarat
    elHaratd = lHarat-np.log10(ldiff)
    elHaratd.loc[elHaratd<0] = 4
    elHaratd.fillna(4,inplace=True)
    print elHaratd.iloc[212]
    elHbratu = np.log10(Hbrat+eHbrat)-lHbrat
    ldiff = Hbrat-eHbrat
    elHbratd = lHbrat-np.log10(ldiff)
    elHbratd.loc[elHbratd<0] = 4
    elHbratd.fillna(4,inplace=True)
    return (lHarat,lHbrat,elHaratu,elHaratd,elHbratu,elHbratd)

#Plotting parameters
ms = 3
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
    if i==0: bpts = getbpt(fluxdata[filt],err_df[filt],Ha,Hb,N2,O3,p=1)
    else: bpts = getbpt(fluxdata[filt],err_df[filt],Ha,Hb,N2,O3)
    plotframe[colname+'x'] = bpts[0]
    plotframe[colname+'y'] = bpts[1]
    plotframe[colname+'exu'] = bpts[2]
    plotframe[colname+'exd'] = bpts[3]
    plotframe[colname+'eyu'] = bpts[4]
    plotframe[colname+'eyd'] = bpts[5]


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
if dup:
    fig,dupax = plt.subplots(1,1,figsize=(8.5,7),sharex=True,sharey=True)
    axarr = [1,2]
else: fig,axarr = plt.subplots(1,3,figsize=(24,7),sharex=True,sharey=True)
fig2,axarr2 = plt.subplots(1,1,figsize=(9,7))#,sharex=True,sharey=True)

#Color wheel if looking at duplicates
prop_iter = iter(plt.rcParams['axes.prop_cycle'])

colorsmask=['r','b','g','y','orange','purple','cyan','lime','gold','brown','pink']
masks=['a6','b6','d6','e6','f6','g6','h6','i6','b7','e7','j7']



#Plot the data with error bars
#Counter
c = 0
for ax in axarr:
    if ax==2:
        break
    if c==0:
        col = 'good'
        color = 'blue'
    elif c==1:
        col = 'low'
        color = 'orange'
    else:
        col = 'bad'
        color = 'red'
    p = plotframe
    if (c==0 and dup):
        ax = dupax
        '''
        for j in range(0,len(colorsmask)):
            ax.plot(-100,-100,color=colorsmask[j],ms=ms,marker=mark,label=masks[j])
        ax.legend()
        '''
    if dup: ax.errorbar(p[col+'x'],p[col+'y'],xerr=np.array([p[col+'exd'],p[col+'exu']]),yerr=np.array([p[col+'eyd'],p[col+'eyu']]),color='grey',ls='None',ms=ms,marker=mark,lw=lw,label=None,zorder=1)
    else: ax.errorbar(p[col+'x'],p[col+'y'],xerr=np.array([p[col+'exd'],p[col+'exu']]),yerr=np.array([p[col+'eyd'],p[col+'eyu']]),color=color,ls='None',ms=ms,marker=mark,lw=lw,label=None)
    #Titles, axes, legends
    if c==0:
        ax.set_ylabel('log(O[III] 5007 / H$\\beta$)',fontsize = axisfont)
    ax.set_xlabel('log(N[II] 6583 / H$\\alpha$)',fontsize = axisfont)
    ax.tick_params(labelsize = ticksize, size=ticks)
    ax.set_xlim(-2,1)
    ax.set_ylim(-1.2,1.5)
    #Plot the bpt line
    ax.plot(xline,yline,color='black',lw=1,ls='--',label='Theoretical, Kewley+ (2001)')
    ax.plot(xlineemp,ylineemp,color='black',lw=1,ls='-',label='Empirical, Kauffmann+ (2003)')
    ax.text(0.06,0.18,'Star-Forming',fontsize = textfont+6, transform=ax.transAxes)
    ax.text(0.7,0.66,'AGN',fontsize = textfont+6, transform=ax.transAxes)
    pdup = p.iloc[dupobjsarr[0].index]
    ax.errorbar(pdup[col+'x'].iloc[1],pdup[col+'y'].iloc[1],xerr=np.array(pdup[col+'exd'].iloc[1],pdup[col+'exu'].iloc[1]),yerr=np.array(pdup[col+'eyd'].iloc[1],pdup[col+'eyu'].iloc[1]),color='blue',ls='-',ms=ms,marker=mark,lw=1,zorder=3,label='Duplicate Observation')
    if c==0: ax.legend(fontsize=legendfont)
    #Draw lines between duplicates:
    if (c==0 and dup):
            for i in range(0,len(dupobjsarr)):
                pdup = p.iloc[dupobjsarr[i].index]
                
                ax.plot(pdup[col+'x'],pdup[col+'y'],color='red',ls='-',ms=ms,marker=mark,lw=1,zorder=2,label=None)
                print pdup
                
                ax.errorbar(pdup[col+'x'].iloc[0],pdup[col+'y'].iloc[0],xerr=np.array(pdup[col+'exu'].iloc[0],pdup[col+'exd'].iloc[0]),yerr=np.array(pdup[col+'eyd'].iloc[0],pdup[col+'eyu'].iloc[0]),color='blue',ls='-',ms=ms,marker=mark,lw=1,zorder=3,label=None)
                ax.errorbar(pdup[col+'x'].iloc[1],pdup[col+'y'].iloc[1],xerr=np.array(pdup[col+'exu'].iloc[1],pdup[col+'exd'].iloc[1]),yerr=np.array(pdup[col+'eyd'].iloc[1],pdup[col+'eyu'].iloc[1]),color='blue',ls='-',ms=ms,marker=mark,lw=1,zorder=3,label=None)
    c = c+1
fig.tight_layout()
if dup: fig.savefig(figout + 'bpt_dup.pdf')
else: fig.savefig(figout + 'bpt.pdf')
plt.close(fig)


#Plot the data with error bars
c = 0
while 1<2:
    ax = axarr2
#for ax in axarr2:
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
    colors = fluxdata[filt]['LMASS']
    p = plotframe
    ax.errorbar(p[filt][col+'x'],p[filt][col+'y'],xerr=np.array([p[filt][col+'exd'],p[filt][col+'exu']]),yerr=np.array([p[filt][col+'eyd'],p[filt][col+'eyu']]),c='black',ls='None',ms=ms,marker=mark,lw=lw,zorder=1,label=None)
    colordata = ax.scatter(p[filt][col+'x'],p[filt][col+'y'],c=colors,marker=mark,zorder=3,label=None)
    #Titles, axes, legends
    if c==0:
        ax.set_ylabel('log(O[III] 5007 / H$\\beta$)',fontsize = axisfont)
    ax.set_xlabel('log(N[II] 6583 / H$\\alpha$)',fontsize = axisfont)
    ax.tick_params(labelsize = ticksize, size=ticks)
    ax.set_xlim(-2,1)
    ax.set_ylim(-1.2,1.5)
    #Plot the bpt line
    ax.plot(xline,yline,color='black',lw=1,ls='--',label='Theoretical, Kewley+ (2001)')
    ax.plot(xlineemp,ylineemp,color='black',lw=1,ls='-',label='Empirical, Kauffmann+ (2003)')
    ax.text(0.06,0.18,'Star-Forming',fontsize = textfont+6, transform=ax.transAxes)
    ax.text(0.7,0.66,'AGN',fontsize = textfont+6, transform=ax.transAxes) 
    if c==0: ax.legend(fontsize=legendfont)
    c = c+1
cb = fig2.colorbar(colordata)
cb.ax.tick_params(labelsize = ticksize)
fig2.text(0.95,0.76,'log(Stellar Mass) ($M_{sun}$)',fontsize = axisfont,rotation=90)      
fig2.tight_layout()
fig2.savefig(figout + 'bpt_mass.pdf')
plt.close(fig2)
    

