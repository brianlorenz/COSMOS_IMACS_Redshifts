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
from scipy.optimize import curve_fit
#import lnr

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

#Read in the sfr file
sfdata = '/Users/blorenz/COSMOS/COSMOSData/sfrs.txt'
sfr_df = ascii.read(sfdata).to_pandas()
fluxdata = pd.merge(fluxdata,sfr_df,on='fluxfile')

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


HbHg = 0
HaHg = 0

if HbHg: lines = ['4861']
elif HaHg: lines = ['6563_fix']
else: lines=['6563_fix']

'''
#Set low objects to an upper limit
for line in lines:
    for i in range(0,len(fluxdata)):
        if (fluxdata.iloc[i][line+'_flux'] == 0) and (dataqual[line+'_low']=='True'):
            print 'Fixing'
            fluxdata.at[line+'_flux',i] = err_df.iloc[i][line+'_err']
'''
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



dupids = [item for item, count in collections.Counter(fluxdata[allgood]['OBJID']).items() if count > 1]
dupobjsarr = []
weight_df = pd.DataFrame()
weight_df['fluxfile'] = fluxdata.fluxfile
for i in range(len(fluxdata)):
    weight_df.at[i,'Weight'] = 1.0
for obj in dupids:
    dupobjsarr.append(fluxdata[fluxdata.OBJID == obj].OBJID)
for i in range(0,len(dupids)):
    ndup = len(dupobjsarr[i])
    for j in range(0,ndup):
        weight_df.at[dupobjsarr[i].index[j],'Weight'] = 1.0/ndup
fluxdata = pd.merge(fluxdata,weight_df,on='fluxfile')



combinemass = 1

paper = 0
showkel = 1
shows = 0
cutofflow = 1
showmeds = 1


filtSFR = fluxdata['SFR']<10000000

ms=12
lwbw=2

notbad = np.logical_not(baddata)

ssfr = 1

if not paper:
    fig = plt.figure(figsize = (19.5,7))
    ax = fig.add_axes((0.15,0.15,0.315,0.8))
else: fig,ax = plt.subplots(figsize = (8,7))
c=0
msbw = 12
lwbw = 3
colormed2 = 'black'
        

    

key = ''
if HbHg: key = '_HbHg'
elif HaHg: key = '_HaHg'
                    
for w in range(0,2):
    #llim1 = (np.log10(fluxdata['sSFR'+key]) > -10.7)
    #if cutofflow: ax.plot((0,20),(-10.7,-10.7),ls='--',color='black',label='sSFR cutoff for completeness')
    if cutofflow: ax.axhspan(-10.7,-15, color='indianred', alpha=0.1,label='Incomplete, discarding for analysis')
    if c in [0,3]:
        col = 'good'
        filt = allgood
        color='cornflowerblue'
        mark2 = 'o'
        label2 = 'Significant H$\\alpha$ detection'
    elif c in [1,4]:
        col = 'low'
        filt = somelow
        color='cornflowerblue'
        mark2 = 'v'
        label2 = '5$\\sigma$ Upper limit on SSFR'
    else:
        col = 'bad'
        filt = baddata
        color='red'
    filt = np.logical_and(filt,filtSFR)
    
    xdata = fluxdata[filt]['LMASS']
    ydata = np.log10(fluxdata[filt]['SFR'+key])
    ax.set_xlabel('log(Stellar Mass) (M$_\odot)$',fontsize = axisfont)
    ax.set_ylabel('log(SFR) (M$_{sun})$/yr)',fontsize = axisfont)
    if ssfr:
        ydata = np.log10(fluxdata[filt]['sSFR'+key])
        #Upper error
        yerru = np.log10(fluxdata[filt]['sSFR']+fluxdata[filt]['ssfr_err_u'])-np.log10(fluxdata[filt]['sSFR'])
        #If lower error is 0 or negative, set it to be very large
        ldiff = fluxdata[filt]['sSFR']-fluxdata[filt]['ssfr_err_d']
        ldiff.loc[ldiff<=0] = 10000
        yerrd = np.abs(np.log10(fluxdata[filt]['sSFR'])-np.log10(ldiff))
        ax.set_ylabel('log(sSFR) (yr$^{-1}$)',fontsize = axisfont)
        kelmodelcolor = 'orange'
        kelmodelw = 4
        kelz = 100
        if (c==0 and showkel): pkel = ax.plot((-100,100),(-9.46,-9.46),color=kelmodelcolor,ls='-',label='Model (Kelson 2014)',zorder=kelz,lw=kelmodelw)
        if (c==0 and showkel): ax.plot((-100,100),(-9.86,-9.86),color=kelmodelcolor,ls='--',label=None,zorder=kelz,lw=kelmodelw)
        if (c==0 and showkel): ax.plot((-100,100),(-9.06,-9.06),color=kelmodelcolor,ls='--',label=None,zorder=kelz,lw=kelmodelw)
        smodelcolor = 'orange'
        smodelw = 4
        sz = 100
        x2 = np.arange(9,9.4,0.01)
        x2b = np.arange(9.4,10,0.01)
        y2 = -0.17*(x2-10)-9.65
        y2b = -0.53*(x2b-10)-9.87
        if (c==0 and shows): psal = ax.plot(x2,y2,color=smodelcolor,ls='-',label='Fit to SDSS z<0.1 Galaxies (Salim+ 2007)',zorder=sz,lw=smodelw)
        if (c==0 and shows): ax.plot(x2b,y2b,color=smodelcolor,ls='-',label=None,zorder=sz,lw=smodelw)
        
        fluxdata['lsSFR'] = np.log10(fluxdata['sSFR'+key])
        mr1 = (fluxdata[notbad]['LMASS']<9.25)
        mr2 = np.logical_and(fluxdata[notbad]['LMASS']>=9.25,fluxdata[notbad]['LMASS']<9.5)
        mr3 = np.logical_and(fluxdata[notbad]['LMASS']>=9.5,fluxdata[notbad]['LMASS']<9.75)
        mr4 = (fluxdata[notbad]['LMASS']>=9.75)
        mrs = [mr1,mr2,mr3,mr4]
        if cutofflow:
                llim = (np.log10(fluxdata[notbad]['sSFR']) > -10.7)
                mrs = [np.logical_and(i,llim) for i in mrs]

        def getWmed(fluxdata, mr):
            sflux = fluxdata[notbad][mr].sort_values('lsSFR')
            cumsum = sflux.Weight.cumsum()
            cutoff = sflux.Weight.sum()/2.0
            median = sflux.lsSFR[cumsum>=cutoff].iloc[0]
            return median

        def geteWmed(fluxdata, mr):
            sflux = fluxdata.sort_values('sSFR')
            cumsum = sflux.Weight.cumsum()
            cutoff = sflux.Weight.sum()/2.0
            median = sflux.sSFR[cumsum>=cutoff].iloc[0]
            fluxdata['absSFR']=np.abs(fluxdata['sSFR']-median)
            sflux = fluxdata[notbad][mr].sort_values('absSFR')
            cumsum = sflux.Weight.cumsum()
            cutoff = sflux.Weight.sum()/2.0
            median = sflux.absSFR[cumsum>=cutoff].iloc[0]
            return median
        
       
        
        meds = np.array([getWmed(fluxdata,i) for i in mrs])
        emeds = 1.49*np.array([geteWmed(fluxdata,i) for i in mrs])
        emeds = (emeds/np.median(10**fluxdata['lsSFR'][notbad]))/2.303
        bins = np.array([9.125,9.375,9.625,9.875])
        #bins=(np.arange(1,17,2)/16.0)+9
        msbw = 12
        lwbw = 3
        colormed2 = 'black'
        if c==0:
                if showmeds: pmed = ax.errorbar(bins,meds,yerr=emeds,marker='o',ms=msbw,lw=lwbw,ls='None',zorder=1000,markerfacecolor='None', markeredgecolor=colormed2,mew=3,ecolor=colormed2,label='Median in bin, log(sSFR)>-10.7')
                def linefit(x,m,b):
                        y=m*x+b
                        return y
                coeff,pcov = curve_fit(linefit,bins,meds,sigma=np.array(emeds)/100)
                print coeff
                perr = np.sqrt(np.diag(pcov))
                sbins = np.array([8.5,10.5])
                if showmeds: ax.plot(sbins,linefit(sbins,coeff[0],coeff[1]),color='red',lw=4,ls='-',label='Fit to median',zorder=4)




        

    if c==0:
        pcir = ax.errorbar(xdata,ydata,yerr=np.array([yerrd,yerru]),color=color,marker=mark2,ms=4,lw=0.5,ls='None',zorder=10,label='Significant H$\\alpha$ detection')
        errfilt = yerru<0.1
        pcirdark = ax.errorbar(xdata[errfilt],ydata[errfilt],yerr=np.array([yerrd[errfilt],yerru[errfilt]]),color='indigo',marker=mark2,ms=4,lw=0.5,ls='None',zorder=11,label='Significant H$\\alpha$ detection (error <0.1 dex)')
    else:
        ptri = ax.plot(xdata,ydata,color=color,marker=mark2,mfc='None',ms=6,lw=0.5,ls='None',zorder=10,label=label2)
        #if HbHg: a = ax.plot((0,0),(0,0),color=color,marker='o',ms=4,lw=0.5,ls='None',zorder=1,label='Significant H$\\beta$ detection')
        #else: a = ax.plot((0,0),(0,0),color=color,marker='o',ms=4,lw=0.5,ls='None',zorder=1,label='Significant H$\\alpha$ detection')
        #b = ax.plot((0,0),(0,0),color=color,marker=mark2,mfc='None',ms=6,lw=0.5,ls='None',zorder=2,label=label2)
        #if showmeds: c1 = ax.errorbar(0,0,yerr=0.4,marker='o',ms=msbw,lw=lwbw,ls='None',zorder=3,markerfacecolor='None', markeredgecolor=colormed2,mew=3,ecolor=colormed2,label='Median in bin, log(sSFR)>-10.7')
        #if showkel: d = ax.plot((-100,0),(-9.46,-9.46),color=kelmodelcolor,ls='-',label='Model (Kelson 2014)',zorder=4,lw=kelmodelw)
        #if shows: e = ax.plot((0,0),(1,1),color=smodelcolor,ls='-',label='Empirical Fit (Salim 2007)',zorder=5,lw=smodelw)
        #if showmeds: f = ax.plot((0,0),(1,1),color='red',lw=4,ls='-',label='Fit to median',zorder=6)
        handles, labels = ax.get_legend_handles_labels()
        if not paper:
            '''
            hand = [a[0],b[0]]
            if showmeds: hand.append(c1[0])
            if showmeds: hand.append(f[0])
            if showkel: hand.append(d[0])
            if shows: hand.append(e[0])
            '''
            if (showmeds) and (showkel or shows): hand = [handles[-1],handles[-2],handles[2],handles[3],handles[5],handles[1],handles[0]]
            elif (showmeds) and (cutofflow): hand = [handles[-1],handles[-2],handles[1],handles[3],handles[4],handles[0]]
            elif (cutofflow): hand = [handles[-1],handles[-2],handles[0],handles[1]]
            else: hand = [handles[-1],handles[-2],handles[0]]
            ax.legend(handles=hand,fontsize=axisfont-2,bbox_to_anchor=(1.01, 0.5))
        else:
            hand = [handles[-1],handles[-2],handles[2],handles[3],handles[5],handles[1],handles[0]]
            ax.legend(handles=hand,fontsize=legendfont-6,loc=1,frameon=False)
    ax.tick_params(labelsize = ticksize, size=ticks)
    ax.set_xlim(8.95,10.05)
    ax.set_ylim(-3,1.6)
    
    if ssfr:
        if not paper:
            if HbHg: ax.set_ylim(-13,-6)
            else: ax.set_ylim(-12,-8)

        else:
            ax.set_ylim(-12,-7)
    c=c+1

fig.tight_layout()
if ssfr:
    if HbHg: fig.savefig(figout + 'sSFR_Mass_HbHg.pdf')
    elif HaHg: fig.savefig(figout + 'sSFR_Mass_HaHg.pdf')
    else: fig.savefig(figout + 'sSFR_Mass.pdf')
    
else: fig.savefig(figout + 'SFR_Mass.pdf')
plt.close(fig)
