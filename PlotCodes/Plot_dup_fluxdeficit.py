import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import sys, os, string
import pandas as pd
from astropy.io import fits
import collections
from astropy.stats import biweight_midvariance
import matplotlib.patches as patches

#Folder to save the figures
figout = '/Users/blorenz/COSMOS/Reports/2018/Images/'

#The location with the file for all of our data
fluxdatapath = '/Users/blorenz/COSMOS/COSMOSData/lineflux.txt'

#File with the error array
errdatapath = '/Users/blorenz/COSMOS/COSMOSData/errs.txt'
#Read in the scale of the lines 
err_df = ascii.read(errdatapath,data_start=1,header_start=0,format='csv').to_pandas()

#The location to store the scale and its stddev of each line
qualdatapath = '/Users/blorenz/COSMOS/COSMOSData/dataqual.txt'
#Read in the scale of the lines 
dataqual = ascii.read(qualdatapath).to_pandas()
d = {'True': True, 'False': False}

#Merge our data with the UVISTA catalog
#The location of the muzzin et al data:
#Location of the equivalent width data
ewdata = '/Users/blorenz/COSMOS/COSMOSData/lineew.txt'

#Location to save scatter in duplicates
biwtout = '/Users/blorenz/COSMOS/COSMOSData/biwt_dup.txt'

#The location to store the scale and its stddev of each line
scaledata = '/Users/blorenz/COSMOS/COSMOSData/scales.txt'
#Read in the scale of the lines 
scale_df = ascii.read(scaledata).to_pandas()

#Read in the ew of the lines 
ew_df = ascii.read(ewdata).to_pandas()

#Amount by which to multiply absorption error
errabs = 0.25



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

fluxdata = pd.merge(fluxdata,ew_df,on='fluxfile')

#All of the lines to plot
lines = ['4861','5007','4959','6563_fix','6583_fix','6548_fix','3727','4340','4102','6563','6583','6548']
lines = np.sort(lines)

#Make the figure




lw=0.25
mark='o'
minlim = 0.01
maxlim = 500

#Finds the duplicate objects and computes the ew for them
def finddups(line):
    #Filter the data
    allgood = dataqual[line+'_good'].map(d)
    #Needs to be bad in any line to be bad
    baddata = dataqual[line+'_bad'].map(d)
    somelow = dataqual[line+'_low'].map(d)
    goodflux = fluxdata[np.logical_not(baddata)]
    #Get the duplicate OBJIDs of only the good objects
    dupids = [item for item, count in collections.Counter(goodflux['OBJID']).items() if count > 1]
    dupobjsarr = []
    lineflux1 = []
    lineflux2 = []
    lineflux3 = []
    lineflux4 = []
    lineflux5 = []
    lineflux6 = []
    #Append the data for the good objects
    for obj in dupids:
        dupobjsarr.append(fluxdata[fluxdata.OBJID == obj])
    #Get the ews for these objects
    for i in range(0,len(dupobjsarr)):
        lineflux1.append(dupobjsarr[i].iloc[0][line+'_flux'])
        lineflux2.append(dupobjsarr[i].iloc[1][line+'_flux'])
        lineflux3.append(dupobjsarr[i].iloc[0][line+'_flux'])
        lineflux4.append(dupobjsarr[i].iloc[1][line+'_flux'])
        lineflux5.append(dupobjsarr[i].iloc[0][line+'_modelabs'])
        lineflux6.append(dupobjsarr[i].iloc[1][line+'_modelabs'])
    lineflux3 = np.array(lineflux3)
    lineflux4 = np.array(lineflux4)
    lineflux5 = np.array(lineflux5)
    lineflux6 = np.array(lineflux6)
    #Make another cut for high S/N
    flagidx = allgood
    #Repeat the process just for the high SN objects
    goodflux = fluxdata[flagidx]
    linefluxsig1 = []
    linefluxsig2 = []
    dupids = [item for item, count in collections.Counter(goodflux['OBJID']).items() if count > 1]
    dupobjsarr = []
    for obj in dupids:
        dupobjsarr.append(fluxdata[fluxdata.OBJID == obj])
    for i in range(0,len(dupobjsarr)):
        linefluxsig1.append(dupobjsarr[i].iloc[0][line+'_flux'])
        linefluxsig2.append(dupobjsarr[i].iloc[1][line+'_flux'])
    fluxdiffmed = np.median(np.abs(lineflux4-lineflux3))
    fluxdiffbi = np.sqrt(biweight_midvariance(np.abs(lineflux4-lineflux3)))
    absdiffmed = np.median(np.abs(lineflux6-lineflux5))
    absdiffbi = np.sqrt(biweight_midvariance(np.abs(lineflux6-lineflux5)))
    diffarr = [fluxdiffmed, fluxdiffbi, absdiffmed, absdiffbi]
    return lineflux1,lineflux2, linefluxsig1, linefluxsig2, diffarr


fig,axarr = plt.subplots(3,4,figsize=(30,20))
axarr = np.reshape(axarr,12)
counter = 0
biwt_df = pd.DataFrame()
biwt_df['fluxfile'] = fluxdata['fluxfile']
#Loop over the lines
for line in lines:
    ax = axarr[counter]
    #Get the ew of the duplicates 
    xdata,ydata,xdatasn,ydatasn,diffarr = finddups(line)
    #Plot the data
    ax.scatter(xdata,ydata,color='orange',marker=mark,label='Below 5 sigma')
    ax.scatter(xdatasn,ydatasn,color='blue',marker=mark, label='Above 5 sigma')
    #Plot a unity line
    ax.plot((0,1000),(0,1000),color='black',ls='--')
    #Titles, axes, legends
    #ax.legend(fontsize=axisfont)
    ax.set_xlabel(line + ' Flux ($10^{-17}$ erg/s/$cm^2$)',fontsize = axisfont)
    ax.set_ylabel(line + ' Flux ($10^{-17}$ erg/s/$cm^2$)',fontsize = axisfont)
    biwt = np.sqrt(biweight_midvariance(np.abs(np.array(xdata)-np.array(ydata))))
    for i in range(len(biwt_df)):
        biwt_df.at[i,line] = biwt
    #Add a line to show the mad of the scatter in ew
    ax.text(0.02,0.96,'sqrt(biweight): ' + str(round(biwt,2)),fontsize = textfont, transform=ax.transAxes)      
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(minlim,maxlim)
    ax.set_ylim(minlim,maxlim)
    ax.tick_params(labelsize = ticksize)
    counter = counter+1
fig.tight_layout()
fig.savefig(figout + 'dup_flux.pdf')
plt.close(fig)
biwt_df.to_csv(biwtout,index=False)
'''

fig2,axarr2 = plt.subplots(3,4,figsize=(30,20))
axarr2 = np.reshape(axarr2,12)
counter=0
for line in lines:
    ax = axarr2[counter]
    #Get the ew of the duplicates 
    good = fluxdata[np.logical_or((fluxdata[line + '_flag'] == 0),(fluxdata[line + '_flag'] == 4))]
    good = good[np.abs(np.log10(fluxdata[line+'_scale'])-scale_df[line+'_medscale'][0]) < 3*scale_df[line+'_sigscale'][0]]
    #Plot the data
    #ax.scatter(low[line+'_flux'],low[line+'_EW'],color='blue',marker=mark,label='Below 3*1.49*mad')
    ax.scatter(good[line+'_flux'],good[line+'_EW'],color='blue',marker=mark,label='Good data')
    mad = mad_df[line+'_mad'][0]
    ax.plot((mad,mad),(-1000,1000),color='black',ls='--',label='1.49*MAD')
    #Titles, axes, legends
    ax.legend(fontsize=axisfont)
    ax.set_title(line + ' Flux vs EW, ' + str(len(good[line+'_flux'])) + ' Objects',fontsize = titlefont)
    ax.set_xlabel(line + ' Flux ($10^{-17}$ erg/s/$cm^2$)',fontsize = axisfont)
    ax.set_ylabel(line + ' EW ($\AA$)',fontsize = axisfont)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(minlim,maxlim)
    ax.set_ylim(minlim,maxlim)
    ax.tick_params(labelsize = ticksize)
    counter = counter+1
fig2.tight_layout()
fig2.savefig(figout + 'ew_flux.pdf')
plt.close(fig2)

fig3,axarr3 = plt.subplots(3,4,figsize=(30,20))
axarr3 = np.reshape(axarr3,12)
counter=0
for line in lines:
    minlim,maxlim = 0.001,500
    ax = axarr3[counter]
    #Plot the data
    ax.scatter(fluxdata[line+'_modelEW'],fluxdata[line+'_EW'],color='blue',marker=mark)
    modelbi = np.sqrt(biweight_midvariance(fluxdata[line+'_modelEW']))
    modelmed = np.median(fluxdata[line+'_modelEW'])
    databi = np.sqrt(biweight_midvariance(fluxdata[line+'_EW']))
    datamed = np.median(fluxdata[line+'_EW'])
    ax.errorbar(modelmed,datamed,xerr=modelbi,color='red',marker='o',ms=14,lw=4,ls='None')
    ax.text(0.02,0.02,'sqrt(model biweight): ' + str(np.round(modelbi,3)),fontsize = axisfont-2, transform=ax.transAxes)
    #Titles, axes, legends
    ax.set_title(line + ' EW vs model EW',fontsize = titlefont)
    ax.set_ylabel(line + ' EW ($\AA$)',fontsize = axisfont)
    ax.set_xlabel(line + ' Model EW ($\AA$)',fontsize = axisfont)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(minlim,maxlim)
    ax.set_ylim(minlim,maxlim)
    ax.tick_params(labelsize = ticksize)
    counter = counter+1
fig3.tight_layout()
fig3.savefig(figout + 'ew_modelew.pdf')
plt.close(fig3)

fig4,axarr4 = plt.subplots(3,4,figsize=(30,20))
axarr4 = np.reshape(axarr4,12)
counter=0
for line in lines:
    minlim,maxlim = 0.001,10
    ax = axarr4[counter]
    good = fluxdata[fluxdata[line+'_modelEW'] > 0]
    good = good[good[line+'_modelabs'] > 0]
    #Plot the data
    ax.scatter(good[line+'_modelabs'],good[line+'_modelEW'],color='blue',marker=mark)
    modelbi = np.sqrt(biweight_midvariance(good[line+'_modelEW']))
    modelmed = np.median(good[line+'_modelEW'])
    absbi = np.sqrt(biweight_midvariance(good[line+'_modelabs']))
    absmed = np.median(good[line+'_modelabs'])
    ax.errorbar(absmed,modelmed,xerr=absbi,yerr=modelbi,color='red',marker='o',ms=14,lw=4,ls='None')
    ax.text(0.02,0.02,'sqrt(model biweight): ' + str(np.round(modelbi,3)),fontsize = axisfont-2, transform=ax.transAxes)
    ax.text(0.02,0.08,'sqrt(absorb biweight): ' + str(np.round(absbi,3)),fontsize = axisfont-2, transform=ax.transAxes)
    #Plot a unity line
    ax.plot((0,1000),(0,1000),color='black',ls='--')
    #Titles, axes, legends
    ax.set_title(line + ' Model Absorption vs EW',fontsize = titlefont)
    ax.set_xlabel(line + ' Model Absorption Flux ($10^{-17}$ erg/s/$cm^2$)',fontsize = axisfont)
    ax.set_ylabel(line + ' Model EW ($\AA$)',fontsize = axisfont)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(minlim,maxlim)
    ax.set_ylim(minlim,maxlim)
    ax.tick_params(labelsize = ticksize)
    counter = counter+1
fig4.tight_layout()
fig4.savefig(figout + 'model_abs_flux.pdf')
plt.close(fig4)
'''
def geterr(pd_df,line):
    err = np.sqrt(pd_df[line+'_usig']**2+(errabs*ew_df.iloc[pd_df.index.values][line+'_modelabs'])**2)
    return err

fig5,axarr5 = plt.subplots(3,4,figsize=(30,20))
axarr5 = np.reshape(axarr5,12)
counter=0

for line in lines:
    minlim,maxlim = 0.001,100
    ax = axarr5[counter]
    #Filter the data
    allgood = dataqual[line+'_good'].map(d)
    #Needs to be bad in any line to be bad
    baddata = dataqual[line+'_bad'].map(d)
    somelow = dataqual[line+'_low'].map(d)
    good = fluxdata[np.logical_not(baddata)]
    highsn = fluxdata[allgood]
    xdata,ydata,xdatasn,ydatasn,diffarr = finddups(line)
    #Get the errors
    gooderry = geterr(good,line)
    gooderrhsn = geterr(highsn,line)
    #Plot the data
    ax.errorbar(good[line+'_modelabs'],good[line+'_flux'],xerr=errabs*good[line+'_modelabs'],yerr=gooderry,color='blue',marker=mark,label='<10 sigma',ls='None',zorder=1)
    ax.errorbar(highsn[line+'_modelabs'],highsn[line+'_flux'],xerr=errabs*highsn[line+'_modelabs'],yerr=gooderrhsn,color='darkblue',marker=mark,label='>10 sigma',ls='None',zorder=2)
    fluxbi = diffarr[1]
    fluxmed = diffarr[0]
    absbi = diffarr[3]
    absmed = diffarr[2]
    #ax.errorbar(absmed,fluxmed,xerr=absbi,yerr=fluxbi,color='red',marker='o',ms=14,lw=4,ls='None')
    ax.plot(absmed,fluxmed,color='red',marker='o',ms=14,lw=4,ls='None',label='Dup center',zorder=4)
    rect = patches.Rectangle((absmed-absbi,fluxmed-fluxbi),2*absmed,2*fluxbi,linewidth=1,edgecolor='r',facecolor='pink',alpha=0.5,label='Dup sqrt(biweight)',zorder=3)
    ax.add_patch(rect)
    ax.text(0.54,0.02,'sqrt(y biweight): ' + str(np.round(fluxbi,3)),fontsize = axisfont-2, transform=ax.transAxes)
    ax.text(0.54,0.08,'sqrt(x biweight): ' + str(np.round(absbi,3)),fontsize = axisfont-2, transform=ax.transAxes)
    #Plot a unity line
    ax.plot((0,1000),(0,1000),color='black',ls='--',zorder=5)
    ax.plot((0,1000),(0,3000),color='grey',ls='--',zorder=6)
    #Titles, axes, legends
    ax.set_title('Model Absorption vs ' + line + ' Flux',fontsize = titlefont)
    ax.set_xlabel(line + ' Model Absorption Flux ($10^{-17}$ erg/s/$cm^2$)',fontsize = axisfont)
    ax.set_ylabel(line + ' Emission Flux ($10^{-17}$ erg/s/$cm^2$)',fontsize = axisfont)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(minlim,maxlim)
    ax.set_ylim(minlim,maxlim)
    ax.tick_params(labelsize = ticksize)
    counter = counter+1
fig5.tight_layout()
fig5.savefig(figout + 'model_flux_deficit.pdf')
plt.close(fig5)

