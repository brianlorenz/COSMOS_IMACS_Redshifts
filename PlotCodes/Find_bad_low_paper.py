#Find which objects are bad and low based on various cuts through the data. Output this as a dataframe containing True False for every object in every line
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import sys, os, string
import pandas as pd
from astropy.io import fits
import collections
from astropy.stats import biweight_midvariance
from matplotlib.font_manager import FontProperties

#Folder to save the figures
figout = '/Users/blorenz/COSMOS/Reports/2018/Images/'

#The location with the file for all of our data
qualout = '/Users/blorenz/COSMOS/COSMOSData/dataqual.txt'

#The location with the file for all of our data
fluxdatapath = '/Users/blorenz/COSMOS/COSMOSData/lineflux.txt'


#The location to store the scale and its stddev of each line
scaledata = '/Users/blorenz/COSMOS/COSMOSData/scales.txt'
#Read in the scale of the lines 
#scale_df = ascii.read(scaledata).to_pandas()

#Location of the equivalent width data
ewdata = '/Users/blorenz/COSMOS/COSMOSData/lineew.txt'
#Read in the ew of the lines 
ew_df = ascii.read(ewdata).to_pandas()

#Read the datafile:
fluxdata = ascii.read(fluxdatapath).to_pandas()

#Division function
def divz(X,Y):
        return X/np.where(Y,Y,Y+1)*np.not_equal(Y,0)
    

sig = 5

lines = ['3727','4102','4340','4861','4959','5007','6548','6548_fix','6563','6563_fix','6583','6583_fix']
lines = np.sort(lines)
fig,axarr = plt.subplots(4,12,figsize=(100,30))


#Plotting parameters
mark='o'
ms=6
ls='--'
#Fontsizes for plotting
axisfont = 18
ticksize = 16
titlefont = 32
legendfont = 16
textfont = 16

med_bi_df = pd.DataFrame()
qualframe = pd.DataFrame()

cullnames = ['_scale','_rchi2','_shift','_stddev']

def readdata(cull,line):
    #Special case for shift - we need to compute this, and not fit in log
    if cull in ['_shift']:
        if len(line)==8: lineval=int(line[0:4])
        else: lineval = int(line)
        #Compute the shift
        culldata = (fluxdata[line+'_mean']-(1+fluxdata['zcc'])*lineval)/2
        culldatapos=culldata
    elif cull in ['_stddev']:
        culldata = fluxdata[line+cull]
        culldatapos = culldata
    else:
        #Set the data to cull
        culldata = np.log10(fluxdata[line+cull])
        culldatapos = np.log10(fluxdata[fluxdata[line+cull]>0][line+cull])
    return culldata,culldatapos
    

#Loop over all lines
for line in lines:
    #Read in the data that won't change
    usig = fluxdata[line+'_usig']
    err = np.sqrt(fluxdata[line+'_usig']**2+(0.25*ew_df[line+'_modelabs'])**2)
    rchi2 = fluxdata[line+'_rchi2']
    flux = fluxdata[line+'_flux']
    usigrat = divz(flux,(err*np.sqrt(rchi2)))
    #Find the low data at a cut of sig signal-to-noise
    low = (usigrat < sig)
    badflag = (fluxdata[line+'_flag'] != 0)
    qualframe[line+'_temp_low'] = low
    #Loop over all culling methods
    for cull in cullnames:
        culldata, culldatapos = readdata(cull,line)
        #We cull in log so that 0s are filtered out
        #Find the median and biweight of the distribution of objects
        med = np.median(culldatapos)
        biwt = np.sqrt(biweight_midvariance(culldatapos))
        #Cut the data 3 sigma from the median
        bad = (np.abs(culldata-med) > (3*biwt))
        bad = np.logical_or(bad,badflag)
        #Set the good, bad, and low for this cut
        low = np.logical_and(low,np.logical_not(bad))
        good = np.logical_and(np.logical_not(low),np.logical_not(bad))
        #Store the med, biwt in the dataframe
        med_bi_df.at[0,line+cull] = med
        med_bi_df.at[1,line+cull] = biwt
        #Store the good, low, and bad into the frame
        qualframe[line+cull+'_good'] = good
        qualframe[line+cull+'_low'] = low
        qualframe[line+cull+'_bad'] = bad


counter = 0
#Now that the frame is complete, combine all of the bads to flag the points
for line in lines:
    cullcount = 0
    #Read in the fluxes
    flux = fluxdata[line+'_flux']
    #Get the bad data across all cuts
    bads = [qualframe[line+name+'_bad'] for name in cullnames]
    #If it is bad in any of the cuts, it is bad data
    badline = np.logical_or.reduce(bads)
    #Find the low data that is not bad
    lowline = np.logical_and(qualframe[line+'_temp_low'],np.logical_not(badline))
    goodline = np.logical_and(np.logical_not(lowline),np.logical_not(badline))
    #Store the low, bad and good for each line across all methods
    qualframe[line+'_bad'] = badline
    qualframe[line+'_low'] = lowline
    qualframe[line+'_good'] = goodline
    #Now, plot the culling for each line:
    for cull in cullnames:
        #Read the data
        culldata, culldatapos = readdata(cull,line)
        #Set the current plot
        ax = axarr[cullcount,counter]
        med = med_bi_df.iloc[0][line+cull]
        biwt = med_bi_df.iloc[1][line+cull] 
        #Plot the cut
        ax.plot(flux[goodline],culldata[goodline],color='blue',marker=mark,ms=ms,label='Detection',ls='None')
        ax.plot(flux[lowline],culldata[lowline],color='orange',marker=mark,ms=ms,label='Low S/N',ls='None')
        ax.plot(flux[badline],culldata[badline],color='red',marker=mark,ms=ms,label='Not usuable',ls='None')
        ax.plot((-100,10000),(med,med),color='black',ls=ls,label='Median')
        #ax.plot((-100,10000),(med+biwt,med+biwt),color='darkgrey',ls=ls,label='1 sigma')
        #ax.plot((-100,10000),(med-biwt,med-biwt),color='darkgrey',ls=ls)
        ax.plot((-100,10000),(med+3*biwt,med+3*biwt),color='grey',ls=ls,label='3 sigma')
        ax.plot((-100,10000),(med-3*biwt,med-3*biwt),color='grey',ls=ls)
        ax.legend(fontsize=legendfont)
        font = FontProperties()
        font.set_weight('bold')
        ax.set_xlabel('H$\\alpha$ Flux ($10^{-17}$ erg/s/${cm}^2/\AA$)',fontsize=axisfont)
        if cull == '_scale':
            ax.set_title(line,fontsize=titlefont)
            ax.set_ylabel('log(Scale)',fontsize=axisfont)
            ax.set_ylim(-2,2)
        if cull == '_rchi2':
            ax.set_ylabel('log($\\chi_\\nu^2$)',fontsize=axisfont)
            ax.set_ylim(-3,3)
        if cull == '_shift':
            ax.set_ylim(-8,8)
            ax.set_ylabel('Shift' +' (pixels)',fontsize=axisfont)
        if cull == '_stddev':
            ax.set_ylim(0,10)
            ax.set_ylabel('Sigma'+' ($\AA$)',fontsize=axisfont)
        ax.tick_params(labelsize = ticksize)
        ax.set_xscale('log')
        ax.set_xlim(0.001,500)
        #ax.set_yscale('log')
        cullcount = cullcount+1
    counter = counter+1
    

#qualframe.to_csv(qualout,index=False)

fig.tight_layout()
fig.savefig(figout + 'data_culling.pdf')
plt.close(fig)
