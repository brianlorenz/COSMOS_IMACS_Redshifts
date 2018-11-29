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
fluxout = '/Users/blorenz/COSMOS/COSMOSData/lineflux_zero.txt'

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

#File with the error array
errdatapath = '/Users/blorenz/COSMOS/COSMOSData/errs.txt'
#Read in the scale of the lines 
err_df = ascii.read(errdatapath,data_start=1,header_start=0,format='csv').to_pandas()


#Read the datafile:
fluxdata = ascii.read(fluxdatapath).to_pandas()

#Division function
def divz(X,Y):
        return X/np.where(Y,Y,Y+1)*np.not_equal(Y,0)
    

sig = 5

skylines = np.array([5577.34,5889.95,5895.92])

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


#Loop over all lines to find bad objects
for line in lines:
    #Read in the data that won't change
    usig = fluxdata[line+'_usig']
    err = err_df[line]
    rchi2 = fluxdata[line+'_rchi2']
    flux = fluxdata[line+'_flux']
    usigrat = divz(flux,(err*np.sqrt(rchi2)))
    #Find the low data at a cut of sig signal-to-noise
    low = (usigrat < sig)
    badflag = (fluxdata[line+'_flag'] != 0)

    #Find where rchi2 is too small (too big is ok since they may not be gaussians)
    lrchi2 = np.log10(fluxdata[line+'_rchi2'])
    lrchi2gt0 = np.log10(fluxdata[fluxdata[line+'_rchi2']>0][line+'_rchi2'])
    lrchi2med = np.median(lrchi2gt0)
    lrchi2biwt = np.sqrt(biweight_midvariance(lrchi2gt0))
    badrchi2 = (lrchi2<(lrchi2med-5*lrchi2biwt))

    #Find where the scale goes bad
    scale = fluxdata[line+'_scale']
    scalemed = np.median(scale)
    scalebwt = np.sqrt(biweight_midvariance(scale))
    scalebad = (scale>1.5*scalemed)+(scale<0.5*scalemed)+(scale<=0)

    #Find where the line hits skylines
    if len(line)==8: lineval=int(line[0:4])
    else: lineval = int(line)
    zline = (lineval*(1+fluxdata['zcc']))
    skybad0 = np.abs(zline-skylines[0])
    skybad1 = np.abs(zline-skylines[1])
    skybad2 = np.abs(zline-skylines[2])
    skybad = (skybad0<16)+(skybad1<16)+(skybad2<16)

    #Find if the usig is too low
    usigbad=(usig<0)
    #Flag if err is too high
    errbad = (err>50)
        

    #Find objects that are bad
    bad = badrchi2+scalebad+skybad+badflag+usigbad+errbad

    #Flag these objects as bad
    qualframe[line+'_bad'] = bad

#Find mean sigma of each galaxy across all lines
allbads = np.asarray([qualframe[line+'_bad'] for line in lines])
allgoods = np.logical_not(allbads)
allsigs = np.asarray([fluxdata[line+'_stddev'] for line in lines])

refsigs = divz(np.add.reduce(allgoods*allsigs),np.add.reduce(allgoods))
for line in lines:
    qualframe[line+'_refstddev'] = refsigs
    qualframe[line+'_ratstddev'] = divz(fluxdata[line+'_stddev'],refsigs)
    
#Loop over all lines to find objects to set to zero
for line in lines:
    #Read in the data that won't change
    usig = fluxdata[line+'_usig']
    err = err_df[line]
    rchi2 = fluxdata[line+'_rchi2']
    flux = fluxdata[line+'_flux']
    usigrat = divz(flux,(err*np.sqrt(rchi2)))
    #Find the low data at a cut of sig signal-to-noise
    low = (usigrat < sig)
    goodrats = np.logical_and(np.logical_not(qualframe[line+'_bad']),(fluxdata[line+'_stddev']>2.001))
    userats = qualframe[goodrats][line+'_ratstddev']
    ratbwt = np.sqrt(biweight_midvariance(userats))
    ratmed = np.median(userats)
    lowsig = (qualframe[line+'_ratstddev']>(ratmed+2.5*ratbwt))

    #Find the shift
    if len(line)==8: lineval=int(line[0:4])
    else: lineval = int(line)
    shift = (fluxdata[line+'_mean']-(1+fluxdata['zcc'])*lineval)/2
    shiftmed = np.median(shift[np.logical_not(qualframe[line+'_bad']+lowsig)])
    shiftcut = 1.5
    lowshift = (shift<shiftmed-shiftcut)+(shift>shiftmed+shiftcut)
    
    setzero = np.logical_or(lowshift,lowsig)
    fluxdata.loc[setzero,line+'_flux'] = 0    

    low = np.logical_or(setzero,low)
    low = np.logical_and(low,np.logical_not(qualframe[line+'_bad']))
    good = np.logical_and(np.logical_not(low),np.logical_not(qualframe[line+'_bad']))
    
    #Store the good, low, and bad into the frame
    qualframe[line+'_good'] = good
    qualframe[line+'_low'] = low
    
cullnames = ['_scale','_rchi2','_shift','_stddev']

counter = 0
#Now that the frame is complete, combine all of the bads to flag the points
for line in lines:
    cullcount = 0
    #Read in the fluxes
    flux = fluxdata[line+'_flux']
    #Get the bad data 
    badline = qualframe[line+'_bad']
    #Find the low data that is not bad
    lowline = qualframe[line+'_low']
    goodline = qualframe[line+'_good']
    #Now, plot the culling for each line:
    for cull in cullnames:
        #Set the current plot
        ax = axarr[cullcount,counter]
        #med = med_bi_df.iloc[0][line+cull]
        #biwt = med_bi_df.iloc[1][line+cull]
        if cull == '_shift':
            if len(line)==8: lineval=int(line[0:4])
            else: lineval = int(line)
            culldata = (fluxdata[line+'_mean']-(1+fluxdata['zcc'])*lineval)/2
        elif cull == '_rchi2':
            culldata = np.log10(fluxdata[line+'_rchi2'])
        elif cull == '_scale':
            culldata = np.log10(fluxdata[line+'_scale'])
        else:
            culldata=fluxdata[line+cull]
        #If the data is zero, plot it at 0.001
        zerofilt = (fluxdata[line+'_flux']==0)
        goodz = np.logical_and(zerofilt,goodline)
        lowz = np.logical_and(zerofilt,lowline)
        badz = np.logical_and(zerofilt,badline)
        zeros = np.zeros(len(fluxdata))+0.002
        #Plot the cut
        ax.plot(flux[goodline],culldata[goodline],color='blue',marker=mark,ms=ms,label='Good Data',ls='None')
        ax.plot(flux[lowline],culldata[lowline],color='orange',marker=mark,ms=ms,label='Low Data',ls='None')
        ax.plot(flux[badline],culldata[badline],color='red',marker=mark,ms=ms,label='Bad Data',ls='None')
        ax.plot(zeros[goodz],culldata[goodz],color='blue',marker=mark,ms=ms,label='Good Data',ls='None')
        ax.plot(zeros[lowz],culldata[lowz],color='orange',marker=mark,ms=ms,label='Low Data',ls='None')
        ax.plot(zeros[badz],culldata[badz],color='red',marker=mark,ms=ms,label='Bad Data',ls='None')
        '''
        ax.plot((-100,10000),(med,med),color='black',ls=ls,label='Median')
        ax.plot((-100,10000),(med+biwt,med+biwt),color='darkgrey',ls=ls,label='1 sigma')
        ax.plot((-100,10000),(med-biwt,med-biwt),color='darkgrey',ls=ls)
        ax.plot((-100,10000),(med+3*biwt,med+3*biwt),color='grey',ls=ls,label='3 sigma')
        ax.plot((-100,10000),(med-3*biwt,med-3*biwt),color='grey',ls=ls)
        '''
        font = FontProperties()
        font.set_weight('bold')
        #ax.text(0.02,0.02,'Median: ' + str(np.round(med,3)),fontsize = axisfont, transform=ax.transAxes,fontproperties=font)
        #ax.text(0.57,0.02,'sqrt(biweight): ' + str(np.round(biwt,3)),fontsize = axisfont, transform=ax.transAxes,fontproperties=font)
        ax.set_xlabel(line+' Flux ($10^{-17}$ erg/s/${cm}^2/\AA$)',fontsize=axisfont)
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
    

qualframe.to_csv(qualout,index=False)
fluxdata.to_csv(fluxout,index=False)

fig.tight_layout()
fig.savefig(figout + 'data_culling.pdf')
plt.close(fig)

