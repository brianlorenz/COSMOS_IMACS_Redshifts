#Finds the ratio between our observed spectra and the ULTRAVista ones so that we can flux calibrate the data

#Usage:      run SED_Ratio.py a6      to find the ratio for the a6 mask. 

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import sys, os, string
import pandas as pd
from astropy.io import fits
from astropy.convolution import convolve, Box1DKernel
from scipy.interpolate import splrep, splev
from scipy.signal import medfilt

#Folder to save the figures
figout = '/Users/blorenz/COSMOS/COSMOSData/SED_Ratio_Images/'

#The location with the file for all of our data
ourdatapath = '/Users/blorenz/COSMOS/COSMOSData/all_c_hasinger.txt'
#Location for data files from the SED fitting code
seddatapath = '/Users/blorenz/COSMOS/COSMOSData/SEDmodels/'
#Location for all of the .fits files
fitsdatapath = '/Users/blorenz/COSMOS/COSMOSData/corFitsFileOut/'
#Where the calibrated spectra are stored (not necessary)
caldatapath = '/Users/blorenz/COSMOS/COSMOSData/flxFitsFileOut/'

letnum = sys.argv[1]

#Read in all of our data
ourdata = ascii.read(ourdatapath).to_pandas()
ourdata = ourdata[ourdata.ImageName.str.contains('feb1' + letnum[1] + '_' + letnum[0]) == True]
ourdata = ourdata[ourdata.Unsure == 0]
ourdata = ourdata[ourdata.Bad == 0]
ourdata = ourdata[ourdata.Flag3 == 0]
ourdata = ourdata[ourdata.Flag1 == 0]

#Division function
def divz(X,Y):
        return X/np.where(Y,Y,Y+1)*np.not_equal(Y,0)

#Polynomial fitting funciton
def svdfit(b,y):
        decomp = np.linalg.svd(b,full_matrices=False)
        sol1 = np.transpose(decomp[2])
        sol2 = divz(1.0,decomp[1])
        sol3 = np.dot(np.transpose(decomp[0]),y)
        if np.sometrue(sol3):
            solr = (sol2*sol3)
            soll = np.dot(sol1,solr)
        else:
            soll = np.zeros(sol3.shape)
        return soll

#Array to store the model ratios
modelrat = []
#Array to store the flats
flatarr = []

#Change to for loop to cycle through the data
for i in range(0,len(ourdata)):
    #Find the identifying pieces of the string. e.g. 'a' '6' '026110'
    let = ourdata.iloc[i].ImageName[17]
    num = ourdata.iloc[i].ImageName[15]
    objid = ourdata.iloc[i].ImageName[4:10]

    #Read the corresponding SED model, e.g. '026110_a6.dat'
    if os.path.exists(seddatapath + objid + '_' + let + num + '.dat'):
        sedmodel = ascii.read(seddatapath + objid + '_' + let + num + '.dat').to_pandas()
        #Read the corresponding spectrum 
        hdu = fits.open(fitsdatapath + ourdata.iloc[i].ImageName)[0]
        head = hdu.header
        fitsdata = hdu.data
        #Calculate the wavelength range for the data
        crval1 = head["crval1"]
        crpix1 = head["crpix1"]
        cdelt1 = head["cdelt1"]
        naxis1 = head["naxis1"]
        dcflag = head["dc-flag"]
        wavelength = (1.0+np.arange(naxis1)-crpix1)*cdelt1 + crval1
        if dcflag: wavelength = np.power(10.0,wavelength)

        #Split the fits data into its compenents
        spec = fitsdata[0]
        orig = fitsdata[1]
        noise = fitsdata[2]
        flat = fitsdata[3]
        sky = fitsdata[4]
        tc = fitsdata[5]
        mask = fitsdata[6]
        fit = fitsdata[7]

        #Ranges where we will check the red, green, blue signal to noise
        blueidx = np.logical_and(wavelength > 5500, wavelength < 6000)
        greenidx = np.logical_and(wavelength > 6500, wavelength < 7500)
        redidx = np.logical_and(wavelength > 8500, wavelength < 9000)

        #Check if the signal to noise ratio [S/N] is good enough
        sn = divz(spec,noise)
        bluesn = np.median(sn[blueidx])
        greensn = np.median(sn[greenidx])
        redsn = np.median(sn[redidx])
        if (redsn > 2 and greensn > 10 and bluesn > 5):
        
            #Clip the model to the size of the spectrum:
            idxarr = np.where(np.logical_and(sedmodel['lambda']>=np.min(wavelength), sedmodel['lambda']<=np.max(wavelength)))
            sedmodel = sedmodel.truncate(np.min(idxarr),np.max(idxarr))
            
            #Make the model the same resolution as the spectrum
            sdtck = splrep(sedmodel["lambda"],sedmodel["model"],s=0,k=3)
            himodel = splev(wavelength,sdtck)
            
            #Divide the spectrum by the model to get the ratio
            modelrat.append(divz(spec,himodel))
            #Add the flat to the array only if everything worked
            flatarr.append(flat)
        else: print 'SN too low for ' + objid + '_' + let + num

    else: print 'Cannot Find File for ' + objid + '_' + let + num

#Stack the appended arrays into a single array
modelrat = np.vstack(modelrat)
flatarr = np.vstack(flatarr)

#We only want the median between the good data, so 6500-7500 Angstroms
#Get the indices
goodidx = np.logical_and(wavelength > 6500, wavelength < 7500)
#Find the median of the flat in the right range
medflat = np.median(flat[goodidx])

#Compute the deblazed ratio for each object
deblazarr = divz(modelrat,flatarr/medflat)
deblaze50 = np.median(deblazarr,axis=0)
deblaze16 = np.percentile(deblazarr,16,axis=0)
deblaze84 = np.percentile(deblazarr,84,axis=0)

#Find the median, 16th percentile, and 84th percentile of the models pixel by pixel
perc50 = np.median(modelrat,axis=0)
perc16 = np.percentile(modelrat,16,axis=0)
perc84 = np.percentile(modelrat,84,axis=0)
flat = np.median(flatarr,axis=0)

#Normalize the models first, then find the three percentiles
#Find the median of each model
#modelratmeds = [np.median(i) for i in modelrat]
medlist = []
demedlist = []
for j in range(0,len(modelrat)):
    #Compute the median over only the good pixel range
    medlist.append(np.median(modelrat[j][goodidx]))
    demedlist.append(np.median(deblazarr[j][goodidx]))
#Make the median lists into arrays
meds = np.asarray(medlist)
demeds = np.asarray(demedlist)
#Divide through by the median lists to normalize each line of the array
nmodelrat = np.transpose(np.divide(np.transpose(modelrat),meds))
ndeblaze = np.transpose(np.divide(np.transpose(deblazarr),demeds))
#Find the 50th, 16th, 84th percentile pixel by pixel in each of the two arrays
normperc50 = np.median(nmodelrat,axis=0)
normperc16 = np.percentile(nmodelrat,16,axis=0)
normperc84 = np.percentile(nmodelrat,84,axis=0)
normdeblaze50 = np.median(ndeblaze,axis=0)
normdeblaze16 = np.percentile(ndeblaze,16,axis=0)
normdeblaze84 = np.percentile(ndeblaze,84,axis=0)

#Find their medians
#med50 = np.median(perc50[goodidx])
#med16 = np.median(perc16[goodidx])
#med84 = np.median(perc84[goodidx])
#Divide through by the median
#normperc50 = divz(perc50,med50)
#normperc16 = divz(perc16,med50)
#normperc84 = divz(perc84,med50)

#Deblazing the spectra:

#Divide the particular flat by its median to normalize it
#normflat = divz(flat,medflat)
#Now divide each of the models by the normalized flat
#deblaze50 = divz(perc50,normflat)
#deblaze16 = divz(perc16,normflat)
#deblaze84 = divz(perc84,normflat)
#Finally, normalize them
#meddeblaze = np.median(deblaze50[goodidx])
#normdeblaze50 = divz(deblaze50,meddeblaze)
#normdeblaze16 = divz(deblaze16,meddeblaze)
#normdeblaze84 = divz(deblaze84,meddeblaze)


#Polynomial fitting to the deblazed, normalized spectrum
order = 4
#only fit to data before 9200 angstroms since it is noisy after that
polyidx = np.where(wavelength < 9200)
#Median filter the spectrum to get rid of the giant peaks
snormdeblaze50 = medfilt(normdeblaze50,51)
#Fit a polynomial to the median filtered spectrum,but only in the disred wavelengths
poly = np.polyfit(wavelength[polyidx],snormdeblaze50[polyidx],order)
#Get the output of this polynomial
polygraph = np.polyval(poly,wavelength)


#Plot the data, smooth, model, and fit
fig, axarr = plt.subplots(2,2,figsize=(14,7),sharex=True)
ax0,ax1,ax2,ax3 = axarr[0,0],axarr[0,1],axarr[1,0],axarr[1,1]
ax0.plot(wavelength,perc50, color='red', label = 'Median')
ax0.plot(wavelength,perc84, color='black', label = '84th Percentile')
ax0.plot(wavelength,perc16, color='gray', label = '16th Percentile')
ax1.plot(wavelength,normperc50, color='red', label = 'Median')
ax1.plot(wavelength,normperc84, color='black', label = '84th Percentile')
ax1.plot(wavelength,normperc16, color='gray', label = '16th Percentile')
ax2.plot(wavelength,deblaze50, color='red', label = 'Median')
ax2.plot(wavelength,deblaze84, color='black', label = '84th Percentile')
ax2.plot(wavelength,deblaze16, color='gray', label = '16th Percentile')
ax3.plot(wavelength,normdeblaze50, color='red', label = 'Median')
ax3.plot(wavelength,normdeblaze84, color='black', label = '84th Percentile')
ax3.plot(wavelength,normdeblaze16, color='gray', label = '16th Percentile')
ax3.plot(wavelength,polygraph, color='blue', label = 'Polynomial Fit ' + str(order))
ax3.plot([9200,9200],[-1000,1000],color='blue')
ax0.set_title('Raw Ratio, Mask ' + let + num)
ax1.set_title('Normalized Raw Ratio')
ax2.set_title('Deblazed Ratio')
ax3.set_title('Normalized Deblazed Ratio')
ax0.set_ylim(0,3*np.median(perc50))
ax1.set_ylim(0,3*np.median(normperc50))
ax2.set_ylim(0,3*np.median(deblaze50))
ax3.set_ylim(0,3*np.median(normdeblaze50))
ax3.set_xlabel('Wavelength (Angstroms)')
ax0.legend()
ax3.legend()
plt.tight_layout()
fig.savefig(figout + letnum + '_SEDratio.png')
plt.close()


#Flux Calibrate the objects

i=9

if 1 == 1:
#for i in range(1,30):
    OBJID = ourdata.iloc[i].ImageName[4:10]
    objid = OBJID
    #Read the corresponding SED model, e.g. '026110_a6.dat'
    if os.path.exists(seddatapath + objid + '_' + let + num + '.dat'):
        sedmodel2 = ascii.read(seddatapath + objid + '_' + let + num + '.dat').to_pandas()
        #Read the corresponding spectrum 
        hdu = fits.open(fitsdatapath + ourdata.iloc[i].ImageName)[0]
        head = hdu.header
        fitsdata = hdu.data
        #Calculate the wavelength range for the data
        crval1 = head["crval1"]
        crpix1 = head["crpix1"]
        cdelt1 = head["cdelt1"]
        naxis1 = head["naxis1"]
        dcflag = head["dc-flag"]
        exptime = head['exptime']
        wavelength = (1.0+np.arange(naxis1)-crpix1)*cdelt1 + crval1
        if dcflag: wavelength = np.power(10.0,wavelength)

        #Split the fits data into its compenents
        spec = fitsdata[0]
        orig = fitsdata[1]
        noise = fitsdata[2]
        flat = fitsdata[3]
        sky = fitsdata[4]
        tc = fitsdata[5]
        mask = fitsdata[6]
        fit = fitsdata[7]

        #Flux calibrate the spectrum
        m = np.median(divz(flat[goodidx],medflat))
        corspec = divz(spec/exptime,polygraph*flat*m/medflat)
    
        #Plot the flux-calibrated corrected spectrum for one object
        fig, axarr = plt.subplots(2,2,figsize=(14,8),sharex=True)
        ax0,ax3,ax1,ax2 = axarr[0,0],axarr[0,1],axarr[1,0],axarr[1,1]
        ax0.plot(wavelength,spec, color='cornflowerblue', label = 'Orginial Spectrum')
        #ax0.plot(wavelength,noise, color='noise', label = 'Noise')
        ax1.plot(wavelength,normdeblaze84, color='lightgray', label = '+$\sigma$')
        Briax1.plot(wavelength,normdeblaze16, color='gray', label = '-$\sigma$')
        ax1.plot(wavelength,normdeblaze50, color='Black', label = 'Median')
        ax1.plot(wavelength,polygraph, color='C1', label = 'Polynomial Fit') #, order ' + str(order))
        ax1.plot([9200,9200],[-1000,1000],color='C1',ls='-')
        #ax2.plot(wavelength,corspec,color='cornflowerblue',label='Corrected Spectrum')
        ax3.plot(wavelength,flat,color='cornflowerblue',label='Flat')
        #Plot Dan's flux corrected spectrum if it exists
        flxfits = caldatapath + 'flx_' + OBJID + '_feb1' + num + '_' + let + 'big.fits'
        if os.path.exists(flxfits):
            flxdata = fits.open(flxfits)[0].data
            flxspec = flxdata[0]
            model = flxdata[3]
            ax2.plot(wavelength,flxspec*10**(-17),color='cornflowerblue',label='Flux-Calibrated Spectrum', alpha=1)
            ax2.plot(wavelength,model*10**(-17),color='red',label='Stellar Population Model',alpha=1)
        #ax2.plot(sedmodel2['lambda'],sedmodel2['model'],color='red',label='Model Spectrum',alpha=0.5)
        titlefont = 18
        axisfont = 14
        ticklabel = 12
        ax0.set_title('Original Spectrum ' + OBJID + '_' + letnum,fontsize=titlefont)
        ax1.set_title('Normalized Deblazed Flux Ratio',fontsize=titlefont)
        ax2.set_title('Flux-Calibrated Spectrum',fontsize=titlefont)
        ax3.set_title('Flatfield Spectrum',fontsize=titlefont)
        scale = 2.5
        ax0.set_ylim(0,scale*np.median(spec))
        ax3.set_ylim(0,scale*np.median(flat))
        ax1.set_ylim(0,scale*np.median(normdeblaze50))
        ax2.set_ylim(0,scale*np.median(flxspec*10**(-17)))
        ax2.legend()
        ax1.legend()
        ax1.set_xlabel('Wavelength ($\AA$)',fontsize = axisfont)
        ax2.set_xlabel('Wavelength ($\AA$)',fontsize = axisfont)
        ax0.set_ylabel('Counts',fontsize = axisfont)
        ax3.set_ylabel('Counts',fontsize = axisfont)
        ax1.set_ylabel('Normalized Counts',fontsize = axisfont)
        ax2.set_ylabel('Flux (erg/s/${cm}^2/\AA$)',fontsize = axisfont)
        ax3.set_xlim(4800,11000)
        ax0.tick_params(labelsize=ticklabel)
        ax1.tick_params(labelsize=ticklabel)
        ax2.tick_params(labelsize=ticklabel)
        ax3.tick_params(labelsize=ticklabel)
        plt.tight_layout()
        plt.show()
        fig.savefig(figout + 'flx_' + OBJID + '_' + letnum + '.png')
        plt.close()

