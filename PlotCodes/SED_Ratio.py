#Finds the ratio between our observed spectra and the ULTRAVista ones so that we can flux calibrate the data

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import sys, os, string
import pandas as pd
from astropy.io import fits
from astropy.convolution import convolve, Box1DKernel


#Folder to save the figures
figout = '/Users/blorenz/COSMOS/COSMOSData/SED_Ratio_Images/'

#The location with the file for all of our data
ourdatapath = '/Users/blorenz/COSMOS/COSMOSData/all_c_hasinger.txt'
#Location for data files from the SED fitting code
seddatapath = '/Users/blorenz/COSMOS/COSMOSData/SEDmodels/'
#Location for all of the .fits files
fitsdatapath = '/Users/blorenz/COSMOS/COSMOSData/corFitsFileOut/'

#Read in all of our data
ourdata = ascii.read(ourdatapath).to_pandas()

#Define the wavelengths in anstroms of the bad skyines and telluric correction region
skylines = [5577.0,5890.0,6300,6364]
telluric = [(7584,7650)]

#Smoothing function to clip our emission lines and smooth
def Smooth(y,good,p=50,h=25):
    m = np.zeros(y.shape,y.dtype)
    for j in range(len(y)):
        a,b = np.clip([j-h,j+h+1],0,len(y)-1)
        u = np.compress(good[a:b],y[a:b])
        if len(u):
            if p == 50:
                m[j] = np.median(u)
            else:
                m[j] = np.sort(u)[len(u)*p/100]
    return m

###Smoothing function to exactly match our data to the model
# target - a list of the wavelengths to center each smooth on
# wavelengths - a list of all of the wavelegnths in the data
# counts - a list of the counts of the data
# good - a list of the points that are not bad data or emission lines
# perc - the percentile to smooth to, rnages from 0 to 100. Suggested 50 or 95. 
# width - the size over which to smooth (in pixels)
def Smoothfinal(target,wavelength,counts,good,perc=50,width=25):
    #Make an array of zeros to populate with the smoothed
    m = np.zeros(target.shape,target.dtype)
    #Loop over the indicies of the points in the array
    for j in range(len(target)):
        #Find the indices of elements within range of the target
        idxarr = np.where(np.logical_and(wavelength>=(target.iloc[j]-width), wavelength<=(target.iloc[j]+width)))
        #Adds together all points within the wavelength range that are not bad data
        u = np.compress(good[np.min(idxarr):np.max(idxarr)],counts[np.min(idxarr):np.max(idxarr)])
        #Checks if u exists. If so, take the median if percentile is 50, otherwise compute the data point at that percentile.
        if len(u):
            if perc == 50:
                m[j] = np.median(u)
            else:
                m[j] = np.sort(u)[len(u)*perc/100]
    #Output the list of new smoothed y values that correspond to the input target
    return m

#Division function
def divz(X,Y):
        return X/np.where(Y,Y,Y+1)*np.not_equal(Y,0)

#Function to fit a polynomial to data
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


i = 10

#Change to for loop to cycle through the data
if 1 == 1:
    #Find the identifying pieces of the string. e.g. 'a' '6' '026110'
    let = ourdata.iloc[i].ImageName[17]
    num = ourdata.iloc[i].ImageName[15]
    objid = ourdata.iloc[i].ImageName[4:10]

    #Read the corresponding SED model, e.g. '026110_a6.dat'
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
    
    #Smooth the spectrum
    bad = 1.*np.logical_or(np.less_equal(spec,1e-4),np.less_equal(noise,1e-4))
    bad = 1.*np.logical_or(np.less_equal(spec,-3*noise),np.less_equal(noise,1e-4))
    for skyline in skylines:
        bad = 1.*np.logical_or(bad,np.greater(wavelength,skyline-3*cdelt1)*np.less(wavelength,skyline+3*cdelt1))
    for absorb in telluric:
        bad = 1.*np.logical_or(bad,np.greater(wavelength,absorb[0])*np.less(wavelength,absorb[-1]))
    bad = 1.*np.greater(convolve(bad,Box1DKernel(8)),0) #astropy boxcar
    sm = Smooth(spec,np.logical_not(bad))
    sr = spec-sm
    ss = 1.49*Smooth(abs(sr),np.logical_not(bad),h=50)
    good = np.logical_not(bad)
    good = good*np.less(sr,2.5*ss) # CLIPS OUT EMISSION LINES AND BAD POSITIVE SKY LINE RESIDUALS

    #This is where the final smooth happens - change h to change the level to which the spectrum is smoothed
    smoothspec = Smooth(spec,np.logical_not(bad),h=20)


    idxarr = np.where(np.logical_and(sedmodel['lambda']>=np.min(wavelength), sedmodel['lambda']<=np.max(wavelength)))
    sedmodel = sedmodel.truncate(np.min(idxarr),np.max(idxarr))
    smoothfinal = Smoothfinal(sedmodel['lambda'],wavelength,spec,np.logical_not(bad),width=60,perc=50)


    #Polynomial fitting
    #order - order of the polynomial
    order = 30
    x = np.arange(0,smoothfinal.shape[0],1,float)
    norm = (x-(x[0]+x[-1])/2)/ ((x[-1]-x[0])/2)
    basis = np.asarray([np.power(norm,i) for i in range(order)])
    sol = svdfit(np.transpose(basis)*self.mask[::,np.newaxis],self.obj*self.mask)
    smoothfit = svdfit()

    
    #Plot the data, smooth, model, and fit
    fig, axarr = plt.subplots(4,1,figsize=(14,8))
    ax0,ax1,ax2,ax3 = axarr[0],axarr[1],axarr[2],axarr[3]
    ax0.plot(wavelength,spec, color='cornflowerblue',label='Spectrum')
    ax0.plot(wavelength,noise,color='orange',label='Noise')
    ax1.plot(wavelength,smoothspec, color='cornflowerblue')
    ax2.plot(sedmodel['lambda'],smoothfinal, color='cornflowerblue')
    ax3.plot(sedmodel['lambda'],sedmodel.model, color='cornflowerblue')
    ax0.set_ylim(ax1.get_ylim())
    ax2.set_ylim(ax1.get_ylim())
    ax3.set_xlim(ax1.get_xlim())
    ax0.set_title('Original Spectrum, ' + ourdata.iloc[i].ImageName)
    ax1.set_title('Old Smooth Spectrum (High Res)')
    ax2.set_title('New Smooth Spectrum (Lower Res)')
    ax3.set_title('SEDmodel Spectrum')
    ax0.legend()
    plt.tight_layout()

    #Find the relevant info for saving
    OBJID = ourdata.iloc[i].ImageName[4:10]
    let = ourdata.iloc[i].ImageName[17]
    num = ourdata.iloc[i].ImageName[15]
    fig.savefig(figout + OBJID + '_' + let + num + '_SEDratio.png')
    plt.show()

    #We have a smooth spectrum, but now maybe we want data points at only the places where the model has data? Need to figure out whether to do this or not. Then we have to find the ratio between the smooth spectrum and the model
