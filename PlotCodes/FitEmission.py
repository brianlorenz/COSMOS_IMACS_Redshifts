#Fits an emission ine with a Gaussian and returns the amplitude, standard deviation, and continuum line

#Usage:      run SED_Ratio.py 4861 6563 [...]      to fit the lines at rest wavelengths 4861,6563 (Ha and HB) for the a6 mask. 

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import sys, os, string
import pandas as pd
from astropy.io import fits
from astropy.convolution import convolve, Box1DKernel
from scipy.interpolate import splrep, splev
from scipy.signal import medfilt
from scipy.optimize import curve_fit

#Folder to save the figures
figout = '/Users/blorenz/COSMOS/COSMOSData/fitEmissionOut/'

#The location with the file for all of our data
ourdatapath = '/Users/blorenz/COSMOS/COSMOSData/all_c_hasinger.txt'
#Where the calibrated spectra are stored
caldatapath = '/Users/blorenz/COSMOS/COSMOSData/flxFitsFileOut/'

#Read in all of our data
ourdata = ascii.read(ourdatapath).to_pandas()
ourdata = ourdata[ourdata.Unsure == 0]
ourdata = ourdata[ourdata.Bad == 0]
ourdata = ourdata[ourdata.Flag3 == 0]
#ourdata = ourdata[ourdata.Flag1 == 0]
ourdata = ourdata[ourdata.Star == 0]

#Find the objid of every object, and it's corresponding letter number combination
#objs[0] - objid
#objs[1] - letter
#objs[2] - number
objs = [(i[4:10],i[17],i[15]) for i in ourdata.ImageName]

#Set up Gaussian Function
#mu - mean value of the gaussian
#sigma - standard deviation
#y0 - constant of continuum line
#y1 - slope of continuum line
#def gauss2(x, mu, sigma, y0, y1):
#    g = amp2(x,mu,sigma,y0,y1)*np.exp(-(x-mu)**2/(2.*sigma**2))+y0+x*y1
#    return g



#Loop the fitting over all objects
for i in range(0,1):
#for i in range(len(objs)):
    zcc = ourdata.iloc[i].z_cc
    #Set the location of the data file
    flxfits = caldatapath + 'flx_' + objs[i][0] + '_feb1' + objs[i][2] + '_' + objs[i][1] + 'big.fits'
    #Read in its datafile if it exists
    if os.path.exists(flxfits):
            flxdata = fits.open(flxfits)[0].data
            flxhead = fits.open(flxfits)[0].header
            #Read in the spectrum and model
            spec = flxdata[0]
            model = flxdata[3]
            #Calculate the wavelength range for the data
            crval1 = flxhead["crval1"]
            crpix1 = flxhead["crpix1"]
            cdelt1 = flxhead["cdelt1"]
            naxis1 = flxhead["naxis1"]
            dcflag = flxhead["dc-flag"]
            exptime = flxhead['exptime']
            wavelength = (1.0+np.arange(naxis1)-crpix1)*cdelt1 + crval1

            #Loop over all of the emission lines to fit:
            for j in range(1, len(sys.argv)):

                line = int(sys.argv[j])
                #Compute the wavelength of the line redshifted to the galaxy
                zline = (1+zcc)*line
                #Set the range over which to look for the line (in angstroms, each pixel is 2A)
                srange = 50
                #Find the indices to crop the spectra around the line
                idx = np.logical_and(wavelength > zline-srange, wavelength < zline+srange)
                #Crop the spectrum to the proper range
                waveline = wavelength[idx]
                specline = spec[idx]
                modelline = model[idx]
                
                #The next line should be commented out, introduces a linear term to the continuum for testing
                #specline = specline*np.arange(len(specline))/50

                #Set up Gaussian Function
                #mu - mean value of the gaussian
                #sigma - standard deviation
                #scale - scaling factor to move the model continuum up to the level of the spectrum
                def gauss2(x, mu, sigma, scale):    #
                    g = amp2(x,mu,sigma,scale)*np.exp(-(x-mu)**2/(2.*sigma**2))+scale*modelline
                    return g

                #Definte the amplitude function (we need specline for this)
                #def amp2(x, mu, sigma, y0, y1):
                #    g = np.exp(-(x-mu)**2/(2.*sigma**2))
                #    t = np.add.reduce((specline-y0-y1*x)*g)
                #    b = np.add.reduce(g*g)
                #    if b: A = t/b
                #    else: A = 0
                #    return A

                def amp2(x, mu, sigma, scale):     #
                    g = np.exp(-(x-mu)**2/(2.*sigma**2))
                    t = np.add.reduce((specline-modelline)*g)
                    b = np.add.reduce(g*g)
                    if b: A = t/b
                    else: A = 0
                    return A
                
                ###Set initial guess parameters
                #find the highest peak, get the wavelength value of it
                guessmu = waveline[np.argmax(specline)]
                #guess = (guessmu,2,np.median(spec),0)
                #guesscurve = gauss2(waveline,guess[0],guess[1],guess[2],guess[3])
                guess = (guessmu,2,1)   #
                guesscurve = gauss2(waveline,guess[0],guess[1],guess[2])   #
                print guess
                #Set the bounds
                #bounds = ([guessmu-10,1,0,-1],[guessmu+10,10,5,1])
                bounds = ([guessmu-10,1,0],[guessmu+10,10,2])   #

                
                #Check if there is a lot of bad data
                if np.count_nonzero(~np.isnan(specline)):
                    try:
                        #Fit the Gaussian
                        coeff2, var_matrix2 = curve_fit(gauss2, waveline, specline, p0=guess, bounds=bounds)
                        #Compute the values of the fit
                        #gausscurve2 = gauss2(waveline,coeff2[0],coeff2[1],coeff2[2],coeff2[3])
                        #amp2 = amp2(waveline,coeff2[0],coeff2[1],coeff2[2],coeff2[3])
                        gausscurve2 = gauss2(waveline,coeff2[0],coeff2[1],coeff2[2])   #
                        amp2 = amp2(waveline,coeff2[0],coeff2[1],coeff2[2])   #
                        mu2 = coeff2[0] 
                        stddev2 = np.abs(coeff2[1])
                        #y0 = coeff2[2]
                        #y1 = coeff2[3]
                        scale = coeff2[2]   #
                        print objs[i], coeff2


                        def mkplot():
                            

                        #Create the plot
                        fig,ax0 = plt.subplots(figsize = (13,7))
                        #Plotting
                        ax0.plot(waveline,specline,color='cornflowerblue',label='Spectrum')
                        ax0.plot(waveline,modelline,color='red',label='Model')
                        ax0.plot(waveline,gausscurve2,color='black',label='Gausiian fit')
                        ax0.plot(waveline,guesscurve,color='orange',label='Initial Guess')
                        #Fontsizes
                        axisfont = 18
                        ticksize = 16
                        titlefont = 24
                        legendfont = 16
                        textfont = 16
                        #Titles, axes, legends
                        ax0.legend(fontsize = legendfont,loc=1)
                        ax0.set_title('Fit for line at rest $\lambda$ of ' + str(line) + ', OBJID ' + str(objs[i][0]),fontsize = titlefont)
                        ax0.set_xlabel('Wavelength ($\AA$)',fontsize = axisfont)
                        ax0.set_ylabel('Flux ($10^{-17}$ erg/s/${cm}^2/\AA$)',fontsize = axisfont)
                        ax0.tick_params(labelsize = ticksize)
                        fig.text(0.14,0.84,'Amplitude: ' + str(round(amp2,3)),fontsize = textfont)
                        fig.text(0.14,0.8,'Std Dev: ' + str(round(stddev2,3)),fontsize = textfont)
                        #fig.text(0.14,0.76,'Continuum: y=' + str(round(y1,3)) + 'x+' + str(round(y0,3)),fontsize = textfont)
                        fig.text(0.14,0.76,'Scale: ' + str(round(scale,3)),fontsize = textfont)      #
                        plt.show()
                        fig.savefig(figout + objs[i][0] + '_' + objs[i][1] + objs[i][2] + '_' + str(line) + '.png')
                    except (RuntimeError):
                        fig,ax0 = plt.subplots(figsize = (13,7))
                        #Plotting
                        ax0.plot(waveline,specline,color='cornflowerblue',label='Spectrum')
                        ax0.plot(waveline,modelline,color='red',label='Model')
                        ax0.plot(waveline,guesscurve,color='orange',label='Initial Guess')
                        #Fontsizes
                        axisfont = 18
                        ticksize = 16
                        titlefont = 24
                        legendfont = 16
                        textfont = 16
                        #Titles, axes, legends
                        ax0.legend(fontsize = legendfont,loc=1)
                        ax0.set_title('Fit for line at rest $\lambda$ of ' + str(line) + ', OBJID ' + str(objs[i][0]),fontsize = titlefont)
                        ax0.set_xlabel('Wavelength ($\AA$)',fontsize = axisfont)
                        ax0.set_ylabel('Flux ($10^{-17}$ erg/s/${cm}^2/\AA$)',fontsize = axisfont)
                        ax0.tick_params(labelsize = ticksize)
                        fig.text(0.14,0.84,'Fitting Failed',fontsize = textfont)
                        plt.show()
                        fig.savefig(figout + objs[i][0] + '_' + objs[i][1] + objs[i][2] + '_' + str(line) + '.png')
                else: print('Bad data at ' + str(line) + ', too many NaN. ' + 'flx_' + objs[i][0] + '_feb1' + objs[i][2] + '_' + objs[i][1] + 'big.fits' )
               

    #If not, give an error but continue
    else: print('Could not read file ' + flxfits)


