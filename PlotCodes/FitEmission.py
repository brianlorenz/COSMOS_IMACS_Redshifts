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
from scipy.optimize import curve_fit,nnls

#Folder to save the figures
figout = '/Users/blorenz/COSMOS/COSMOSData/fitEmissionOut/'

#The location with the file for all of our data
ourdatapath = '/Users/blorenz/COSMOS/COSMOSData/all_c_hasinger.txt'
#Where the calibrated spectra are stored
caldatapath = '/Users/blorenz/COSMOS/COSMOSData/flxFitsFileOut/'
#File for all of the emission/absorption features of the galaxy (to mask out other features when fitting)
linedata = '/Users/blorenz/COSMOS/COSMOSData/corFitsFileOut/galaxylines.dat'

#Read in the spectral lines for masking
gallines = ascii.read(linedata).to_pandas()

#Division function
def divz(X,Y):
        return X/np.where(Y,Y,Y+1)*np.not_equal(Y,0)

#Fontsizes for plotting
axisfont = 18
ticksize = 16
titlefont = 24
legendfont = 16
textfont = 16

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
for i in range(0,10):
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
            noise = flxdata[2] #?
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
                #Set the short range to try to find the peak
                shrange = 6
                #Find the indices to crop the spectra around the line
                idx = np.logical_and(wavelength > zline-srange, wavelength < zline+srange)
                idx2 = np.logical_and(wavelength > zline-shrange, wavelength < zline+shrange)
                #Crop the spectrum to the proper range
                waveline = wavelength[idx]
                specline = spec[idx]
                shspecline = spec[idx2]
                modelline = model[idx]
                noiseline = noise[idx]

                #Set up Gaussian Function
                #mu - mean value of the gaussian
                #sigma - standard deviation
                #Model continuum
                m = modelline
                #Get the weights so we can downweight by noise
                w = divz(1,noiseline)
                
                def gauss3(x, mu, sigma):
                    A,B = amp3(x,mu,sigma)
                    g = np.exp(-0.5*(x-mu)**2/sigma**2)/np.sqrt(2*np.pi*sigma**2) #NORMALIZED GAUSSIAN
                    s = A*g + B*m
                    print A,B
                    return s

                #A is ?, B is the scale factor of the continuum
                def amp3(x, mu, sigma):
                    g = np.exp(-0.5*(x-mu)**2/sigma**2)/np.sqrt(2*np.pi*sigma**2) #NORMALIZED GAUSSIAN
                    A,B = nnls(np.transpose([g,m])*w[::,np.newaxis],specline*w)[0]
                    return A,B

               

                #Definte the amplitude function (we need specline for this)
                #def amp2(x, mu, sigma, y0, y1):
                #    g = np.exp(-(x-mu)**2/(2.*sigma**2))
                #    t = np.add.reduce((specline-y0-y1*x)*g)
                #    b = np.add.reduce(g*g)
                #    if b: A = t/b
                #    else: A = 0
                #    return A
                
                def gauss2(x, mu, sigma, scale):    #
                    g = amp2(x,mu,sigma,scale)*np.exp(-(x-mu)**2/(2.*sigma**2))+scale+modelline
                    return g

                def amp2(x, mu, sigma, scale):     #
                    g = np.exp(-(x-mu)**2/(2.*sigma**2))
                    t = np.add.reduce((specline-modelline-scale)*g)
                    b = np.add.reduce(g*g)
                    if b: A = t/b
                    else: A = 0
                    return A
                
                ###Set initial guess parameters
                #find the highest peak, get the wavelength value of it
                guessmu = waveline[np.argmax(shspecline)+srange/2-shrange/2]
                #guess = (guessmu,2,np.median(spec),0)
                #guesscurve = gauss2(waveline,guess[0],guess[1],guess[2],guess[3])
                #guess = (guessmu,2,0)   #
                #guesscurve = gauss2(waveline,guess[0],guess[1],guess[2])   #
                guess3 = (guessmu,2)   #
                guesscurve3 = gauss3(waveline,guess3[0],guess3[1])   #
                #Set the bounds
                #bounds = ([guessmu-10,1,0,-1],[guessmu+10,10,5,1])
                #bounds = ([guessmu-10,1,-1],[guessmu+10,7,1])   #
                bounds3 = ([guessmu-10,1],[guessmu+10,7])


                #Mask out the spectral lines with this function
                #data - the data to mask out
                #line - the line to keep (others are masked)
                def droplines(wavedrop=waveline,specdrop=specline,zline=zline):
                    #We first find the line that you are fitting so we don't mask it
                    #Compute the differenc between the current line and every line in the data
                    linediff = zgallines - zline
                    #Find the index of the closest value to 0. There may be negatives
                    closelineidx = np.abs(linediff).idxmin()
                    #Drop the closest line fromt he table so that we mask the others
                    otherlines = zgallines.drop(closelineidx)
                    #Find the other lines that are around the current line, as integers
                    rounded = [np.round(i) for i in otherlines if (i > zline-srange and i < zline+srange)]
                    #Make them even if they are odd to match up with wavelengths
                    centers = [int(i)+(int(i)&1) for i in rounded]
                    #Arrays for the pixels on either side of each center
                    centerrange = [np.arange(i-shrange+2,i+shrange,2) for i in centers]
                    #Find the indices where the arrays match (we will drop these)
                    dropidx = [np.nonzero(np.in1d(wavedrop,i))[0] for i in centerrange]
                    #Drop the values at those indices from both wavelength and spectrum
                    newwave = np.delete(wavedrop,dropidx)
                    newspec = np.delete(specdrop,dropidx)
                    return newwave,newspec

                #Redshift the lines to the current galaxy
                zgallines = gallines.col1*(1+zcc)
                dropwaveline,dropspecline = droplines()
                    

                
                #Check if there is a lot of bad data
                if np.count_nonzero(~np.isnan(specline)):
                    try:
                        #Fit the Gaussian
                        #coeff2, var_matrix2 = curve_fit(gauss2, waveline, specline, p0=guess, bounds=bounds)
                        #coeff3, var_matrix3 = curve_fit(gauss3, waveline, specline, p0=guess3, bounds=bounds3)
                        coeff3, var_matrix3 = curve_fit(gauss3, dropwaveline, dropspecline, p0=guess3, bounds=bounds3)
                        #Compute the values of the fit
                        #gausscurve2 = gauss2(waveline,coeff2[0],coeff2[1],coeff2[2],coeff2[3])
                        #amp2 = amp2(waveline,coeff2[0],coeff2[1],coeff2[2],coeff2[3])
                        #gausscurve2 = gauss2(waveline,coeff2[0],coeff2[1],coeff2[2])   #
                        #amp2 = amp2(waveline,coeff2[0],coeff2[1],coeff2[2])   #
                        gausscurve3 = gauss3(waveline,coeff3[0],coeff3[1])   #
                        amp3 = amp3(waveline,coeff3[0],coeff3[1])   #
                        #mu2 = coeff2[0]
                        mu3 = coeff3[0] 
                        #stddev2 = np.abs(coeff2[1])
                        stddev3 = np.abs(coeff3[1])
                        #y0 = coeff2[2]
                        #y1 = coeff2[3]
                        #scale2 = coeff2[2]   #
                        #scale2fac = np.median(divz(modelline+scale2,modelline))
                        print objs[i], coeff3


                        def mkplot():
                            #Create the plot    
                            fig,ax0 = plt.subplots(figsize = (13,7))
                            #Plotting
                            ax0.plot(waveline,specline,color='cornflowerblue',label='Spectrum')
                            ax0.plot(waveline,modelline,color='red',label='Model')
                            ax0.plot(waveline,guesscurve3,color='orange',label='Initial Guess')
                            #Titles, axes, legends
                            ax0.legend(fontsize = legendfont,loc=1)
                            ax0.set_title('Fit for line at rest $\lambda$ of ' + str(line) + ', OBJID ' + objs[i][0] + '_' + objs[i][1] + objs[i][2] ,fontsize = titlefont)
                            ax0.set_xlabel('Wavelength ($\AA$)',fontsize = axisfont)
                            ax0.set_ylabel('Flux ($10^{-17}$ erg/s/${cm}^2/\AA$)',fontsize = axisfont)
                            ax0.tick_params(labelsize = ticksize)
                            return fig,ax0
                            

                        fig,ax0 = mkplot()
                        fig.text(0.14,0.76,'A:             ' + str(round(amp3[0],3)),fontsize = textfont)
                        fig.text(0.14,0.8,'Std Dev:   ' + str(round(stddev3,3)),fontsize = textfont)
                        #fig.text(0.14,0.76,'Continuum: y=' + str(round(y1,3)) + 'x+' + str(round(y0,3)),fontsize = textfont)
                        #fig.text(0.14,0.76,'Scale: ' + str(round(scale2,3)),fontsize = textfont)      #
                        fig.text(0.14,0.72,'Scale:       ' + str(round(amp3[1],3)),fontsize = textfont)      #
                        fig.text(0.14,0.84,'Mean:       ' + str(int(mu3)),fontsize = textfont)      #
                        ax0.plot(waveline,gausscurve3,color='black',label='Gausiian fit')
                        ax0.legend(fontsize = legendfont,loc=1)
                        plt.show()
                        fig.savefig(figout + objs[i][0] + '_' + objs[i][1] + objs[i][2] + '_' + str(line) + '.png')
                        plt.close()
                    except (RuntimeError):
                        fig.text(0.14,0.84,'Fitting Failed',fontsize = textfont)
                        #plt.show()
                        fig.savefig(figout + objs[i][0] + '_' + objs[i][1] + objs[i][2] + '_' + str(line) + '.png')
                        plt.close()
                else: print('Bad data at ' + str(line) + ', too many NaN. ' + 'flx_' + objs[i][0] + '_feb1' + objs[i][2] + '_' + objs[i][1] + 'big.fits' )
               

    #If not, give an error but continue
    else: print('Could not read file ' + flxfits)
