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

#Location of output data file
dataout = '/Users/blorenz/COSMOS/COSMOSData/lineflux.txt'

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

#Read the datafile (if there is one), then create a blank one to write to:
if os.path.exists(dataout):
    outarr = ascii.read(dataout).to_pandas()
else: outarr = pd.DataFrame()
    

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


#Loop the fitting over all objects
for i in range(15,17):
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
            noise = flxdata[1] #?
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

                #Mask out the spectral lines with this function
                #data - the data to mask out
                #line - the line to keep (others are masked)
                def droplines(wavedrop=waveline,specdrop=specline,modeldrop=modelline,noisedrop = noiseline,zline=zline):
                    #We first find the line that you are fitting so we don't mask it
                    #Compute the differenc between the current line and every line in the data
                    linediff = zgallines - zline
                    #Find the index of the closest value to 0. There may be negatives
                    closelineidx = np.abs(linediff).idxmin()
                    #Save the name of the line for later
                    linename = gallines.iloc[closelineidx].col3
                    #Drop the closest line fromt he table so that we mask the others
                    otherlines = zgallines.drop(closelineidx)
                    #Find the other lines that are around the current line, as integers
                    rounded = [np.round(i) for i in otherlines if (i > zline-srange and i < zline+srange)]
                    #Make them even if they are odd to match up with wavelengths
                    centers = [int(i)+(int(i)&1) for i in rounded]
                    #Arrays for the pixels on either side of each center
                    centerrange = [np.arange(i-shrange,i+shrange+2,2) for i in centers]
                    #Find the indices where the arrays match (we will drop these)
                    dropidx = [np.nonzero(np.in1d(wavedrop,i))[0] for i in centerrange]
                    #Drop the values at those indices from both wavelength and spectrum
                    newwave = np.delete(wavedrop,dropidx)
                    newspec = np.delete(specdrop,dropidx)
                    newmodel = np.delete(modeldrop,dropidx)
                    newnoise = np.delete(noisedrop,dropidx)
                    return newwave,newspec,newmodel,newnoise,dropidx,linename


                #Redshift the lines to the current galaxy
                zgallines = gallines.col1*(1+zcc)
                dropwaveline,dropspecline,dropmodelline,dropnoiseline,dropidx,linename = droplines()

                #Model continuum
                m = dropmodelline
                #Get the weights so we can downweight by noise
                w = divz(1,dropnoiseline)

                #Set up Gaussian Function
                #mu - mean value of the gaussian
                #sigma - standard deviation
                def gauss3(x, mu, sigma):
                    A,B = amp3(x,mu,sigma)
                    g = np.exp(-0.5*(x-mu)**2/sigma**2)/np.sqrt(2*np.pi*sigma**2) #NORMALIZED GAUSSIAN
                    s = A*g + B*m
                    return s

                #A is area under Gauss curve, B is the scale factor of the continuum
                def amp3(x, mu, sigma):
                    g = np.exp(-0.5*(x-mu)**2/sigma**2)/np.sqrt(2*np.pi*sigma**2) #NORMALIZED GAUSSIAN
                    A,B = nnls(np.transpose([g,m])*w[::,np.newaxis],dropspecline*w)[0]
                    return A,B

                ###Set initial guess parameters
                #find the highest peak, get the wavelength value of it
                guessmu = waveline[np.argmax(shspecline)+srange/2-shrange/2]
                guess3 = (guessmu,2)   
                guesscurve3 = gauss3(dropwaveline,guess3[0],guess3[1])   
                #Set the bounds
                bounds3 = ([guessmu-10,1],[guessmu+10,7])

                
                #Check if there is a lot of bad data
                if np.count_nonzero(~np.isnan(specline)):
                    try:
                        #Fit the Gaussian
                        #coeff3, var_matrix3 = curve_fit(gauss3, waveline, specline, p0=guess3, bounds=bounds3)
                        coeff3, var_matrix3 = curve_fit(gauss3, dropwaveline, dropspecline, p0=guess3, bounds=bounds3)
                        #Compute the values of the fit
                        gausscurve3 = gauss3(dropwaveline,coeff3[0],coeff3[1])   #
                        amp3 = amp3(dropwaveline,coeff3[0],coeff3[1])   #
                        mu3 = coeff3[0] 
                        stddev3 = np.abs(coeff3[1])
                        flux3 = amp3[0]
                        scale3 = amp3[1]
                        
                        #Compute chi^2 statistics in the range of the line
                        #Degrees of freedom: mu, sigma, area, scale
                        dof = 4
                        #Set the lower and upper bounds for the region to find chi2
                        chilb = mu3-3*stddev3
                        chiub = mu3+3*stddev3
                        #Get only the indices in that region
                        cidx = np.logical_and(dropwaveline > chilb, dropwaveline < chiub)
                        arrchi2 = divz((dropspecline[cidx]-gausscurve3[cidx]),dropnoiseline[cidx])**2
                        chi2 = np.add.reduce(arrchi2)
                        rchi2 = divz(chi2,len(dropwaveline[cidx])-dof)


                        def mkplot():
                            #Create the plot    
                            fig,ax0 = plt.subplots(figsize = (13,7))
                            #Plotting
                            ax0.plot(waveline,specline,color='cornflowerblue',label='Spectrum')
                            #ax0.plot(dropwaveline,dropspecline,color='darkblue',label='Masked Spectrum')
                            #plt.axvspan(8920,8940, color='grey', alpha=0.5)
                            [ax0.axvspan(np.min(waveline[j]),np.max(waveline[j]), color='indianred', alpha=0.1) for j in dropidx]
                            ax0.plot(waveline,modelline,color='red',label='Model')
                            #ax0.plot(dropwaveline,guesscurve3,color='orange',label='Initial Guess')
                            ax0.plot(dropwaveline,dropnoiseline,color='orange',label='Noise')
                            #Titles, axes, legends
                            ax0.set_title('Fit for line at rest $\lambda$ of ' + str(line) + ', OBJID ' + objs[i][0] + '_' + objs[i][1] + objs[i][2] + ', z=' + str(np.around(zcc,4)),fontsize = titlefont)
                            ax0.legend(fontsize = legendfont,loc=1)
                            ax0.set_xlabel('Wavelength ($\AA$)',fontsize = axisfont)
                            ax0.set_ylabel('Flux ($10^{-17}$ erg/s/${cm}^2/\AA$)',fontsize = axisfont)
                            ax0.tick_params(labelsize = ticksize)
                            return fig,ax0

                        
                            

                        fig,ax0 = mkplot()
                        fig.text(0.14,0.76,'Flux:         ' + str(round(amp3[0],2)),fontsize = textfont)
                        fig.text(0.14,0.8,'Std Dev:   ' + str(round(stddev3,2)),fontsize = textfont)
                        fig.text(0.14,0.72,'Scale:       ' + str(round(amp3[1],2)),fontsize = textfont)      
                        fig.text(0.14,0.84,'Mean:       ' + str(round(mu3,2)),fontsize = textfont)      
                        fig.text(0.14,0.68,'Chi2:         ' + str(round(chi2,2)),fontsize = textfont)      
                        fig.text(0.14,0.64,'rChi2:       ' + str(round(rchi2,2)),fontsize = textfont)
                        #fig.text(0.14,0.60,'Redshift:   ' + str(round(zcc,4)),fontsize = textfont)
                        #fig.text(0.14,0.60,'Luminosity (erg/s):   ' + str(round(lumin,2)),fontsize = textfont)      
                        ax0.plot(dropwaveline,gausscurve3,color='black',label='Gausiian fit')
                        #These lines plot where the chi squared is fitting
                        #ax0.plot((np.min(dropwaveline[cidx]),np.min(dropwaveline[cidx])),(-.5,.5),color='black',label='Gausiian fit')
                        #ax0.plot((np.max(dropwaveline[cidx]),np.max(dropwaveline[cidx])),(-.5,.5),color='black',label='Gausiian fit')
                        ax0.legend(fontsize = legendfont,loc=1)
                        #plt.show()
                        fig.savefig(figout + objs[i][0] + '_' + objs[i][1] + objs[i][2] + '_' + str(line) + '.png')
                        plt.close()
                        #Store the results to the output array:
                        #First we find the index with a matching objid
                        midx = np.where(outarr.OBJID.astype(float)-float(objs[i][0])==0)[0]
                        #We make sure outarr has correct column types
                        if os.path.exists(dataout): 
                            outarr.OBJID = outarr.OBJID.astype(str)
                            outarr.Mask = outarr.Mask.astype(str)
                            outarr.fluxfile = outarr.fluxfile.astype(str)
                        #We check to make sure there is only one.
                        #If there are none, we append a new row onto outarr
                        if len(midx)>1:
                            print('Error, check text document for duplicates')
                        elif len(midx)==0:
                            #Makes the index the length of the array, which will add a new row at the bottom
                            midx = len(outarr)
                            #Store the info that doesn't change
                            outarr.at[midx,'OBJID'] = objs[i][0]
                            outarr.at[midx,'Mask'] = objs[i][1]+objs[i][2]
                            outarr.at[midx,'fluxfile'] = 'flx_' + objs[i][0] + '_feb1' + objs[i][2] + '_' + objs[i][1] + 'big.fits'
                            outarr.at[midx,'zcc'] = zcc
                        #Write in the new info from the fit
                        if (linename == 'Halpha'):
                            outarr.at[midx,'Ha_mean'] = mu3
                            outarr.at[midx,'Ha_stddev'] = stddev3
                            outarr.at[midx,'Ha_flux'] = flux3
                            outarr.at[midx,'Ha_scale'] = scale3
                            outarr.at[midx,'Ha_chi2'] = chi2
                            outarr.at[midx,'Ha_rchi2'] = rchi2
                        elif (linename == 'Hbeta'):
                            outarr.at[midx,'HB_mean'] = mu3
                            outarr.at[midx,'HB_stddev'] = stddev3
                            outarr.at[midx,'HB_flux'] = flux3
                            outarr.at[midx,'HB_scale'] = scale3
                            outarr.at[midx,'HB_chi2'] = chi2
                            outarr.at[midx,'HB_rchi2'] = rchi2
                        else: print('Could not find line name, writing nothing. Check gallines')
                        
    
                        
                        
                    except (RuntimeError):
                        fig.text(0.14,0.84,'Fitting Failed',fontsize = textfont)
                        #plt.show()
                        fig.savefig(figout + objs[i][0] + '_' + objs[i][1] + objs[i][2] + '_' + str(line) + '.png')
                        plt.close()
                else: print('Bad data at ' + str(line) + ', too many NaN. ' + 'flx_' + objs[i][0] + '_feb1' + objs[i][2] + '_' + objs[i][1] + 'big.fits' )
               

    #If not, give an error but continue
    else: print('Could not read file ' + flxfits)

###Write the file
outarr = outarr.sort_values('OBJID')
outarr = outarr.reindex(sorted(outarr.columns), axis=1)
outarr.to_csv(dataout,index=False) 


#Need to figure out how to output the data table outarr
#Also need to figure out how ot add strings to pandas array
#Then take the output file and see if it can be read as input
