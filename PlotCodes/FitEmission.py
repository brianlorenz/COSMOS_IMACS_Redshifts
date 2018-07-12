#Fits an emission ine with a Gaussian and returns the amplitude, standard deviation, and continuum line

#Usage:      run SED_Ratio.py a6 4861   to fit the lines at rest wavelengths 6563 (Ha) for the a6 mask. 



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
viewdataout = '/Users/blorenz/COSMOS/COSMOSData/lineflux_view.txt'

#Folder to save the figures
figout = '/Users/blorenz/COSMOS/COSMOSData/fitEmissionOut/'

#The location with the file for all of our data
ourdatapath = '/Users/blorenz/COSMOS/COSMOSData/all_c_hasinger.txt'
#Where the calibrated spectra are stored
caldatapath = '/Users/blorenz/COSMOS/COSMOSData/flxFitsFileOut/'
#File for all of the emission/absorption features of the galaxy (to mask out other features when fitting)
linedata = '/Users/blorenz/COSMOS/COSMOSData/corFitsFileOut/galaxylines.dat'
#File for the MAD of the difference in flux of duplicates in each line (to flag low S/N lines)
maddatapath = '/Users/blorenz/COSMOS/COSMOSData/linemad.txt'

#Read in the spectral lines for masking
gallines = ascii.read(linedata).to_pandas()
#Remove all absoption lines
gallines = gallines[gallines.col2==1]
gallines = gallines.reset_index()

#Read in the mad of the lines 
maddata = ascii.read(maddatapath).to_pandas()


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

#Set the letnum 
letnum = sys.argv[1]

#Read in all of our data
ourdata = ascii.read(ourdatapath).to_pandas()
ourdata = ourdata[ourdata.ImageName.str.contains('feb1' + letnum[1] + '_' + letnum[0]) == True]
ourdata = ourdata[ourdata.Unsure == 0]
ourdata = ourdata[ourdata.Bad == 0]
ourdata = ourdata[ourdata.Flag3 == 0]
#ourdata = ourdata[ourdata.Flag1 == 0]
ourdata = ourdata[ourdata.Star == 0]

#Function to make the mask before the gaussian
def getMask(modelspec,sigspec,spectrum):
    #Model continuum
    m = modelspec
    #Find all of the pixels where the flux goes to 0 or negative, and set those to 0
    maskline = (spectrum > 0)
    #Get the weights so we can downweight by noise
    w = divz(1,sigspec)*maskline
    return m,w
    


#Find the objid of every object, and it's corresponding letter number combination
#objs[0] - objid
#objs[1] - letter
#objs[2] - number
objs = [(i[4:10],i[17],i[15]) for i in ourdata.ImageName]

#Start two counters to run along the plot
plt1 = 0
plt10 = 0
#Set the gridsize, so 12 means a 12x12 grid
gridsize = 12
#Start the plot before the loop:
fig,axarr = plt.subplots(gridsize,gridsize,figsize = (150,80))

#Loop the fitting over all objects
#for i in range(16,20):
for i in range(len(objs)):
    #Mark the data as good
    fitflag = 0     #Good data
    #Get the redshift
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
            #for j in range(1, len(sys.argv)):
            #Changed to only fitting one line at a time, don't want to unindent everything
            if 1==1:
                #line = int(sys.argv[j])
                line = int(sys.argv[2])
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
                shnoiseline = noise[idx2]

                #Redshift the lines to the current galaxy
                zgallines = gallines.col1*(1+zcc)

            
                #Mask out the spectral lines with this function
                #data - the data to mask out
                #line - the line to keep (others are masked)
                def droplines(wavedrop=waveline,specdrop=specline,modeldrop=modelline,noisedrop = noiseline,zline=zline,peakwave=0,zcc=zcc):
                    #We first find the line that you are fitting so we don't mask it
                    #Compute the differenc between the current line and every line in the data
                    linediff = zgallines - zline
                    #Find the index of the closest value to 0. There may be negatives
                    closelineidx = np.abs(linediff).idxmin()
                    #Save the name of the line for later
                    linename = gallines.iloc[closelineidx].col3
                    restwave = gallines.iloc[closelineidx].col1
                    #Drop the closest line fromt he table so that we mask the others
                    otherlines = zgallines.drop(closelineidx)
                    #Find the other lines that are around the current line, as integers
                    rounded = [np.round(i) for i in otherlines if (i > zline-srange and i < zline+srange)]
                    #Make them even if they are odd to match up with wavelengths
                    centers = [int(i)+(int(i)&1) for i in rounded]
                    #Find offset from expected
                    lineval = gallines.iloc[closelineidx].col1
                    zlineval = lineval*(1+zcc)
                    if peakwave:
                        waveoffset = peakwave-zline
                        #Round it and make it even
                        waveoffset = np.floor(waveoffset)
                        waveoffset = int(waveoffset)+(int(waveoffset)&1)
                        centers = [i+waveoffset for i in centers]
                    #Arrays for the pixels on either side of each center
                    centerrange = [np.arange(i-shrange,i+shrange+2,2) for i in centers]
                    #Find the indices where the arrays match (we will drop these)
                    dropidx = [np.nonzero(np.in1d(wavedrop,i))[0] for i in centerrange]
                    #Save this version for plotting
                    pdropidx = dropidx 
                    #Drop the values at those indices from both wavelength and spectrum
                    #Fixes a bug when they are not the same length -happens if line is on an edge
                    if len(dropidx) == 2:
                        dropidx = np.append(dropidx[0],dropidx[1])
                    #Drop the lines
                    newwave = np.delete(wavedrop,dropidx)
                    newspec = np.delete(specdrop,dropidx)
                    newmodel = np.delete(modeldrop,dropidx)
                    newnoise = np.delete(noisedrop,dropidx)
                    return newwave,newspec,newmodel,newnoise,dropidx,linename,restwave,pdropidx


                #Mask the other emission lines
                dropwaveline,dropspecline,dropmodelline,dropnoiseline,dropidx,linename,restwave,pdropidx = droplines()

                m,w = getMask(dropmodelline, dropnoiseline, dropspecline)
                #Model continuum
                #m = dropmodelline
                #Get the weights so we can downweight by noise
                #w = divz(1,dropnoiseline)

                #Set up Gaussian Function
                #mu - mean value of the gaussian
                #sigma - standard deviation
                def gauss3(x, mu, sigma):
                    A,B = amp3(x,mu,sigma)
                    g = np.exp(-0.5*(x-mu)**2/(np.e**sigma)**2)/np.sqrt(2*np.pi*(np.e**sigma)**2) #NORMALIZED GAUSSIAN
                    s = A*g + B*m
                    return s

                #A is area under Gauss curve, B is the scale factor of the continuum
                def amp3(x, mu, sigma):
                    g = np.exp(-0.5*(x-mu)**2/(np.e**sigma)**2)/np.sqrt(2*np.pi*(np.e**sigma)**2) #NORMALIZED GAUSSIAN
                    A,B = nnls(np.transpose([g,m])*w[::,np.newaxis],dropspecline*w)[0]
                    return A,B

                ###Set initial guess parameters
                #find the highest peak, get the wavelength value of it
                #Index of highest peak
                pkidx = np.argmax(shspecline)+srange/2-shrange/2
                #Wavelength of peak
                peakwave = waveline[pkidx]
                guess3 = (peakwave,np.log(2))   
                guesscurve3 = gauss3(dropwaveline,guess3[0],guess3[1])   
                #Set the bounds
                bounds3 = ([peakwave-10,np.log(2)],[peakwave+10,np.log(10)])

                
                #Check if there is a lot of bad data
                if np.count_nonzero(~np.isnan(specline)):
                    try:
                        #Fit the Gaussian
                        #coeff3, var_matrix3 = curve_fit(gauss3, waveline, specline, p0=guess3, bounds=bounds3)
                        coeff3, var_matrix3 = curve_fit(gauss3, dropwaveline, dropspecline, p0=guess3, bounds=bounds3)
                        #Fit again with a proper mask
                        #Mask the other emission lines
                        peakwave = coeff3[0]
                        dropwaveline,dropspecline,dropmodelline,dropnoiseline,dropidx,linename,restwave,pdropidx = droplines(peakwave=peakwave)
                        guess3 = (peakwave,coeff3[1])

                        
                        #Redefine the gauss functions since now the model and noise have changed
                        m,w = getMask(dropmodelline, dropnoiseline, dropspecline)
                        #Model continuum
                        #m = dropmodelline
                        #Get the weights so we can downweight by noise
                        #w = divz(1,dropnoiseline)

                        #Set up Gaussian Function
                        #mu - mean value of the gaussian
                        #sigma - log(standard deviation)
                        def gauss3(x, mu, sigma):
                            A,B = amp3(x,mu,sigma)
                            g = np.exp(-0.5*(x-mu)**2/(np.e**sigma)**2)/np.sqrt(2*np.pi*(np.e**sigma)**2) #NORMALIZED GAUSSIAN
                            s = A*g + B*m
                            return s

                        #A is area under Gauss curve, B is the scale factor of the continuum
                        def amp3(x, mu, sigma):
                            g = np.exp(-0.5*(x-mu)**2/(np.e**sigma)**2)/np.sqrt(2*np.pi*(np.e**sigma)**2) #NORMALIZED GAUSSIAN
                            A,B = nnls(np.transpose([g,m])*w[::,np.newaxis],dropspecline*w)[0]
                            return A,B
            
                        coeff3, var_matrix3 = curve_fit(gauss3, dropwaveline, dropspecline, p0=guess3, bounds=bounds3)
                        #Compute the values of the fit
                        gausscurve3 = gauss3(dropwaveline,coeff3[0],coeff3[1])   #
                        amp3 = amp3(dropwaveline,coeff3[0],coeff3[1])   #
                        mu3 = coeff3[0] 
                        stddev3 = np.e**np.abs(coeff3[1])
                        flux3 = amp3[0]
                        scale3 = amp3[1]
                        
                        #Compute chi^2 statistics in the range of the line
                        #Degrees of freedom: mu, sigma, area, scale
                        dof = 4
                        #Set the lower and upper bounds for the region to find chi2
                        chilb = mu3-2*stddev3
                        chiub = mu3+2*stddev3
                        #Get only the indices in that region
                        cidx = np.logical_and(dropwaveline > chilb-2, dropwaveline < chiub+2)
                        arrchi2 = divz((dropspecline[cidx]-gausscurve3[cidx]),dropnoiseline[cidx])**2
                        chi2 = np.add.reduce(arrchi2)
                        rchi2 = divz(chi2,len(dropwaveline[cidx])-dof)

                        #Now compute the weigthed error
                        #Gaussian curve with area=1
                        g = np.exp(-0.5*(dropwaveline-mu3)**2/stddev3**2)/np.sqrt(2*np.pi*stddev3**2) #NORMALIZED GAUSSIAN
                        wsig = np.sqrt(np.sum(g*(dropnoiseline**2))*np.sqrt(2*np.pi*(stddev3**2)))
                        usig = np.sqrt(np.add.reduce(shnoiseline**2))

                        #Get the string of the nearest wavelength to the line. Used for saving everythign
                        linestr = (str(int(np.round(restwave))))
                        
                        #Check if more than half of the spectrum is masked - if so, throw it out
                        if (len(np.where(w<=0)[0])>(len(dropwaveline)/3)):
                            fitflag = 1 #Marks bad data
                        #Check if the width of the line hit the bounds
                        elif ((stddev3 < 2.1) or (stddev3 > 8)):
                            fitflag = 2 #Marks bad sigma
                        #Check if the flux of the detection is above 3sigma over the mad
                        elif (flux3 < 3*maddata[linestr+'_mad'][0]):
                            fitflag = 3 #Marks low S/N
                        
                            


                        def mkplot(plt10,plt1,gridsize):
                            #Create the plot    
                            #fig,ax0 = plt.subplots(figsize = (13,7))
                            #Set the axis to the correct number
                            ax0 = axarr[plt10,plt1]
                            #Increment the counters for next time
                            plt1 = plt1 + 1
                            if plt1 == gridsize:
                                plt1 = 0
                                plt10 = plt10 + 1
                            #Plotting
                            ax0.plot(waveline,specline,color='cornflowerblue',label='Spectrum')
                            #ax0.plot(dropwaveline,dropspecline,color='darkblue',label='Masked Spectrum')
                            #This will break if one of the lines has an empty array, the except statement fixes it. This is only for plotting
                            try: [ax0.axvspan(np.min(waveline[j]),np.max(waveline[j]), color='indianred', alpha=0.1) for j in pdropidx]
                            except: [ax0.axvspan(np.min(waveline[j]),np.max(waveline[j]), color='indianred', alpha=0.1) for j in dropidx]
                            #Check if any weights were set to 0 - if so, plot the mask for those
                            if np.where(w<=0)[0].any():
                                [ax0.plot(dropwaveline[j],dropspecline[j], marker='o', color='red', alpha=0.7) for j in np.where(w<=0)[0]]
                            #Plot the region over which we fit chi2
                            ax0.axvspan(np.min(dropwaveline[cidx]),np.max(dropwaveline[cidx]), color='grey', alpha=0.2, label='chi2 region')
                            ax0.plot(waveline,modelline,color='red',label='Model')
                            #ax0.plot(dropwaveline,guesscurve3,color='orange',label='Initial Guess')
                            ax0.plot(dropwaveline,dropnoiseline,color='orange',label='Noise')
                            #Titles, axes, legends
                            ax0.set_title('Fit for line at rest $\lambda$ of ' + str(int(np.round(restwave))) + ', OBJID ' + objs[i][0] + '_' + objs[i][1] + objs[i][2] + ', z=' + str(np.around(zcc,4)),fontsize = titlefont)
                            ax0.legend(fontsize = legendfont,loc=1)
                            ax0.set_xlabel('Wavelength ($\AA$)',fontsize = axisfont)
                            ax0.set_ylabel('Flux ($10^{-17}$ erg/s/${cm}^2/\AA$)',fontsize = axisfont)
                            ax0.tick_params(labelsize = ticksize)
                            return ax0, plt10, plt1

                        
                            

                        ax0,plt10,plt1 = mkplot(plt10,plt1,gridsize)
                        ax0.text(0.02,0.95,'Mean:       ' + str(round(mu3,2)),fontsize = textfont, transform=ax0.transAxes)      
                        ax0.text(0.02,0.90,'Std Dev:   ' + str(round(stddev3,2)),fontsize = textfont, transform=ax0.transAxes)
                        ax0.text(0.02,0.80,'Scale:       ' + str(round(amp3[1],2)),fontsize = textfont, transform=ax0.transAxes)      
                        ax0.text(0.02,0.85,'Flux:         ' + str(round(amp3[0],2)),fontsize = textfont, transform=ax0.transAxes)
                        ax0.text(0.02,0.75,'Chi2:        ' + str(round(chi2,2)),fontsize = textfont, transform=ax0.transAxes)      
                        ax0.text(0.02,0.70,'rChi2:       ' + str(round(rchi2,2)),fontsize = textfont, transform=ax0.transAxes)
                        ax0.text(0.02,0.65,'wsig:        ' + str(round(wsig,3)),fontsize = textfont, transform=ax0.transAxes)
                        ax0.text(0.02,0.60,'usig:         ' + str(round(usig,3)),fontsize = textfont, transform=ax0.transAxes)
                        if fitflag:
                            ax0.text(0.02,0.55,'flag:          ' + str(fitflag),fontsize = textfont, transform=ax0.transAxes)
                        #fig.text(0.14,0.60,'Redshift:   ' + str(round(zcc,4)),fontsize = textfont)
                        #fig.text(0.14,0.60,'Luminosity (erg/s):   ' + str(round(lumin,2)),fontsize = textfont)      
                        ax0.plot(dropwaveline,gausscurve3,color='black',label='Gausiian fit')
                        ax0.legend(fontsize = legendfont,loc=1)
                        #plt.show()

                        
                        #Store the results to the output array:
                        #First we find the index with a matching objid
                        #midx = np.where((outarr.OBJID.astype(float)-float(objs[i][0])==0) and (outarr.Mask == (objs[i][1]+objs[i][2])))[0]
                        #Get the array of trues and falses where the OBJID and mask both match
                        tfarr = (outarr.OBJID.astype(float)-float(objs[i][0])==0) & (outarr.Mask == (objs[i][1]+objs[i][2]))
                        #Get the index of the matching element
                        midx = outarr.index[tfarr]
                        #We make sure outarr has correct column types
                        if os.path.exists(dataout): 
                            #outarr.OBJID = outarr.OBJID.astype(str)
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
                            
                        #Write in the new info from the fit. outarr.at auto generates new columns if needed
                        outarr.at[midx,linestr + '_mean'] = mu3
                        outarr.at[midx,linestr + '_stddev'] = stddev3
                        outarr.at[midx,linestr + '_flux'] = flux3
                        outarr.at[midx,linestr + '_scale'] = scale3
                        outarr.at[midx,linestr + '_chi2'] = chi2
                        outarr.at[midx,linestr + '_rchi2'] = rchi2
                        outarr.at[midx,linestr + '_wsig'] = wsig
                        outarr.at[midx,linestr + '_usig'] = usig
                        outarr.at[midx,linestr + '_flag'] = fitflag
                        '''
                        Flag values:
                        1 - too many zeros, we threw out the fit
                        2 - sigma >8 or <2.1, so it hit the bounds. 
                        3 - less than 3 sigma above MAD
                        '''
            
                        
                    except (RuntimeError):
                        ax0.text(0.14,0.84,'Fitting Failed',fontsize = textfont, transform=ax0.transAxes)
                        #plt.show()
                else: print('Bad data at ' + str(line) + ', too many NaN. ' + 'flx_' + objs[i][0] + '_feb1' + objs[i][2] + '_' + objs[i][1] + 'big.fits' )
               

    #If not, give an error but continue
    else: print('Could not read file ' + flxfits)

###Write the file
outarr = outarr.sort_values('OBJID')
outarr = outarr.reindex(sorted(outarr.columns), axis=1)
outarr = outarr.fillna(value = -99.999999999999)
#outarr = outarr.drop('Ha_mean',axis=1)
outarr.to_csv(dataout,index=False)


#Save the figure
#plt.show()
plt.tight_layout()
fig.savefig(figout + str(int(np.round(restwave))) + '_' + linename + '_' + letnum + '.pdf')
plt.close()
                        
