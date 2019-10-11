#Fits an emission ine with a Gaussian and returns the amplitude, standard deviation, and continuum line

#Usage:      run FitEmission.py 'a6' 4861   to fit the lines at rest wavelengths 6563 (Ha) for the a6 mask.
#Typing  run FitEmission.py 'a6' 'HaNII' will fit all three lines around Ha simulaaneously  



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
dataout = '/Users/blorenz/COSMOS/COSMOSData/lineflux_tot.txt'
viewdataout = '/Users/blorenz/COSMOS/COSMOSData/lineflux_view.txt'

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
#Remove all absoption lines
gallines = gallines[gallines.col2==1]
gallines = gallines.reset_index()


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
#ourdata = ourdata[ourdata.Bad == 0]
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
plt1b = 0
plt10b = 0
#Set the gridsize, so 12 means a 12x12 grid
gridsize = 12
#Start the plot before the loop:
fig,axarr = plt.subplots(gridsize,gridsize,figsize = (100,80))
figb,axarrb = plt.subplots(gridsize,gridsize,figsize = (100,80))

#Loop the fitting over all objects
#for i in range(16,20):
for i in range(len(objs)):
    #Mark the data as good
    fitflag = 0     #Good data
    #Set that we are not looking at the lines around Ha
    HaNII = False
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
            noise = flxdata[1] 
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
                line = sys.argv[2]
                #Check if we are fitting the Ha and NII lines toether:
                if line == 'HaNII':
                    line = 6563
                    #Variable to know that we are fitting three lines
                    HaNII = True
                    #Dataframe that we will store everything in
                    HaNIIdat = pd.DataFrame()
                    #Set up the rest wavelengths for the lines
                    HaNIIdat.at[0,'restwave'] = 6548.1
                    HaNIIdat.at[1,'restwave'] = 6562.8
                    HaNIIdat.at[2,'restwave'] = 6583.0
                    
                else: line = int(line)
                #Compute the wavelength of the line redshifted to the galaxy
                zline = (1+zcc)*line
                #Set the range over which to look for the line (in angstroms, each pixel is 2A)
                srange = 150
                #Set the short range to try to find the peak
                shrange = 6
                #Find the indices to crop the spectra around the line
                idx = np.logical_and(wavelength > zline-srange, wavelength < zline+srange)
                idx2 = np.logical_and(wavelength > zline-shrange, wavelength < zline+shrange)
                #Special case for OII doublet if it isn't redshifted into view:
                if zline < 4910:
                    idx = np.arange(0,srange)
                    idx2 = np.arange(0,shrange)
                    fitflag = 5 #Flagged for not in view
                    
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
                def droplines(wavedrop=waveline,specdrop=specline,modeldrop=modelline,noisedrop = noiseline,zline=zline,peakwave=0,zcc=zcc,HaNII = HaNII):
                    #Mark that we plot the dropped region
                    pdrop = 1
                    #We first find the line that you are fitting so we don't mask it
                    #Compute the differenc between the current line and every line in the data
                    linediff = zgallines - zline
                    #Find the index of the closest value to 0. There may be negatives
                    closelineidx = np.abs(linediff).idxmin()
                    #Save the name of the line for later
                    linename = gallines.iloc[closelineidx].col3
                    restwave = gallines.iloc[closelineidx].col1
                    #Drop the closest line from the table so that we mask the others
                    otherlines = zgallines.drop(closelineidx)
                    #Special case for OII doublet, since it should find 3726.2, then also drop 3728.9
                    if linename == '[OII]':
                        otherlines = otherlines.drop(closelineidx+1)
                        restwave = 3727
                    #Special case for Ha three lines, since it should find Ha, then also drop NII on either side of it
                    if HaNII:
                        otherlines = otherlines.drop(closelineidx-1)
                        otherlines = otherlines.drop(closelineidx+1)
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
                    elif not dropidx:
                        #Variable to say whether or not to plot the dropidx
                        pdrop = 0
                    #Drop the lines
                    newwave = np.delete(wavedrop,dropidx)
                    newspec = np.delete(specdrop,dropidx)
                    newmodel = np.delete(modeldrop,dropidx)
                    newnoise = np.delete(noisedrop,dropidx)
                    return newwave,newspec,newmodel,newnoise,dropidx,linename,restwave,pdropidx,pdrop


                #Mask the other emission lines
                dropwaveline,dropspecline,dropmodelline,dropnoiseline,dropidx,linename,restwave,pdropidx,pdrop = droplines()

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

                def gaussHa(x, z, sigma48, sigma63, sigma83):
                    A48,A63,A83,B  = ampHa(x, z, sigma48, sigma63, sigma83)
                    g48 = np.exp(-0.5*(x-(6548.1*(1+z)))**2/(np.e**sigma48)**2)/np.sqrt(2*np.pi*(np.e**sigma48)**2)
                    g63 = np.exp(-0.5*(x-(6562.8*(1+z)))**2/(np.e**sigma63)**2)/np.sqrt(2*np.pi*(np.e**sigma63)**2)
                    g83 = np.exp(-0.5*(x-(6583.0*(1+z)))**2/(np.e**sigma83)**2)/np.sqrt(2*np.pi*(np.e**sigma83)**2)
                    s = A48*g48 + A63*g63 + A83*g83 + B*m
                    return s

                #A is area under Gauss curve, B is the scale factor of the continuum
                def ampHa(x, z, sigma48, sigma63, sigma83):
                    g48 = np.exp(-0.5*(x-(6548.1*(1+z)))**2/(np.e**sigma48)**2)/np.sqrt(2*np.pi*(np.e**sigma48)**2)
                    g63 = np.exp(-0.5*(x-(6562.8*(1+z)))**2/(np.e**sigma63)**2)/np.sqrt(2*np.pi*(np.e**sigma63)**2)
                    g83 = np.exp(-0.5*(x-(6583.0*(1+z)))**2/(np.e**sigma83)**2)/np.sqrt(2*np.pi*(np.e**sigma83)**2)
                    A48,A63,A83,B = nnls(np.transpose([g48,g63,g83,m])*w[::,np.newaxis],dropspecline*w)[0]
                    return A48,A63,A83,B 

                ###Set initial guess parameters
                #find the highest peak, get the wavelength value of it
                #Index of highest peak
                pkidx = np.argmax(shspecline)+srange/2-shrange/2
                #Wavelength of peak
                peakwave = waveline[pkidx]
                guess3 = (peakwave,np.log(2))   
                guesscurve3 = gauss3(dropwaveline,guess3[0],guess3[1])
                #Set the bounds, from expected position of the line +- 4 pixels, and sigma from 2 to 10
                bounds3 = ([restwave*(1+zcc)-8,np.log(2)],[restwave*(1+zcc)+8,np.log(10)])
                #Special case for OII doublet
                if (line == 3727):
                    guess3 = (peakwave,np.log(4))
                    #Set the bounds
                    bounds3 = ([peakwave-8,np.log(2)],[peakwave+8,np.log(15)])
                    guesscurve3 = gauss3(dropwaveline,guess3[0],guess3[1])   
                    #Special case for Ha lines, need to set for all three gaussians
                if HaNII:
                    guessHa = (zcc,np.log(2),np.log(2),np.log(2))
                    guesscurveHa = gaussHa(dropwaveline,guessHa[0],guessHa[1],guessHa[2],guessHa[3])
                    boundsHa =  ([zcc-0.0012,np.log(2),np.log(2),np.log(2)],[zcc+0.0012,np.log(10),np.log(10),np.log(10)])
                    

                
                #Check if there is a lot of bad data
                if np.count_nonzero(~np.isnan(specline)):
                    try:
                        #Fit the Gaussian
                        #coeff3, var_matrix3 = curve_fit(gauss3, waveline, specline, p0=guess3, bounds=bounds3)
                        if not HaNII:
                            coeff3, var_matrix3 = curve_fit(gauss3, dropwaveline, dropspecline, p0=guess3, bounds=bounds3)
                        else:
                            coeffHa, var_matrixHa = curve_fit(gaussHa, dropwaveline, dropspecline, p0=guessHa, bounds=boundsHa)
                        #Fit again with a proper mask
                        #Mask the other emission lines
                        if not HaNII:
                            peakwave = coeff3[0]
                            dropwaveline,dropspecline,dropmodelline,dropnoiseline,dropidx,linename,restwave,pdropidx,pdrop = droplines(peakwave=peakwave)
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

                        #Only fit if you're not doing HaNII, otherwise nothing is masked so we don't need to fit again
                        if not HaNII:
                            coeff3, var_matrix3 = curve_fit(gauss3, dropwaveline, dropspecline, p0=guess3, bounds=bounds3)
                        #Compute the values of the fit
                        if not HaNII:
                            gausscurve3 = gauss3(dropwaveline,coeff3[0],coeff3[1])   #
                            amp3 = amp3(dropwaveline,coeff3[0],coeff3[1])   #
                            mu3 = coeff3[0] 
                            stddev3 = np.e**np.abs(coeff3[1])
                            flux3 = amp3[0]
                            scale3 = amp3[1]
                        else:
                            gausscurveHa = gaussHa(dropwaveline,coeffHa[0],coeffHa[1],coeffHa[2],coeffHa[3])
                            ampHa = ampHa(dropwaveline,coeffHa[0],coeffHa[1],coeffHa[2],coeffHa[3])
                            #Fit redshift
                            zgauss = coeffHa[0]
                            #Mean of each line
                            for num in np.arange(0,3):
                                HaNIIdat.at[num,'mu'] = HaNIIdat.iloc[num]['restwave']*(1+zgauss)
                                HaNIIdat.at[num,'sig'] = np.e**np.abs(coeffHa[num+1])
                                HaNIIdat.at[num,'flux'] = ampHa[num]
                                HaNIIdat.at[num,'scale'] = ampHa[3]
                            mu3 = HaNIIdat.iloc[1]['mu']
                            stddev3 = HaNIIdat.iloc[1]['sig']
                            flux3 = HaNIIdat.iloc[1]['flux']
                            scale3 = HaNIIdat.iloc[1]['scale']
                            
                        #Compute chi^2 statistics in the range of the line
                        if not HaNII:
                            #Degrees of freedom: mu, sigma, area, scale
                            dof = 4
                            #Set the lower and upper bounds for the region to find chi2
                            chilb = mu3-2*stddev3
                            chiub = mu3+2*stddev3
                            #Get only the indices in that region
                            cidx = np.logical_and(dropwaveline > chilb-2, dropwaveline < chiub+2)
                            cidx[np.where(dropspecline[cidx]<=0)[0]] = ~cidx[np.where(dropspecline[cidx]<=0)[0]]
                            arrchi2 = divz((dropspecline[cidx]-gausscurve3[cidx]),dropnoiseline[cidx])**2
                            chi2 = np.add.reduce(arrchi2)
                            rchi2 = divz(chi2,len(dropwaveline[cidx])-dof)
                            
                            #Compute the sum of the fluxes in the line in the same region
                            sumflux = 2*np.add.reduce(dropspecline[cidx]-dropmodelline[cidx])
                        else:
                            #Degrees of freedom: z, scale, sigma (x3, for each line), area (x3, for each line)
                            dof = 8
                            cidxarr = []
                            #Set the lower and upper bounds for the region to find chi2
                            for num in np.arange(0,3):
                                HaNIIdat.at[num,'chilb'] = (1+zgauss)*HaNIIdat.iloc[num]['restwave']-2*HaNIIdat.iloc[num]['sig']
                                HaNIIdat.at[num,'chiub'] = (1+zgauss)*HaNIIdat.iloc[num]['restwave']+2*HaNIIdat.iloc[num]['sig']
                                cidxarr.append(np.logical_and(dropwaveline > HaNIIdat.iloc[num]['chilb']-2, dropwaveline < HaNIIdat.iloc[num]['chiub']+2))
                                #Chi2 just in this line
                                arrchi2 = divz((dropspecline[cidxarr[num]]-gausscurveHa[cidxarr[num]]),dropnoiseline[cidxarr[num]])**2
                                HaNIIdat.at[num,'chi2'] = np.add.reduce(arrchi2)
                                HaNIIdat.at[num,'rchi2'] = divz(HaNIIdat.iloc[num]['chi2'],len(dropwaveline[cidxarr[num]])-4)
                                #Compute the sum of the fluxes in the line in the same region
                                HaNIIdat.at[num,'sumflux'] = 2*np.add.reduce(dropspecline[cidxarr[num]]-dropmodelline[cidxarr[num]])
                               
                                zrestline = HaNIIdat.iloc[num]['restwave']*(1+zcc)
                                idx3 = np.logical_and(waveline > zrestline-shrange, waveline < zrestline+shrange)
                                HaNIIdat.at[num,'usig'] = np.sqrt(np.add.reduce(noiseline[idx3]**2))
                            #wsig for each line
                            #Masks out the other two lines, %3 is %3 is mod3
                            for num in np.arange(0,3):
                                wsigidx = np.logical_not(np.logical_or(cidxarr[(num+1)%3],cidxarr[(num+2)%3]))
                                g = np.exp(-0.5*(dropwaveline[wsigidx]-HaNIIdat.iloc[num]['mu'])**2/HaNIIdat.iloc[num]['sig']**2)/np.sqrt(2*np.pi*HaNIIdat.iloc[num]['sig']**2)
                                HaNIIdat.at[num,'wsig'] = np.sqrt(np.sum(g*(dropnoiseline[wsigidx]**2))*np.sqrt(2*np.pi*(HaNIIdat.iloc[num]['sig']**2)))
                            #Chi2 over the whole region
                            cidxtot = np.logical_or(np.logical_or(cidxarr[0],cidxarr[1]),cidxarr[2])
                            arrchi2tot = divz((dropspecline[cidxtot]-gausscurveHa[cidxtot]),dropnoiseline[cidxtot])**2
                            chi2tot = np.add.reduce(arrchi2tot)
                            rchi2tot = divz(chi2tot,len(dropwaveline[cidxtot])-dof)

                            

                            
                        #Now compute the weigthed error
                        #Gaussian curve with area=1
                        if not HaNII:
                            g = np.exp(-0.5*(dropwaveline-mu3)**2/stddev3**2)/np.sqrt(2*np.pi*stddev3**2) #NORMALIZED GAUSSIA
                            wsig = np.sqrt(np.sum(g*(dropnoiseline**2))*np.sqrt(2*np.pi*(stddev3**2)))
                            usig = np.sqrt(np.add.reduce(shnoiseline**2))
                            #Get the string of the nearest wavelength to the line. Used for saving everything
                            linestr = (str(int(np.round(restwave))))
                        else:
                            wsig = HaNIIdat.iloc[1]['wsig']
                            usig = HaNIIdat.iloc[1]['usig']
                            linestr = 'HaNII'
                                        
                        ###Set flags
                        #Make sure the flag isn't 5 (out of view). if it is, don't flag it otherwise
                        if fitflag ==5:
                            pass
                        #Check if more than half of the spectrum is masked - if so, throw it out
                        elif (len(np.where(w<=0)[0])>(len(dropwaveline)/3)):
                            fitflag = 1 #Marks bad data
                        #Check if the width of the line hit the bounds
                        elif (stddev3 > 7.0):
                            fitflag = 2 #Marks bad sigma
                    

                        #Check the flag for each line when fitting HaNII
                        if HaNII:
                            for num in np.arange(0,3):
                                if fitflag == 1: HaNIIdat.at[num,'flag'] = 1
                                elif (HaNIIdat.iloc[num]['sig'] > 7.0):
                                    HaNIIdat.at[num,'flag'] = 2
                                elif ((HaNIIdat.iloc[num]['scale'] < 0.7) or (HaNIIdat.iloc[num]['scale'] > 1.3)):
                                    HaNIIdat.at[num,'flag'] = 4
                                else:
                                    HaNIIdat.at[num,'flag'] = 0
                                    
                                
                        def mkplot(plt10,plt1,plt10b,plt1b,gridsize):
                            #Create the plot    
                            #fig,ax0 = plt.subplots(figsize = (13,7))
                            
                            #Set the axis to the correct number - check if it is flagged or not
                            if fitflag:
                                ax0 = axarrb[plt10b,plt1b]
                                #Increment the counters for next time
                                plt1b = plt1b + 1
                                if plt1b == gridsize:
                                    plt1b = 0
                                    plt10b = plt10b + 1
                            else:
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
                            if pdrop:
                                if dropidx[0].size > 0:
                                    pass
                                    #try: [ax0.axvspan(np.min(waveline[j]),np.max(waveline[j]), color='indianred', alpha=0.1) for j in pdropidx]
                                    #except: [ax0.axvspan(np.min(waveline[j]),np.max(waveline[j]), color='indianred', alpha=0.1) for j in dropidx]
                            #Check if any weights were set to 0 - if so, plot the mask for those
                            if np.where(w<=0)[0].any():
                                [ax0.plot(dropwaveline[j],dropspecline[j], marker='o', color='red', alpha=0.7) for j in np.where(w<=0)[0]]
                            #Plot the region over which we fit chi2
                            if not HaNII:
                                pass
                                #ax0.axvspan(np.min(dropwaveline[cidx]),np.max(dropwaveline[cidx]), color='grey', alpha=0.2, label='chi2 region')
                            else:
                                [ax0.axvspan(np.min(dropwaveline[cidxarr[num]]),np.max(dropwaveline[cidxarr[num]]), color='grey', alpha=0.2, label='chi2 region') for num in np.arange(0,3)]
                            ax0.plot(waveline,modelline*scale3,color='red',label='Model')
                            #ax0.plot(waveline,modelline,color='red',ls='--')
                            #ax0.plot(dropwaveline,guesscurve3,color='orange',label='Initial Guess')
                            ax0.plot(dropwaveline,dropnoiseline,color='orange',label='Noise')
                            ax0.set_ylim(-0.01,np.max(gausscurve3)*1.1)
                            #ax0.plot((zline,zline),(-100,1000),color='black',label='Expected Mean', ls='--')
                            #Titles, axes, legends
                            #ax0.set_title('H$\\beta$, OBJID ' + objs[i][0] + '_' + objs[i][1] + objs[i][2] + ', z=' + str(np.around(zcc,4)),fontsize = titlefont)
                            ax0.legend(fontsize = legendfont,loc=1)
                            ax0.set_xlabel('Wavelength ($\AA$)',fontsize = axisfont)
                            ax0.set_ylabel('Flux ($10^{-17}$ erg/s/${cm}^2/\AA$)',fontsize = axisfont)
                            ax0.tick_params(labelsize = ticksize)
                            return ax0, plt10, plt1, plt10b, plt1b

                        
                        usigrat = divz(flux3,(usig*np.sqrt(rchi2)))

                        ax0,plt10,plt1,plt10b,plt1b = mkplot(plt10,plt1,plt10b,plt1b,gridsize)
                        if not HaNII:
                            ax0.text(0.02,0.95,'z:              ' + str(round(zcc,4)),fontsize = textfont, transform=ax0.transAxes)
                            ax0.text(0.02,0.90,'Mean:       ' + str(round(mu3,2)),fontsize = textfont, transform=ax0.transAxes)      
                            ax0.text(0.02,0.85,'Width:      ' + str(round(stddev3,2)),fontsize = textfont, transform=ax0.transAxes)
                            ax0.text(0.02,0.80,'Scale:       ' + str(round(amp3[1],2)),fontsize = textfont, transform=ax0.transAxes)      
                            ax0.text(0.02,0.75,'Flux:         ' + str(round(amp3[0],2)),fontsize = textfont, transform=ax0.transAxes)
                            #ax0.text(0.02,0.75,'Sumflux:   ' + str(round(sumflux,2)),fontsize = textfont, transform=ax0.transAxes)
                            #ax0.text(0.02,0.70,'Chi2:        ' + str(round(chi2,2)),fontsize = textfont, transform=ax0.transAxes)      
                            #ax0.text(0.02,0.65,'rChi2:       ' + str(round(rchi2,2)),fontsize = textfont, transform=ax0.transAxes)
                            #ax0.text(0.02,0.60,'wsig:        ' + str(round(wsig,3)),fontsize = textfont, transform=ax0.transAxes)
                            #ax0.text(0.02,0.55,'usig:         ' + str(round(usig,3)),fontsize = textfont, transform=ax0.transAxes)
                            #ax0.text(0.02,0.50,'usigrat:     ' + str(round(usigrat,3)),fontsize = textfont, transform=ax0.transAxes)
                        else:
                            ax0.text(0.02,0.95,'Mean:    ' + str(round(HaNIIdat.iloc[0]['mu'],2)),fontsize = textfont, transform=ax0.transAxes)
                            ax0.text(0.24,0.95, str(round(HaNIIdat.iloc[1]['mu'],2)),fontsize = textfont, transform=ax0.transAxes)
                            ax0.text(0.40,0.95, str(round(HaNIIdat.iloc[2]['mu'],2)),fontsize = textfont, transform=ax0.transAxes)
                            ax0.text(0.02,0.90,'Std Dev: ' + str(round(HaNIIdat.iloc[0]['sig'],2)),fontsize = textfont, transform=ax0.transAxes)
                            ax0.text(0.24,0.90, str(round(HaNIIdat.iloc[1]['sig'],2)),fontsize = textfont, transform=ax0.transAxes)
                            ax0.text(0.40,0.90, str(round(HaNIIdat.iloc[2]['sig'],2)),fontsize = textfont, transform=ax0.transAxes)
                            ax0.text(0.02,0.85,'Flux:   ' + str(round(HaNIIdat.iloc[0]['flux'],2)),fontsize = textfont, transform=ax0.transAxes)
                            ax0.text(0.24,0.85, str(round(HaNIIdat.iloc[1]['flux'],2)),fontsize = textfont, transform=ax0.transAxes)
                            ax0.text(0.40,0.85, str(round(HaNIIdat.iloc[2]['flux'],2)),fontsize = textfont, transform=ax0.transAxes)
                            ax0.text(0.02,0.80,'Flag:   ' + str(int(HaNIIdat.iloc[0]['flag'])),fontsize = textfont, transform=ax0.transAxes)
                            ax0.text(0.15,0.80, str(int(HaNIIdat.iloc[1]['flag'])),fontsize = textfont, transform=ax0.transAxes)
                            ax0.text(0.28,0.80, str(int(HaNIIdat.iloc[2]['flag'])),fontsize = textfont, transform=ax0.transAxes)
                            ax0.text(0.02,0.75, 'Scale:  ' + str(round(HaNIIdat.iloc[2]['scale'],2)),fontsize = textfont, transform=ax0.transAxes)
                            ax0.text(0.02,0.70, 'zfit:   ' + str(round(zgauss,4)),fontsize = textfont, transform=ax0.transAxes)
                            ax0.text(0.02,0.65, 'chi2tot: ' + str(round(chi2tot,2)),fontsize = textfont, transform=ax0.transAxes)
                            ax0.text(0.02,0.60, 'rchi2tot: ' + str(round(rchi2tot,2)),fontsize = textfont, transform=ax0.transAxes)
                            
                        
                        if fitflag:
                            pass
                            #ax0.text(0.02,0.50,'flag:          ' + str(fitflag),fontsize = textfont, transform=ax0.transAxes)
                        #fig.text(0.14,0.60,'Redshift:   ' + str(round(zcc,4)),fontsize = textfont)
                        #fig.text(0.14,0.60,'Luminosity (erg/s):   ' + str(round(lumin,2)),fontsize = textfont)
                        if not HaNII:
                            ax0.plot(dropwaveline,gausscurve3,color='black',label='Gaussian fit')
                        else:
                            ax0.plot(dropwaveline,gausscurveHa,color='black',label='Gaussian fit')
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
                        if not HaNII:
                            outarr.at[midx,linestr + '_mean'] = mu3
                            outarr.at[midx,linestr + '_stddev'] = stddev3
                            outarr.at[midx,linestr + '_flux'] = flux3
                            outarr.at[midx,linestr + '_scale'] = scale3
                            outarr.at[midx,linestr + '_chi2'] = chi2
                            outarr.at[midx,linestr + '_rchi2'] = rchi2
                            outarr.at[midx,linestr + '_sumflux'] = sumflux
                            outarr.at[midx,linestr + '_wsig'] = wsig
                            outarr.at[midx,linestr + '_usig'] = usig
                            outarr.at[midx,linestr + '_flag'] = fitflag
                        else:
                            linearr = ['6548','6563','6583']
                            counter = 0
                            for linestr in linearr:
                                outarr.at[midx,linestr + '_mean'] = HaNIIdat.iloc[counter]['mu']
                                outarr.at[midx,linestr + '_stddev'] = HaNIIdat.iloc[counter]['sig']
                                outarr.at[midx,linestr + '_flux'] = HaNIIdat.iloc[counter]['flux']
                                outarr.at[midx,linestr + '_scale'] = HaNIIdat.iloc[counter]['scale']
                                outarr.at[midx,linestr + '_chi2'] = HaNIIdat.iloc[counter]['chi2']
                                outarr.at[midx,linestr + '_rchi2'] = HaNIIdat.iloc[counter]['rchi2']
                                outarr.at[midx,linestr + '_sumflux'] = HaNIIdat.iloc[counter]['sumflux']
                                outarr.at[midx,linestr + '_wsig'] = HaNIIdat.iloc[counter]['wsig']
                                outarr.at[midx,linestr + '_usig'] = HaNIIdat.iloc[counter]['usig']
                                outarr.at[midx,linestr + '_flag'] = HaNIIdat.iloc[counter]['flag']
                                counter = counter + 1
                            outarr.at[midx,'6563_chi2tot'] = chi2tot
                            outarr.at[midx,'6563_rchi2tot'] = rchi2tot
                            outarr.at[midx,'6563_zgauss'] = zgauss
                       
                        '''
                        Flag values:
                        1 - too many zeros, we threw out the fit
                        2 - sigma >8, so it hit the bounds. 
                        4 - scale >1.3 or <0.7, probably something wrong with spectrum in the region
                        5 - the line is not redshifted enough to be in view (e.g. 3727 OII)
                        '''
            
                        
                    except (RuntimeError):
                        ax0.text(0.14,0.84,'Fitting Failed',fontsize = textfont, transform=ax0.transAxes)
                        #plt.show()
                else: print('Bad data at ' + str(line) + ', too many NaN. ' + 'flx_' + objs[i][0] + '_feb1' + objs[i][2] + '_' + objs[i][1] + 'big.fits' )
               

    #If not, give an error but continue
    else: print('Could not read file ' + flxfits)

###Editing the datafile
#Sort by OBJID
outarr = outarr.sort_values('OBJID')
#Sort the columns so the lines are next to each other
outarr = outarr.reindex(sorted(outarr.columns), axis=1)
#Remove all NaN and replace them with -99
outarr = outarr.fillna(value = -99.999999999999)
#Remove columns with this, then take it back out
#outarr = outarr.drop('Ha_chi2',axis=1)

#Write the file
#outarr.to_csv(dataout,index=False)


#Save the figure
#plt.show()
fig.tight_layout()
figb.tight_layout()
if HaNII: linename = 'HaNII'
fig.savefig(figout + str(int(np.round(restwave))) + '_' + linename + '_' + letnum + '.pdf')
figb.savefig(figout + str(int(np.round(restwave))) + '_' + linename + '_' + letnum + '_flagged.pdf')
plt.close(fig)
plt.close(figb)
                        


'''
Make a bpt diagram, look at spectra of possible AGN
'''
