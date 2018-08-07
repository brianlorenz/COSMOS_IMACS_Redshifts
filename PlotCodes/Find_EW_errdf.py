import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import sys, os, string
import pandas as pd
from astropy.io import fits
import collections

#Folder to save the figures
figout = '/Users/blorenz/COSMOS/Reports/2018/Images/'

#The location with the file for all of our data
fluxdatapath = '/Users/blorenz/COSMOS/COSMOSData/lineflux.txt'

#File for the MAD of the difference in flux of duplicates in each line
maddatapath = '/Users/blorenz/COSMOS/COSMOSData/linemad.txt'

#File to write the ew array to
dataout = '/Users/blorenz/COSMOS/COSMOSData/lineew.txt'

#File to write the error array to
errout = '/Users/blorenz/COSMOS/COSMOSData/errs.txt'

#Merge our data with the UVISTA catalog
#The location of the muzzin et al data:
mdatapath = '/Users/blorenz/COSMOS/muzzin_data/UVISTA_final_colors_sfrs_v4.1.dat'

#Where the calibrated spectra are stored
caldatapath = '/Users/blorenz/COSMOS/COSMOSData/flxFitsFileOut/'

#Read in the mad of the lines 
mad_df = ascii.read(maddatapath).to_pandas()


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

        
lines = ['4861','5007','4959','6563_fix','6583_fix','6548_fix','3727','4340','4102','6563','6583','6548']
lines = np.sort(lines)
ew_df = pd.DataFrame()
for i in range(0,len(fluxdata)):
        flxfits = caldatapath+fluxdata.iloc[i]['fluxfile']
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
        for l in range(0,len(lines)):
                if len(lines[l]) == 8: linename = lines[l][0:4]
                else: linename = lines[l]
                lineint = int(linename)
                linez = lineint*(1+fluxdata.iloc[i]['zcc'])
                #Get the flux and scale
                flux = fluxdata.iloc[i][lines[l]+'_flux']
                scale = fluxdata.iloc[i][lines[l]+'_scale']
                stddev = fluxdata.iloc[i][lines[l]+'_stddev']
                mean = fluxdata.iloc[i][lines[l]+'_mean']
                nonzeroidx = (spec != 0)
                #Compute Ews
                #Set the lower and upper bounds for the line
                lb = mean-2*stddev
                ub = mean+2*stddev
                #Get only the indices in that region
                idx = np.logical_and(wavelength > lb-2, wavelength < ub+2)
                idx = np.logical_and(idx,nonzeroidx)
                if len(spec[idx] > 4):
                        ratio = divz(spec[idx],model[idx])
                        lineew = np.abs(np.sum(1-ratio))
                else: lineew = -99.999999999
                #Now get the model ew
                #Left side of view window for line
                lineidxleft = np.logical_and(wavelength>(linez-50),wavelength<(linez-20))
                #Right side
                lineidxright = np.logical_and(wavelength<(linez+50),wavelength>(linez+20))
                #All indices except the line
                lineidx = np.logical_or(lineidxleft,lineidxright)
                lineidx = np.logical_and(lineidx,nonzeroidx)
                #Model continuum
                modelcont = np.median(model[lineidx])
                modelcont = np.nan_to_num(modelcont)
                if (modelcont != 0):
                        ratiom = divz(model[idx],modelcont)
                        modelew = np.abs(np.sum(1-ratiom))
                        modelabs = np.sum(model[idx]-modelcont)*2
                else:
                        modelew = -99.999999999
                        modelabs = -99.999999999
                #Save the computed EWs
                ew_df.at[i,lines[l]+'_EW'] = lineew
                ew_df.at[i,lines[l]+'_modelEW'] = modelew
                ew_df.at[i,lines[l]+'_modelabs'] = np.abs(modelabs)
                
        ew_df.at[i,'fluxfile'] = fluxdata.iloc[i]['fluxfile']
ew_df.to_csv(dataout,index=False)

#Also, compute the errors for each of the lines and save those
err_df=pd.DataFrame()
err_df['fluxfile'] = fluxdata['fluxfile']
for line in lines:
    err = np.sqrt(fluxdata[line+'_usig']**2+(0.25*ew_df[line+'_modelabs'])**2)
    err_df[line] = err
err_df.to_csv(errout,index=False)
