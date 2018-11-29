#Creates a UVJ diagram split by good, low, and bad measurements for specified lines
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import sys, os, string
import pandas as pd
from astropy.io import fits
import collections
from astropy.cosmology import WMAP9 as cosmo
from astropy.stats import biweight_midvariance

#Folder to save the figures
figout = '/Users/blorenz/COSMOS/Reports/2018/Images/'

#The location with the file for all of our data
fluxdatapath = '/Users/blorenz/COSMOS/COSMOSData/lineflux_red.txt'

#The location with the file for all of our data
sfout = '/Users/blorenz/COSMOS/COSMOSData/sfrs.txt'

#Location of the equivalent width data
ewdata = '/Users/blorenz/COSMOS/COSMOSData/lineew.txt'
#Read in the ew of the lines 
ew_df = ascii.read(ewdata).to_pandas()

#The location to store the scale and its stddev of each line
qualdatapath = '/Users/blorenz/COSMOS/COSMOSData/dataqual.txt'
#Read in the scale of the lines 
dataqual = ascii.read(qualdatapath).to_pandas()
d = {'True': True, 'False': False}

#File with the error array
errdatapath = '/Users/blorenz/COSMOS/COSMOSData/errs.txt'
#Read in the scale of the lines 
err_df = ascii.read(errdatapath,data_start=1,header_start=0,format='csv').to_pandas()


#File with the error array
errreddatapath = '/Users/blorenz/COSMOS/COSMOSData/errs_red.txt'
#Read in the scale of the lines 
err_dfred = ascii.read(errreddatapath,data_start=1,header_start=0,format='csv').to_pandas()

#Read the datafile:
fluxdata = ascii.read(fluxdatapath).to_pandas()

#File with the structural properties
spropdatapath = '/Users/blorenz/COSMOS/COSMOSData/struct_prop.txt'
#Read in the scale of the lines 
sprop_df = ascii.read(spropdatapath).to_pandas()
sprop_df = sprop_df.rename(columns={'id':'OBJID'})
fluxdata = pd.merge(fluxdata,sprop_df)

#The location with the file for the filter data
filtdatapath = '/Users/blorenz/COSMOS/COSMOSData/all_c_hasinger.txt'
#Read in the data
filtdata = ascii.read(filtdatapath).to_pandas()
cordata = filtdata[['id','Ks','eKs','Ks_tot','eKs_tot']]
cordata = cordata.rename(columns={'id':'OBJID'})

fluxdata = pd.merge(fluxdata,cordata,on='OBJID',how='inner')
fluxdata = fluxdata.drop_duplicates()
fluxdata = fluxdata.reset_index()

#Fontsizes for plotting
axisfont = 24
ticksize = 18
ticks = 8
titlefont = 24
legendfont = 16
textfont = 16

#Division function
def divz(X,Y):
        return X/np.where(Y,Y,Y+1)*np.not_equal(Y,0)


lines = ['6563_fix']


sfr_df = pd.DataFrame()
sfr_df['fluxfile'] = fluxdata['fluxfile']

sfr_df['sfr_err_u'] = err_dfred['6563_fix_u']*1e-17*4*np.pi*((cosmo.luminosity_distance(fluxdata['zcc'])*3.086*1e24)**2)*10**(-41.27)*(fluxdata.Ks_tot/fluxdata.Ks)
sfr_df['ssfr_err_u'] = divz(sfr_df['sfr_err_u'],10**fluxdata.LMASS)

sfr_df['sfr_err_d'] = err_dfred['6563_fix_d']*1e-17*4*np.pi*((cosmo.luminosity_distance(fluxdata['zcc'])*3.086*1e24)**2)*10**(-41.27)*(fluxdata.Ks_tot/fluxdata.Ks)
sfr_df['ssfr_err_d'] = divz(sfr_df['sfr_err_d'],10**fluxdata.LMASS)

sfr_df['SFR'] = fluxdata['6563_fix_flux_red']*1e-17*4*np.pi*((cosmo.luminosity_distance(fluxdata['zcc'])*3.086*1e24)**2)*(fluxdata.Ks_tot/fluxdata.Ks)*10**(-41.27) #Kennicutt and Evans (2012)
 
sfr_df['SFR'] = sfr_df['SFR'].astype(float)
sfr_df['sSFR'] = divz(sfr_df['SFR'],10**fluxdata['LMASS'])


sfr_df.to_csv(sfout,index=False)
