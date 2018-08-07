#Creates a BPT diagram for all objects, and a second figure that shows objects for which single lines are low
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

#Location of the equivalent width data
ewdata = '/Users/blorenz/COSMOS/COSMOSData/lineew.txt'
#Read in the ew of the lines 
ew_df = ascii.read(ewdata).to_pandas()

#The location to store the scale and its stddev of each line
qualdatapath = '/Users/blorenz/COSMOS/COSMOSData/dataqual.txt'
#Read in the scale of the lines 
dataqual = ascii.read(qualdatapath).to_pandas()

#File with the error array
errdatapath = '/Users/blorenz/COSMOS/COSMOSData/errs.txt'
#Read in the scale of the lines 
err_df = ascii.read(errdatapath,data_start=1,header_start=0,format='csv').to_pandas()


#Read the datafile:
fluxdata = ascii.read(fluxdatapath).to_pandas()

#Division function
def divz(X,Y):
        return X/np.where(Y,Y,Y+1)*np.not_equal(Y,0)

#The location of the muzzin et al data:
mdatapath = '/Users/blorenz/COSMOS/muzzin_data/UVISTA_final_colors_sfrs_v4.1.dat'
#Read in the muzzin data
mdata = ascii.read(mdatapath).to_pandas()
mdata = mdata.rename(columns={'ID':'OBJID'})

fluxdata = pd.merge(fluxdata,mdata)

#Location of the reddening data
reddata = '/Users/blorenz/COSMOS/COSMOSData/reddenings.txt'
#Read in the ew of the lines 
red_df = ascii.read(reddata).to_pandas()
fluxdata = pd.merge(fluxdata,red_df,on='fluxfile')

fluxdata = fluxdata[fluxdata.av>0]

av_med_df = pd.DataFrame()

av_med_df.at[0,'mr1'] = np.median(fluxdata[fluxdata.LMASS<9.25].av)
av_med_df.at[0,'mr2'] = np.median(fluxdata[np.logical_and(fluxdata.LMASS>9.25,fluxdata.LMASS<9.5)].av)
av_med_df.at[0,'mr3'] = np.median(fluxdata[np.logical_and(fluxdata.LMASS>9.5,fluxdata.LMASS<9.75)].av)
av_med_df.at[0,'mr4'] = np.median(fluxdata[fluxdata.LMASS>9.75].av)

av_med_df.to_csv('/Users/blorenz/COSMOS/COSMOSData/av_med_df.txt',index=False)

