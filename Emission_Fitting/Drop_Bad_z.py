#Gets ride of objects outside of our redshift range
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
fluxdatapath = '/Users/blorenz/COSMOS/COSMOSData/lineflux_tot.txt'

#Output file
dataout = '/Users/blorenz/COSMOS/COSMOSData/lineflux.txt'

#Read the datafile:
fluxdata = ascii.read(fluxdatapath).to_pandas()

#Select on the rows with redshifts in the correct range
fluxdata = fluxdata[np.logical_not(np.logical_or((fluxdata.zcc<0.3),(fluxdata.zcc>0.4)))]

#Save the datafile without those redshifts
fluxdata.to_csv(dataout,index=False)
