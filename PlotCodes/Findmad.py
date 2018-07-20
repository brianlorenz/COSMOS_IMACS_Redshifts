import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import sys, os, string
import pandas as pd
from astropy.io import fits
import collections

#The location with the file for all of our data
fluxdatapath = '/Users/blorenz/COSMOS/COSMOSData/lineflux.txt'
#The location to store the MAD of each line
madout = '/Users/blorenz/COSMOS/COSMOSData/linemad.txt'

#Read the datafile
fluxdata = ascii.read(fluxdatapath).to_pandas()

#Division function
def divz(X,Y):
        return X/np.where(Y,Y,Y+1)*np.not_equal(Y,0)

#Finds the OBJIDS of all duplicates
dupobjids = [item for item, count in collections.Counter(fluxdata.OBJID).items() if count > 1]
#Pulls out the errors of the duplications
duprows = [fluxdata[(fluxdata.OBJID == i)] for i in dupobjids]

#Get the strings of each line
lines = ['4861','4959','5007','6563','4340','4102','6548','6583','3727','6563_fix','6548_fix','6583_fix']
#Set up a df to store them
mad_df = pd.DataFrame()

#Loop over every line
for line in lines:
    #Compute the difference for all duplicates in the line
    diff = [np.abs(divz(i.iloc[0][line + '_flux'],i.iloc[0][line + '_scale'])-divz(i.iloc[1][line + '_flux'],i.iloc[1][line + '_scale'])) for i in duprows if ((i.iloc[0][line + '_flag']==0) and (i.iloc[1][line + '_flag']==0))]
    #Compute 1.49*mad
    mad = 1.49*np.median(diff)
    mad = mad/np.sqrt(2)
    #Store the result to the df
    mad_df.at[0,line + '_mad'] = mad


#Use this to set a new mad, then take it out
#mad_df.at[0,'6548_mad'] = 0.5

#Sort the df by linename
mad_df = mad_df.reindex(sorted(mad_df.columns), axis=1)

#Save the df
mad_df.to_csv(madout,index=False)
