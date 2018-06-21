#Adds the missing redshift column to our data

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import sys, os, string
import pandas as pd

#The file for all of our data
fulldataout = '/Users/blorenz/COSMOS/COSMOSData/all_data_UVISTA.txt'

#Read in all of our data
ourdata = ascii.read(fulldataout).to_pandas()

z_cc = []
for i in range(0,len(ourdata)):
    temp =  ourdata.iloc[i].temp
    if temp == 23:
        z_cc.append(ourdata.iloc[i,40])
    elif temp == 24:
        z_cc.append(ourdata.iloc[i,41])
    elif temp == 25:
        z_cc.append(ourdata.iloc[i,42])
    elif temp == 26:
        z_cc.append(ourdata.iloc[i,43])
    elif temp == 27:
        z_cc.append(ourdata.iloc[i,44])
    else:
        z_cc.append(None)


z_ccdf = pd.DataFrame({'z_cc':z_cc})
ourdata = ourdata.join(z_ccdf)

ourdata.to_csv(fulldataout.replace('data_UVISTA','data_UVISTA_z'),index=False)

