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

#Search each row for the temp number, find the redshift associated with that number and append it to a list
z_cc = []
for i in range(0,len(ourdata)):
    temp =  ourdata.iloc[i].temp
    if temp == 23:
        z_cc.append(ourdata.iloc[i,1])
    elif temp == 24:
        z_cc.append(ourdata.iloc[i,7])
    elif temp == 25:
        z_cc.append(ourdata.iloc[i,13])
    elif temp == 26:
        z_cc.append(ourdata.iloc[i,19])
    elif temp == 27:
        z_cc.append(ourdata.iloc[i,25])
    else:
        z_cc.append(None)

        #
OBJID = []
for i in range(0,len(ourdata)):
    OBJID.append(ourdata.iloc[i].ImageName[4:10])

#Turn the list into a df and join it with the data
z_ccdf = pd.DataFrame({'z_cc':z_cc})
ourdata = ourdata.join(z_ccdf)

OBJIDdf = pd.DataFrame({'OBJID':OBJID}).astype(str)
ourdata = ourdata.join(OBJIDdf)

ourdata.to_csv(fulldataout.replace('data_UVISTA','data_UVISTA_z'),index=False)

