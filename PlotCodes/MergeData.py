#Merges all of the separate data files into one

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import sys, os, string
import pandas as pd

#Set where our measured data is stored
#ourdatapath - the folder where the text files are stored
#ourdatafiles - a list of the names of all of the data text files
ourdatapath = '/Users/blorenz/COSMOS/COSMOSData/corFitsFileOut/crossCorOut/'
ourdatafiles = ['verb_a6.txt','verb_b6.txt','verb_b7.txt','verb_d6.txt','verb_e6.txt','verb_e7.txt','verb_f6.txt','verb_g6.txt','verb_h6.txt','verb_i6.txt','verb_j7.txt']

#The output file for all of our data
fulldataout = '/Users/blorenz/COSMOS/COSMOSData/all_data.txt'

#Read in all of our data
firsttime = 1
for i in ourdatafiles:
    tempdata = ascii.read(ourdatapath+i).to_pandas()
    if firsttime == 1:
        ourdata = tempdata
        firsttime = 0
    else:
        ourdata = ourdata.append(tempdata)

#Write to a csv
ourdata.to_csv(fulldataout,index=False)

#Merge our data with the UVISTA catalog
#The location of the muzzin et al data:
mdatapath = '/Users/blorenz/COSMOS/muzzin_data/UVISTA_final_v4.1.cat.txt'

#Read in all of our data
ourdata = ascii.read(fulldataout).to_pandas()

ourdata.rename(columns={'OBJID':'id'}, inplace=True)

#Read in the dataset from Muzzin et al
mdata = ascii.read(mdatapath).to_pandas()
ourdata = pd.merge(ourdata, mdata, on='id', how='inner')
ourdata.to_csv(fulldataout.replace('data','data_UVISTA'),index=False)


