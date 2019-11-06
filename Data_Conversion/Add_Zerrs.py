# Adds the missing redshift uncertainties column to our data

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import sys
import os
import string
import pandas as pd

# The file for all of our data
fulldataout = '/Users/galaxies-air/COSMOS/COSMOSData/all_c_hasinger.txt'

# Read in all of our data
ourdata = ascii.read(fulldataout).to_pandas()

# Search each row for the temp number, find the redshift associated with that number and append it to a list
template_int = ourdata.temp.astype(int)


dzhis = []
dzlos = []
for i in range(0, len(ourdata)):
    # Make sure the template is nonzero for the strange cases where it is
    current_temp = template_int.iloc[i]
    if current_temp > 0:
        dzhis.append(ourdata[f'dzhi{current_temp}'].iloc[i])
        dzlos.append(ourdata[f'dzlo{current_temp}'].iloc[i])
    else:
        dzhis.append(None)
        dzlos.append(None)

        #
OBJID = []
for i in range(0, len(ourdata)):
    OBJID.append(ourdata.iloc[i].ImageName[4:10])

ourdata['dzhi'] = dzhis
ourdata['dzlo'] = dzlos

ourdata.to_csv(fulldataout.replace(
    'all_c_hasinger', 'all_c_hasinger_dz'), index=False)
