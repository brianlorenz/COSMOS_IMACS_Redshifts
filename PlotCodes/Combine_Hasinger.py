#Combines our dataset with Hasinger et al (2018) using the nearest object
#Requires the dataset from https://arxiv.org/abs/1803.09251

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import sys, os, string
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy import units as u


#The location with the file for all of our data
ourdatapath = '/Users/blorenz/COSMOS/COSMOSData/all_data_UVISTA_z.txt'
#The location of the Hasigner data
hdatapath = '/Users/blorenz/COSMOS/deimos_10K_March2018/deimos_redshifts.txt'
#The location to output the comnbined data
outpath = '/Users/blorenz/COSMOS/COSMOSData/all_c_hasinger.txt'

#Read in all of our data
ourdata = ascii.read(ourdatapath)

#Read in the dataset from Hasinger et al
hdata = ascii.read(hdatapath)

#Tutorial here http://www.astropy.org/astropy-tutorials/Coordinates.html
coo_ourdata = SkyCoord(ourdata['ra']*u.deg, ourdata['dec']*u.deg)
coo_hdata = SkyCoord(hdata['Ra']*u.deg, hdata['Dec']*u.deg)
idx_match, d2d_match, d3d_match = coo_ourdata.match_to_catalog_sky(coo_hdata)

#Plot the separation histogram
fig,ax = plt.subplots()
ax.hist(d2d_match.arcsec,bins=60)
fig.savefig('sep_hist.pdf')
plt.close(fig)

#Convert the separation to a pandas file
d2d_match_df = pd.DataFrame({'d2d':d2d_match})
d3d_match_df = pd.DataFrame({'d3d':d3d_match})

#Attach the dataframes
joindata = ourdata.to_pandas().join(hdata[idx_match].to_pandas())
#Add on the new separation information
joindata = joindata.join(d2d_match_df)
joindata = joindata.join(d3d_match_df)

#Write the output
joindata.to_csv(outpath,index=False)

