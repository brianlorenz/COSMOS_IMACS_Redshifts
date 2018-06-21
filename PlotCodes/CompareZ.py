#Compares our measured redshifts to the measurements made by Hasinger et al (2018)
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
ourdatapath = '/Users/blorenz/COSMOS/COSMOSData/all_c_hasinger.txt'

#Read in all of our data
ourdata = ascii.read(ourdatapath).to_pandas()
print len(ourdata)

#Remove all of the galaxies where we are unsure
ourdata = ourdata[ourdata.Unsure==0]
ourdata = ourdata[ourdata.Star==0]
print len(ourdata)

#Separate into three distance bins
ourdata_close = ourdata[ourdata.d2d < .000556]
ourdata_med = ourdata[(ourdata.d2d < .00139) & (ourdata.d2d > .000556)]
ourdata_far = ourdata[ourdata.d2d >= .00139]


#Find the differenc ein our measurements compared to hasinger's
zdiff_close = ourdata_close.z_cc - ourdata_close.zspec
zdiff_med = ourdata_med.z_cc - ourdata_med.zspec
zdiff_far = ourdata_far.z_cc - ourdata_far.zspec

print(len(zdiff_close),len(zdiff_med),len(zdiff_far))
#Remove rows with missing values
zdiff_close = zdiff_close.dropna()
zdiff_med = zdiff_med.dropna()
zdiff_far = zdiff_far.dropna()

#Set the number of histogram bins
bins=30

#Plot the difference
fig,ax = plt.subplots()
#ax.hist([zdiff_close,zdiff_med,zdiff_far],bins=bins,label=['2 arcsec (' + str(len(zdiff_close)) + ')','5 arcsec (' + str(len(zdiff_med)) + ')','Far (' + str(len(zdiff_far)) + ')'])
ax.scatter(ourdata_close.z_cc,ourdata_close.zspec)
plt.legend()
fig.savefig('z_compare_2.pdf')
plt.close(fig)

