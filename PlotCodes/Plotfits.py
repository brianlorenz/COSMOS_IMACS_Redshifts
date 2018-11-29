#Plot a .fits file

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii, fits
import sys, os, string
import pandas as pd


fitsfile = sys.argv[1]

data = fits.open(fitsfile)[0].data
head = fits.open(fitsfile)[0].header

d0 = data[0]
d1 = data[1]
d2 = data[2]
d3 = data[3]
d4 = data[4]
#d5 = data[5]
#d6 = data[6]
#d7 = data[7]
#d8 = data[8]


crval1 = head["crval1"]
crpix1 = head["crpix1"]
cdelt1 = head["cdelt1"]
naxis1 = head["naxis1"]
dcflag = head["dc-flag"]
exptime = head['exptime']
wavelength = (1.0+np.arange(naxis1)-crpix1)*cdelt1 + crval1


fig,axarr = plt.subplots(figsize=(13,7))
#ax0,ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8 = axarr[0,0],axarr[0,1],axarr[0,2],axarr[1,0],axarr[1,1],axarr[1,2],axarr[2,0],axarr[2,1],axarr[2,2]
axarr.plot(wavelength,data[0])
#ax1.plot(wavelength,data[1])
#ax2.plot(wavelength,data[2])
#ax3.plot(wavelength,data[3])
#ax4.plot(wavelength,data[4])
#ax5.plot(wavelength,data[5])
#ax6.plot(wavelength,data[6])
#ax7.plot(wavelength,data[7])
#ax8.plot(wavelength,data[8])
plt.show()
