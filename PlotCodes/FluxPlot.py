#Computes a plot F_lambda vs log(lambda) for all of our objects, saves them to the folder figout

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import sys, os, string
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy import units as u

#Folder to save the figures
figout = '/Users/blorenz/COSMOS/COSMOSData/FluxPlot_Images/'

#The location with the file for all of our data
ourdatapath = '/Users/blorenz/COSMOS/COSMOSData/all_c_hasinger.txt'

#Read in all of our data
ourdata = ascii.read(ourdatapath).to_pandas()

#Read in the list of filters and wavelengths
filterdatapath = '/Users/blorenz/COSMOS/COSMOS_IMACS_Redshifts/PlotCodes/filterlist.txt'
filterlist = ascii.read(filterdatapath).to_pandas()

#Select the columns with data in order of increasing wavelength
ourfilters = ourdata[['u','IB427','B','IB464','gp','IA484','IB505','IA527','V','IB574','IA624','rp','IA679','IB709','IA738','IA767','ip','IB827','zp','Y','J','H','Ks']]

#Select the error columns
ourerrors = ourdata[['eu','eIB427','eB','eIB464','egp','eIA484','eIB505','eIA527','eV','eIB574','eIA624','erp','eIA679','eIB709','eIA738','eIA767','eip','eIB827','ezp','eY','eJ','eH','eKs']]


for i in range(0,len(ourfilters)):
    OBJID = ourdata.iloc[i].ImageName[4:10]
    fig,ax = plt.subplots()
    ax.errorbar(filterlist.Wavelength,np.divide(ourfilters.iloc[i],filterlist.Wavelength),yerr = np.divide(ourerrors.iloc[i],filterlist.Wavelength),
                color='cornflowerblue',ls='none',marker='.')
    ax.set_xscale('log')
    ax.set_xlabel('Wavelength (Angstroms)')
    ax.set_ylabel('F_lambda (units still off)')
    ax.set_title('ULTRAVista Photometry for our objects, objid ' + OBJID)
    fig.savefig(figout + OBJID + '_FluxPlot.png')
    plt.close()

