import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import sys, os, string
import pandas as pd
from astropy.io import fits
import collections

#Folder to save the figures
figout = '/Users/blorenz/COSMOS/COSMOSData/fluxImages/'

#The location with the file for all of our data
fluxdatapath = '/Users/blorenz/COSMOS/COSMOSData/lineflux.txt'
#The location to store the MAD of each line
madout = '/Users/blorenz/COSMOS/COSMOSData/linemad.txt'

#Read the datafile (if there is one), then create a blank one to write to:
fluxdata = ascii.read(fluxdatapath).to_pandas()

#Fontsizes for plotting
axisfont = 18
ticksize = 16
titlefont = 24
legendfont = 16
textfont = 16

#Division function
def divz(X,Y):
        return X/np.where(Y,Y,Y+1)*np.not_equal(Y,0)

#Finds the OBJIDS of all duplicates
dupobjids = [item for item, count in collections.Counter(fluxdata.OBJID).items() if count > 1]
#Pulls out the errors of the duplications
duprows = [fluxdata[(fluxdata.OBJID == i)] for i in dupobjids]

#Get the strings of each line
lines = ['4861','4959','5007','6548','6563','6583']
#Set up a df to store them
mad_df = pd.DataFrame()

#Loop over every line
for line in lines:
    #Compute the difference for all duplicates in the line
    diff = [np.abs(divz(i.iloc[0][line + '_flux'],i.iloc[0][line + '_scale'])-divz(i.iloc[1][line + '_flux'],i.iloc[1][line + '_scale'])) for i in duprows]
    #Compute 1.49*mad
    mad = 1.49*np.median(diff)
    #Store the result to the df
    mad_df.at[0,line + '_mad'] = mad

#Save the df
mad_df.to_csv(madout,index=False)
'''
HBdiff = [np.abs(i.iloc[0].HB_flux-i.iloc[1].HB_flux) for i in duprows]
HB_std = 1.49*np.median(HBdiff)
d4959diff = [np.abs(i.iloc[0]['4959_flux']-i.iloc[1]['4959_flux']) for i in duprows]
d4959_std = 1.49*np.median(d4959diff)
d5007diff = [np.abs(i.iloc[0]['5007_flux']-i.iloc[1]['5007_flux']) for i in duprows]
d5007_std = 1.49*np.median(d5007diff)
d6548diff = [np.abs(i.iloc[0]['6548_flux']-i.iloc[1]['6548_flux']) for i in duprows]
d6548_std = 1.49*np.median(d6548diff)
d6583diff = [np.abs(i.iloc[0]['6583_flux']-i.iloc[1]['6583_flux']) for i in duprows]
d6583_std = 1.49*np.median(d6583diff)

std_df = pd.DataFrame()
std_df.at[0,'Ha_std'] = Ha_std
std_df.at[0,'HB_std'] = HB_std
std_df.at[0,'4959_std'] = d4959_std
std_df.at[0,'5007_std'] = d5007_std
std_df.at[0,'6548_std'] = d6548_std
std_df.at[0,'6583_std'] = d6583_std
'''

sys.exit()

fig,ax = plt.subplots(figsize=(13,7))
for i in duprows:
    ax.errorbar(i.iloc[0].Ha_flux,i.iloc[1].Ha_flux,i.iloc[0].Ha_wsig,i.iloc[1].Ha_wsig,color='cornflowerblue',marker='o')


#Titles, axes, legends
ax.set_title('Uncertainty Estimates in Ha for Repeated Objects',fontsize = titlefont)
ax.legend(fontsize = legendfont,loc=1)
ax.set_xlabel('First Ha',fontsize = axisfont)
ax.set_ylabel('Second Ha',fontsize = axisfont)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(0,120)
ax.set_ylim(0,120)
ax.tick_params(labelsize = ticksize)
plt.show()
