import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import sys, os, string
import pandas as pd
from astropy.io import fits
import collections
from astropy.stats import biweight_midvariance

#The location with the file for all of our data
fluxdatapath = '/Users/blorenz/COSMOS/COSMOSData/lineflux.txt'
#The location to store the MAD of each line
madout = '/Users/blorenz/COSMOS/COSMOSData/linemad.txt'
#The location to store the scale and its stddev of each line
scaleout = '/Users/blorenz/COSMOS/COSMOSData/scales.txt'
#Folder to save the figures
figout = '/Users/blorenz/COSMOS/Reports/2018/Images/'


#Read the datafile
fluxdata = ascii.read(fluxdatapath).to_pandas()

#Get the strings of each line
lines = ['4861','4959','5007','6563','4340','4102','6548','6583','3727','6563_fix','6548_fix','6583_fix']
lines = np.sort(lines)
#Set up a df to store them
mad_df = pd.DataFrame()
scale_df = pd.DataFrame()



#Fontsizes for plotting
axisfont = 18
ticksize = 16
titlefont = 24
legendfont = 16
textfont = 16

fig,axarr5 = plt.subplots(3,4,figsize=(30,20))
axarr5 = np.reshape(axarr5,12)
counter=0
for line in lines:
    #Compute the scale dataframe and make a plot to display it
    ax = axarr5[counter]
    medscale = np.median(np.log10(fluxdata[fluxdata[line+'_scale']>0][line+'_scale']))
    sigscale = np.sqrt(biweight_midvariance(np.log10(fluxdata[fluxdata[line+'_scale']>0][line+'_scale'])))
    scale_df.at[0,line + '_medscale'] = medscale
    scale_df.at[0,line + '_sigscale'] = sigscale
    #Make the plot
    bins = np.log10(np.arange(0.05,4,0.05))
    ax.hist(np.log10(fluxdata[fluxdata[line+'_scale']>0][line+'_scale']),bins=bins,color='grey',label=None)
    #Plot the median, 1sig, and 3 sig stddevs
    ax.plot((medscale,medscale),(-100,1000),color='red',ls='-',label='Median')
    ax.plot((medscale-sigscale,medscale-sigscale),(-100,1000),color='pink',ls='-',label='1 sigma')
    ax.plot((medscale+sigscale,medscale+sigscale),(-100,1000),color='pink',ls='-')
    ax.plot((medscale-3*sigscale,medscale-3*sigscale),(-100,1000),color='thistle',ls='-',label='1 sigma')
    ax.plot((medscale+3*sigscale,medscale+3*sigscale),(-100,1000),color='thistle',ls='-')
    #Titles, axes, legends
    #ax.set_xscale('log')
    ax.set_title(line + ' Scale Histogram',fontsize = titlefont)
    ax.set_xlabel(line + ' Scale',fontsize = axisfont)
    ax.set_ylabel('Counts',fontsize = axisfont)
    ax.set_xlim(np.min(bins),np.max(bins))
    ax.set_ylim(0,500)
    ax.tick_params(labelsize = ticksize)
    ax.legend(fontsize=axisfont)
    counter = counter+1
fig.tight_layout()
fig.savefig(figout + 'scale_hist.pdf')
plt.close(fig)


#Division function
def divz(X,Y):
        return X/np.where(Y,Y,Y+1)*np.not_equal(Y,0)

#Finds the OBJIDS of all duplicates
dupobjids = [item for item, count in collections.Counter(fluxdata.OBJID).items() if count > 1]
#Pulls out the errors of the duplications
duprows = [fluxdata[(fluxdata.OBJID == i)] for i in dupobjids]
    
#Loop over every line
for line in lines:
    #Compute the difference for all duplicates in the line
    diff = [np.abs(divz(i.iloc[0][line + '_flux'],i.iloc[0][line + '_scale'])-divz(i.iloc[1][line + '_flux'],i.iloc[1][line + '_scale'])) for i in duprows if (((i.iloc[0][line + '_flag'] in [0,4]) and (i.iloc[1][line + '_flag'] in [0,4])) and (np.abs(i.iloc[0][line+'_scale']-scale_df[line+'_medscale'][0]) <  (3*scale_df[line+'_sigscale'][0])) and (np.abs(i.iloc[1][line+'_scale']-scale_df[line+'_medscale'][0]) <  (3*scale_df[line+'_sigscale'][0])))]
    #Compute 1.49*mad
    mad = 1.49*np.median(diff)
    mad = mad/np.sqrt(2)
    #Store the result to the df
    mad_df.at[0,line + '_mad'] = mad


#Use this to set a new mad, then take it out
#mad_df.at[0,'6548_mad'] = 0.5

#Sort the df by linename
mad_df = mad_df.reindex(sorted(mad_df.columns), axis=1)
scale_df = scale_df.reindex(sorted(scale_df.columns), axis=1)

#Save the df
mad_df.to_csv(madout,index=False)
scale_df.to_csv(scaleout,index=False)
