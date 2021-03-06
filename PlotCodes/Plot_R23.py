#Creates an R23 diagram - see Kewley and Ellison (2008)

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import sys, os, string
import pandas as pd
from astropy.io import fits
import collections

#Folder to save the figures
figout = '/Users/blorenz/COSMOS/Reports/2018/Images/'

#The location with the file for all of our data
fluxdatapath = '/Users/blorenz/COSMOS/COSMOSData/lineflux.txt'

#Location of the equivalent width data
ewdata = '/Users/blorenz/COSMOS/COSMOSData/lineew.txt'
#Read in the ew of the lines 
ew_df = ascii.read(ewdata).to_pandas()

#The location to store the scale and its stddev of each line
qualdatapath = '/Users/blorenz/COSMOS/COSMOSData/dataqual.txt'
#Read in the scale of the lines 
dataqual = ascii.read(qualdatapath).to_pandas()
d = {'True': True, 'False': False}

#File with the error array
errdatapath = '/Users/blorenz/COSMOS/COSMOSData/errs.txt'
#Read in the scale of the lines 
err_df = ascii.read(errdatapath,data_start=1,header_start=0,format='csv').to_pandas()


#Read the datafile:
fluxdata = ascii.read(fluxdatapath).to_pandas()

#Fontsizes for plotting
axisfont = 24
ticksize = 18
ticks = 8
titlefont = 24
legendfont = 16
textfont = 16

#Division function
def divz(X,Y):
        return X/np.where(Y,Y,Y+1)*np.not_equal(Y,0)

#Lines that are needed to compute R23
lines = ['3727','4861','5007','4959','6583_fix','6563_fix']
O2 = lines[0]
Hb = lines[1]
O3_5007 = lines[2]
O3_4959 = lines[3]
N2 = lines[4]
Ha = lines[5]

#Filter the data
goodlines = [dataqual[line+'_good'].map(d) for line in lines]
#Needs to be good in all lines to be good
allgood = np.logical_and.reduce(goodlines)
#Needs to be bad in any line to be bad
badlines = [dataqual[line+'_bad'].map(d) for line in lines]
baddata = np.logical_or.reduce(badlines)
lowlines = [dataqual[line+'_low'].map(d) for line in lines]
#Needs to be low in any line to be low, and also not bad in a line
somelow = np.logical_and(np.logical_or.reduce(lowlines),np.logical_not(baddata))

#Function to divide the line by its scale to get the true flux
def scaleline(dataframe,line):
    newflux = divz(dataframe[line+'_flux'],dataframe[line+'_scale'])
    return newflux

def getR23(pd_df,err_df,Ha,Hb,N2,O2,O3_5007,O3_4959):
    scales = pd.DataFrame()
    for line in [Ha,Hb,N2,O2,O3_5007,O3_4959]:
        scales[line] = scaleline(pd_df,line)
    R23num = scales[O3_5007]+scales[O3_4959]+scales[O2]
    R23 = np.log10(divz(R23num,scales[Hb]))
    HaNII = np.log10(divz(scales[N2],scales[Ha]))
    eR23num = err_df[O3_5007]+err_df[O3_4959]+err_df[O2]
    eR23 = (1/np.log(10))*divz(scales[Hb],R23num)*np.sqrt((divz(1,scales[Hb]) * eR23num)**2 + (divz(-R23num,(scales[Hb]**2)) * err_df[Hb])**2)
    eHaNII = (1/np.log(10))*divz(scales[Ha],scales[N2])*np.sqrt((divz(1,scales[Ha]) * err_df[N2])**2 + (divz(-scales[N2],(scales[Ha]**2)) * err_df[Ha])**2)
    return (HaNII,R23,eHaNII,eR23)


plotframe = pd.DataFrame()
plotframe['fluxfile'] = fluxdata['fluxfile']

#Get the x and y axes for the good, low, and bad data
for i in range(0,3):
    if i==0:
        filt = allgood
        colname = 'good'
    elif i==1:
        filt = somelow
        colname = 'low'
    else:
        filt = baddata
        colname = 'bad'
    R23s = getR23(fluxdata[filt],err_df[filt],Ha,Hb,N2,O2,O3_5007,O3_4959)
    plotframe[colname+'x'] = R23s[0]
    plotframe[colname+'y'] = R23s[1]
    plotframe[colname+'ex'] = R23s[2]
    plotframe[colname+'ey'] = R23s[3]


#More plot properties
lw=0.5
mark='o'
ms=3

#Make the figures
fig,axarr = plt.subplots(1,3,figsize=(24,7),sharex=True,sharey=True)
#Plot the data with error bars
c = 0
for ax in axarr:
    if c==0:
        col = 'good'
        color = 'blue'
    elif c==1:
        col = 'low'
        color = 'orange'
    else:
        col = 'bad'
        color = 'red'
    p = plotframe
    ax.errorbar(p[col+'x'],p[col+'y'],xerr=p[col+'ex'],yerr=p[col+'ey'],ls='None',ms=ms,marker=mark,lw=lw,color=color)
    #Titles, axes, legends
    if c==0:
        ax.set_ylabel('log(R_{23})',fontsize = axisfont)
    ax.set_xlabel('log(N[II] / H$\\alpha$)',fontsize = axisfont)
    ax.tick_params(labelsize = ticksize, size=ticks)
    ax.set_xlim(-2.5,0)
    ax.set_ylim(-0.2,1.2)
    #Plot lines that divide the upper and lower R23 branches
    ax.plot((-1.1,-1.1),(-10,10),c='black',ls='--')
    ax.plot((-1.3,-1.3),(-10,10),c='black',ls='--')
    c = c+1
fig.tight_layout()
fig.savefig(figout + 'R23.pdf')
plt.close(fig)

