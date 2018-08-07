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
fluxdatapath = '/Users/blorenz/COSMOS/COSMOSData/lineflux_red.txt'

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
#Read in the errors of the lines 
err_df = ascii.read(errdatapath,data_start=1,header_start=0,format='csv').to_pandas()

errreddatapath = '/Users/blorenz/COSMOS/COSMOSData/errs_red.txt'
#Read in the errors of the reddened lines 
err_dfred = ascii.read(errreddatapath,data_start=1,header_start=0,format='csv').to_pandas()



#Read the datafile:
fluxdata = ascii.read(fluxdatapath).to_pandas()

#Check to reddening correct:
red = 0
if len(sys.argv)>1:
    red = sys.argv[1]
    alllines = ['4861','5007','4959','6563_fix','6583_fix','6548_fix','3727','4340','4102','6563','6583','6548']
    
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
lines = ['3727','4861','5007','4959']
O2 = lines[0]
Hb = lines[1]
O3_5007 = lines[2]
O3_4959 = lines[3]

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
def scaleline(dataframe,line,usered):
    if usered:
        newflux = divz(dataframe[line+'_flux_red'],dataframe[line+'_scale'])
    else:
        newflux = divz(dataframe[line+'_flux'],dataframe[line+'_scale'])
    return newflux

def getO32(pd_df,err_df,Hb,O2,O3_5007,O3_4959,usered):
    scales = pd.DataFrame()
    for line in [Hb,O2,O3_5007,O3_4959]:
        scales[line] = scaleline(pd_df,line,usered)
    R23num = scales[O3_5007]+scales[O3_4959]+scales[O2]
    R23 = np.log10(divz(R23num,scales[Hb]))
    O32num = scales[O3_5007]+scales[O3_4959]
    O32 = np.log10(divz(O32num,scales[O2]))
    eR23num = err_df[O3_5007]+err_df[O3_4959]+err_df[O2]
    eR23 = (1/np.log(10))*divz(scales[Hb],R23num)*np.sqrt((divz(1,scales[Hb]) * eR23num)**2 + (divz(-R23num,(scales[Hb]**2)) * err_df[Hb])**2)
    eO32num = err_df[O3_5007]+err_df[O3_4959]
    eO32 = (1/np.log(10))*divz(scales[O2],O32num)*np.sqrt((divz(1,scales[O2]) * eO32num)**2 + (divz(-O32num,(scales[O2]**2)) * err_df[O2])**2)
    return (R23,O32,eR23,eO32)


plotframe = pd.DataFrame()
plotframe['fluxfile'] = fluxdata['fluxfile']
if red:
    pred = pd.DataFrame()
    pred['fluxfile'] = fluxdata['fluxfile']

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
    data = getO32(fluxdata[filt],err_df[filt],Hb,O2,O3_5007,O3_4959,0)
    plotframe[colname+'x'] = data[0]
    plotframe[colname+'y'] = data[1]
    plotframe[colname+'ex'] = data[2]
    plotframe[colname+'ey'] = data[3]
    if red:
        datared = getO32(fluxdata[filt],err_dfred[filt],Hb,O2,O3_5007,O3_4959,1)
        pred[colname+'x'] = datared[0]
        pred[colname+'y'] = datared[1]
        pred[colname+'ex'] = datared[2]
        pred[colname+'ey'] = datared[3]


#More plot properties
lw=0.5
mark='o'
ms=3

#Make the figures
if red:
    fig,axarrs = plt.subplots(2,3,figsize=(24,16),sharex=True,sharey=True)
    axarr = axarrs.reshape(6)
else:
    fig,axarr = plt.subplots(1,3,figsize=(24,7),sharex=True,sharey=True)
#Plot the data with error bars
c = 0
for ax in axarr:
    if (c==0 or c==3):
        col = 'good'
        color = 'blue'
    elif (c==1 or c==4):
        col = 'low'
        color = 'orange'
    else:
        col = 'bad'
        color = 'red'
    p = plotframe
    if c>2:
        ax.errorbar(pred[col+'x'],pred[col+'y'],xerr=pred[col+'ex'],yerr=p[col+'ey'],ls='None',ms=ms,marker=mark,lw=lw,color=color)
        ax.text(0.02,0.96,'Reddening Corrected',fontsize = textfont, transform=ax.transAxes)
    else:
        ax.errorbar(p[col+'x'],p[col+'y'],xerr=p[col+'ex'],yerr=p[col+'ey'],ls='None',ms=ms,marker=mark,lw=lw,color=color)
    #Titles, axes, legends
    if c==0:
        ax.set_ylabel('log($O_{32}$)',fontsize = axisfont)
    if c==3:
        ax.set_ylabel('log($O_{32}$)',fontsize = axisfont)
    ax.set_xlabel('log($R_{23}$)',fontsize = axisfont)
    ax.tick_params(labelsize = ticksize, size=ticks)
    ax.set_xlim(-0.5,1.5)
    ax.set_ylim(-1.25,1.5)
    #ax.plot((-1.1,-1.1),(-10,10),c='black',ls='--')
    c = c+1
fig.tight_layout()
fig.savefig(figout + 'O32_R23.pdf')
plt.close(fig)

