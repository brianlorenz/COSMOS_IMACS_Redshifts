#Creates a BPT diagram for all objects, and a second figure that shows objects for which single lines are low
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import sys, os, string
import pandas as pd
from astropy.io import fits
import collections
from cycler import cycler
from matplotlib import cm


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


#BPT
#Strings of all of the lines needed for bpt
lines = ['6563_fix','6583_fix']
Ha = lines[0]
N2 = lines[1]
Hb = 0
O3 = 0


#Set low lines to upper limits
#for line in lines:
#    fluxdata.loc[dataqual[line+'_low'],line+'_flux'] = err_df[dataqual[line+'_low']][line+'_err']


#fig2,axarr2 = plt.subplots(2,2,figsize=(15,12))
#ax1,ax2,ax3,ax4 = axarr2[0,0],axarr2[0,1],axarr2[1,0],axarr2[1,1]

#Takes the dataframe and the four lines to combine into the bpt
def getbpt(pd_df,err_df,Ha,Hb,N2,O3):
    errHa = err_df[Ha]
    errN2 = err_df[N2]
    #Divide by the scale to calibrate the flux
    calHa = divz(pd_df[Ha+'_flux'],pd_df[Ha+'_scale'])
    calN2 = divz(pd_df[N2+'_flux'],pd_df[N2+'_scale'])
    #Find the ratios
    Harat = (divz(calN2,calHa))
    #Find the errors
    eHarat = divz(calHa,calN2)*np.sqrt((divz(1,calHa) * errN2)**2 + (divz(-calN2,(calHa**2)) * errHa)**2)
    #Get their logs
    lHarat = np.log10(Harat)
    #Get the log errors
    elHaratu = np.log10(Harat+eHarat)-lHarat
    ldiff = Harat-eHarat
    elHaratd = lHarat-np.log10(ldiff)
    elHaratd.loc[elHaratd<0] = 4
    elHaratd.fillna(4,inplace=True)
    return (lHarat,0,elHaratu,0,elHaratd,0)

#Plotting parameters
ms = 3
lw=0.5
mark='o'

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


dupids = [item for item, count in collections.Counter(fluxdata[allgood]['OBJID']).items() if count > 1]
dupobjsarr = []
for obj in dupids:
    dupobjsarr.append(fluxdata[fluxdata.OBJID == obj].OBJID)
droparr = []
b = np.ones(1216,dtype=bool)
for i in range(0,len(dupobjsarr)):
    #fluxdata = fluxdata.drop(dupobjsarr[i].index.values[1:][0])
    droparr.append(dupobjsarr[i].index.values[1:][0])
b[droparr] = False
allgood = np.logical_and(b,allgood)


plotframe = pd.DataFrame()
plotframe['fluxfile'] = fluxdata['fluxfile']
plotframe['OBJID'] = fluxdata['OBJID']
plotframe['Mask'] = fluxdata['Mask']

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
    bpts = getbpt(fluxdata[filt],err_df[filt],Ha,Hb,N2,O3)
    plotframe[colname+'x'] = bpts[0]
    plotframe[colname+'y'] = bpts[1]
    plotframe[colname+'exu'] = bpts[2]
    plotframe[colname+'eyu'] = bpts[3]
    plotframe[colname+'exd'] = bpts[4]
    plotframe[colname+'eyd'] = bpts[5]


    
#Line that defines the agn cutoff
xline = np.arange(-3.0,0.469,0.001)
yline = 0.61/(xline-0.47)+1.19 #Kewley (2001)
xlineemp = np.arange(-3.0,0.049,0.001)
ylineemp = 0.61/(xlineemp-0.05)+1.3 #Kauffman (2003)





#If duplicates, set to 1, otherwise 0
dup = 0






#Make the figures
fig,axarr = plt.subplots(1,1,figsize=(8,7),sharex=True,sharey=True)
fig2,axarr2 = plt.subplots(1,1,figsize=(9,7),sharex=True,sharey=True)

#Color wheel if looking at duplicates
prop_iter = iter(plt.rcParams['axes.prop_cycle'])

colorsmask=['r','b','g','y','orange','purple','cyan','lime','gold','brown','pink']
masks=['a6','b6','d6','e6','f6','g6','h6','i6','b7','e7','j7']




#Plot the data with error bars
c = 0
#for ax in axarr2:
if 1==1:
    ax=axarr2
    if c==0:
        col = 'good'
        filt = allgood
    elif c==1:
        col = 'low'
        filt = somelow
        sys.exit()
    else:
        col = 'bad'
        filt = baddata
    #filt = np.logical_and(filt,plotframe['goodex']<0.5)
    #filt = np.logical_and(filt,plotframe['goodey']<0.5)
    p = plotframe
    #ax.errorbar(fluxdata[somelow]['LMASS'],p[somelow]['lowx'],yerr=np.array([p[somelow]['lowexd'],p[somelow]['lowexu']]),ls='None',ms=ms,marker=mark,lw=lw,zorder=1,label='Low S/N',color='grey')
    ax.errorbar(fluxdata[filt]['LMASS'],p[filt][col+'x'],yerr=np.array([p[filt][col+'exd'],p[filt][col+'exu']]),ls='None',ms=ms,marker=mark,lw=lw,zorder=1,label='High S/N',color='blue')
    print len(p[filt][col+'x'])
    print len(p[somelow]['lowx'])
    #Titles, axes, legends
    ax.set_xlabel('log(Stellar Mass) (M$_\odot)$',fontsize = axisfont)
    ax.set_ylabel('log(N[II] 6583 / H$\\alpha$)',fontsize = axisfont)
    ax.tick_params(labelsize = ticksize, size=ticks)
    ax.set_ylim(-1.5,0.5)
    ax.set_xlim(8.95,10.05)
    #if c==0: ax.legend(fontsize=legendfont)
    c = c+1
fig2.tight_layout()
fig2.savefig(figout + 'NIIHa_M.pdf')
plt.close(fig2)

