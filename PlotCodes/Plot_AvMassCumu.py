#Examines 24um flux for our objects to investigate the location of dust
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import sys, os, string
import pandas as pd
from astropy.io import fits
import collections
from astropy.cosmology import WMAP9 as cosmo
from astropy.stats import biweight_midvariance

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

#File with the structural properties
spropdatapath = '/Users/blorenz/COSMOS/COSMOSData/struct_prop.txt'
#Read in the scale of the lines 
sprop_df = ascii.read(spropdatapath).to_pandas()
sprop_df = sprop_df.rename(columns={'id':'OBJID'})
fluxdata = pd.merge(fluxdata,sprop_df)

#The location with the file for the filter data
filtdatapath = '/Users/blorenz/COSMOS/COSMOSData/all_c_hasinger.txt'
#Read in the data
filtdata = ascii.read(filtdatapath).to_pandas()
mips24data = filtdata[['id','mips24','emips24','Ks','eKs']]
mips24data = mips24data.rename(columns={'id':'OBJID'})

fluxdata = pd.merge(fluxdata,mips24data,on='OBJID',how='inner')
fluxdata = fluxdata.drop_duplicates()
fluxdata = fluxdata.reset_index()


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

lines=['6563_fix','4861']

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




#Plot the data with error bars
#Counter
c = 0
ylabel = ''



ms=12
lwbw=2

notbad = np.logical_not(baddata)

'''
fig,axarr = plt.subplots(1,3,figsize=(24,7),sharex=True,sharey=True)

for ax in axarr:
    if c in [0,3]:
        col = 'good'
        filt = allgood
        color='blue'
    elif c in [1,4]:
        col = 'low'
        filt = somelow
        color='orange'
    else:
        col = 'bad'
        filt = baddata
        color='red'
    filt = np.logical_and(filt,filtmuz)
    #k24 = fluxdata[filt]['Ks'] - fluxdata[filt]['mips24']
    k24 = fluxdata[filt]['Ksmag'] - fluxdata[filt]['mips24mag']
    ax.plot(k24,fluxdata[filt]['UmV'],color=color,marker='o',ms=4,lw=0.5,ls='None')
    #ax.plot((-100,0.69),(1.3,1.3),color='black')
    #ax.plot((1.5,1.5),(2.01,100),color='black')
    #xline = np.arange(0.69,1.5,0.001)
    #yline = xline*0.88+0.69
    #ax.plot(xline,yline,color='black')
    #Titles, axes, legends
    ax.set_ylabel('U-V',fontsize = axisfont)
    ax.set_xlabel('Ks (mag) - [24]',fontsize = axisfont)
    ax.tick_params(labelsize = ticksize, size=ticks)
    ax.set_xlim(-3,3)
    ax.set_ylim(0,2.5)
    c=c+1

fig.tight_layout()
fig.savefig(figout + 'UV_Ks24.pdf')
plt.close(fig)
'''


cutvar = 'LMASS'
cutvarname = 'LMASS'
cutoff = 9.5
cutoff2 = 0
ycutname = 'LMASS'





fig,axarr = plt.subplots(2,1,figsize=(9,15),sharex=False)
axarr 

c=0
for ax in axarr:
    color='dodgerblue'
    color2='darkblue'
    color3= 'blue'
    if c in [0,1]:
        col = 'good'
        filt = notbad
    elif c in [1,4]:
        col = 'low'
        filt = somelow
    else:
        col = 'bad'
        filt = baddata
    filtup = np.logical_and(filt,fluxdata[cutvar]>=cutoff)
    filtdown = np.logical_and(filt,fluxdata[cutvar]<cutoff)
    if cutoff2 != 0:
        filtmid = np.logical_and(fluxdata[cutvar]<cutoff,fluxdata[cutvar]>=cutoff2)
        filtmid = np.logical_and(filt,filtmid)
        filtdown = np.logical_and(filt,fluxdata[cutvar]<cutoff2)
    if c in [0]:
        ax.errorbar(fluxdata[filtup]['LMASS'],fluxdata[filtup]['av'],yerr=fluxdata[filtup]['dav1'],color=color,marker='o',ms=4,lw=0.5,ls='None')
        ax.errorbar(fluxdata[filtdown]['LMASS'],fluxdata[filtdown]['av'],yerr=fluxdata[filtdown]['dav1'],color=color2,marker='o',ms=4,lw=0.5,ls='None')
        if cutoff2 != 0:
            ax.errorbar(fluxdata[filtmid]['av'],fluxdata[filtmid]['k24'],yerr=fluxdata[filtmid]['ek24'],xerr=fluxdata[filtmid]['dav1'],color=color3,marker='o',ms=4,lw=0.5,ls='None')
            ax.plot((-100,100),(cutoff2,cutoff2),color='black',ls='--')
        ax.set_xlim(8.95,10.05)
        ax.set_ylim(0,4.5)
        ax.set_xlabel('log(Mass) (M$_{sun}$)',fontsize = axisfont)
        ax.set_ylabel('Av (mag)',fontsize = axisfont)
        ax.plot((-100,100),(cutoff,cutoff),color='black',ls='--')
    else:
        ydistup = np.arange(len(fluxdata['av'][filtup]))/float(len(fluxdata['av'][filtup]))
        ydistdown = np.arange(len(fluxdata['av'][filtdown]))/float(len(fluxdata['av'][filtdown]))
        if cutoff2 !=0:
            ydistmid = np.arange(len(fluxdata['av'][filtmid]))/float(len(fluxdata['av'][filtmid]))
            xdistmid = np.sort(fluxdata['av'][filtmid])
        xdistup = np.sort(fluxdata['av'][filtup])
        xdistdown = np.sort(fluxdata['av'][filtdown])
        ax.plot(xdistup,ydistup,color=color,lw=2,label=ycutname + ' > ' + str(cutoff))
        if cutoff2 !=0:  ax.plot(xdistmid,ydistmid,color=color3,lw=2,label=str(cutoff2) + ' < ' + ycutname + ' < ' + str(cutoff))
        ax.plot(xdistdown,ydistdown,color=color2,lw=2,label=ycutname + ' < ' + str(cutoff))
        
        ax.set_ylabel('Cumulative Distribution',fontsize = axisfont)
        ax.set_xlabel('Av (mag)',fontsize = axisfont)
        ax.set_ylim(0,1)
        ax.set_xlim(0,4.5)
        ax.legend(loc=4,fontsize=axisfont-2)
    #ax.plot((-100,0.69),(1.3,1.3),color='black')
    #ax.plot((1.5,1.5),(2.01,100),color='black')
    #xline = np.arange(0.69,1.5,0.001)
    #yline = xline*0.88+0.69
    #ax.plot(xline,yline,color='black')
    #Titles, axes, legends
    ax.tick_params(labelsize = ticksize, size=ticks)
    
    c=c+1

fig.tight_layout()
fig.savefig(figout + 'Av_mass_cumu.pdf')
plt.close(fig)

