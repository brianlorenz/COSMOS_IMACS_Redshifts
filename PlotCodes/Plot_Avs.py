#Plots the Av magnitude due to the balmer decerment
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import sys, os, string
import pandas as pd
from astropy.io import fits
import collections
from astropy.stats import biweight_midvariance

#Folder to save the figures
figout = '/Users/blorenz/COSMOS/Reports/2018/Images/'

#The location with the file for all of our data
fluxdatapath = '/Users/blorenz/COSMOS/COSMOSData/lineflux_red.txt'

#The location of the muzzin et al data:
mdatapath = '/Users/blorenz/COSMOS/muzzin_data/UVISTA_final_colors_sfrs_v4.1.dat'

#The location to store the scale and its stddev of each line
scaledata = '/Users/blorenz/COSMOS/COSMOSData/scales.txt'
#Read in the scale of the lines 
scale_df = ascii.read(scaledata).to_pandas()

#Location of the equivalent width data
ewdata = '/Users/blorenz/COSMOS/COSMOSData/lineew.txt'
#Read in the ew of the lines 
ew_df = ascii.read(ewdata).to_pandas()

#File with the error array
errdatapath = '/Users/blorenz/COSMOS/COSMOSData/errs.txt'
#Read in the scale of the lines 
err_df = ascii.read(errdatapath,data_start=1,header_start=0,format='csv').to_pandas()

#The location to store the scale and its stddev of each line
qualdatapath = '/Users/blorenz/COSMOS/COSMOSData/dataqual.txt'
#Read in the scale of the lines 
dataqual = ascii.read(qualdatapath).to_pandas()


#Read the datafile:
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


#Create the figure
fig,axarr = plt.subplots(3,3,figsize=(25,22))

#Plotting parameters
ms = 3
lw=0.5
mark='o'


#Set the Rv value - this one is taken for Calzetti
Rv = 4.05 #+-0.80

#Reddening law from Calzetti et al (2000)
def Calzetti_k(wave):
    waveum = wave*.0001
    if ((waveum >= 0.63) and (waveum <= 2.2)):
        k = 2.659*(-1.857+divz(1.040,waveum))+Rv
    if ((waveum < 0.63) and (waveum >= 0.12)):
        k = 2.659*(-2.156+divz(1.509,waveum)-divz(0.198,waveum**2)+divz(0.011,waveum**3))+Rv
    return k

#Finds the ratio and errors in that ratio of any two lines
def getAv(pd_df,err_df,L1,L2,bdec):
    #Calibrate the fluxes by dividing by scale
    calL1 = divz(pd_df[L1+'_flux'],pd_df[L1+'_scale'])
    calL2 = divz(pd_df[L2+'_flux'],pd_df[L2+'_scale'])
    #Find the ratio
    rat = divz(calL1,calL2)
    #Find the error in the ratio
    erat = np.sqrt((divz(1,calL2) * err_df[L1])**2 + (divz(-calL1,(calL2**2)) * err_df[L2])**2)
    #Get the integer of the line
    if len(L1)==8: iL1=int(L1[0:4])
    else: iL1 = int(L1)
    if len(L2)==8: iL2=int(L2[0:4])
    else: iL2 = int(L2)
    #Get the k value for each line
    L1k = Calzetti_k(iL1)
    L2k = Calzetti_k(iL2)
    #Compute the Av
    Av = divz(np.log10(rat/bdec),(0.4*((L2k/Rv)-(L1k/Rv))))
    #And its error
    eAv = divz((1/np.log(10))*divz((erat/bdec),rat),(0.4*((L2k/Rv)-(L1k/Rv))))
    return Av,eAv


d = {'True': True, 'False': False}

lines0 = ['6563_fix','4861']
lines1 = ['4861','4340']
lines2 = ['6563_fix','4340']

c = 0
Av_df = pd.DataFrame()
Av_df['fluxfile'] = fluxdata['fluxfile']
for lines in [lines0,lines1,lines2]:
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
    if c==0:
        bdec=2.86
    elif c==1:
        bdec=2.137
    else:
        bdec=2.86*2.137
    Av,eAv = getAv(fluxdata[allgood],err_df[allgood],lines[0],lines[1],bdec)
    Av_df['Av'+str(c)] = Av
    Av_df['eAv'+str(c)] = eAv
    c=c+1

for i in range(0,3):
    #Make a histogram of the Avs
    bins = (np.arange(30)-15)*0.5
    #Set the axis
    if i==0:
        lab = 'H$\\alpha$/H$\\beta$ Av'
        lab2 = 'H$\\beta$/H$\\gamma$ Av'
        lines = lines0
    elif i==1:
        lab = 'H$\\beta$/H$\\gamma$ Av'
        lab2 = 'H$\\alpha$/H$\\gamma$ Av'
        lines = lines1
    else:
        lab2 = 'H$\\alpha$/H$\\gamma$ Av'
        lab = 'H$\\alpha$/H$\\beta$ Av'
        lines = lines2
    ax = axarr[0,i]
    ax.hist(Av_df['Av'+str(i)].dropna(),bins=bins,color='grey')
    ax.set_xlabel(lab,fontsize=axisfont)
    ax.set_ylabel('Counts',fontsize = axisfont)
    ax.tick_params(labelsize = ticksize)
    ax.set_xlim(-7.5,7.5)

    ax = axarr[1,i]
    if i==2:
        ax.errorbar(Av_df['Av'+str((i+1)%3)],Av_df['Av'+str(i)],yerr=Av_df['eAv'+str(i)],xerr=Av_df['eAv'+str((i+1)%3)],ls='None',lw=lw,marker=mark,ms=ms,color='blue')
    else:
        ax.errorbar(Av_df['Av'+str(i)],Av_df['Av'+str((i+1)%3)],xerr=Av_df['eAv'+str(i)],yerr=Av_df['eAv'+str((i+1)%3)],ls='None',lw=lw,marker=mark,ms=ms,color='blue')
    ax.plot((-1000,1000),(-1000,1000),color='black',ls='--')
    ax.set_xlabel(lab,fontsize=axisfont)
    ax.set_ylabel(lab2,fontsize=axisfont)
    ax.tick_params(labelsize = ticksize)
    ax.set_xlim(-7,7)
    ax.set_ylim(-7,7)

    ax = axarr[2,i]
    ax.errorbar(fluxdata['LMASS'],Av_df['Av'+str(i)],yerr=Av_df['eAv'+str(i)],color='blue',marker='o',ms=ms,lw=lw,ls='None')

    #Filter the data
    goodlines = [dataqual[line+'_good'].map(d) for line in lines]
    #Needs to be good in all lines to be good
    allgood = np.logical_and.reduce(goodlines)
    
    #Find the mass_weighted medians
    def plotbwt(ratio,ax):
        mr1 = (fluxdata[allgood]['LMASS']<9.25)
        mr2 = np.logical_and(fluxdata[allgood]['LMASS']>=9.25,fluxdata[allgood]['LMASS']<9.5)
        mr3 = np.logical_and(fluxdata[allgood]['LMASS']>=9.5,fluxdata[allgood]['LMASS']<9.75)
        mr4 = (fluxdata[allgood]['LMASS']>=9.75)
        med1 = np.median(ratio[allgood][mr1])
        med2 = np.median(ratio[allgood][mr2])
        med3 = np.median(ratio[allgood][mr3])
        med4 = np.median(ratio[allgood][mr4])
        emed1 = np.sqrt(biweight_midvariance(ratio[allgood][mr1]))
        emed2 = np.sqrt(biweight_midvariance(ratio[allgood][mr2]))
        emed3 = np.sqrt(biweight_midvariance(ratio[allgood][mr3]))
        emed4 = np.sqrt(biweight_midvariance(ratio[allgood][mr4]))
        msbw = 12
        lwbw = 4
        ax.errorbar(9.125,med1,yerr=emed1,color='black',marker='o',ms=msbw,lw=lwbw,ls='None')
        ax.errorbar(9.375,med2,yerr=emed2,color='black',marker='o',ms=msbw,lw=lwbw,ls='None')
        ax.errorbar(9.625,med3,yerr=emed3,color='black',marker='o',ms=msbw,lw=lwbw,ls='None')
        ax.errorbar(9.875,med4,yerr=emed4,color='black',marker='o',ms=msbw,lw=lwbw,ls='None')

        coeff = np.polyfit(fluxdata[allgood]['LMASS'],ratio[allgood],1)
        x = np.arange(9,10,0.001)
        fit = coeff[0]*x+coeff[1]
        ax.plot(x,fit,c='red',lw=4)
    plotbwt(Av_df['Av'+str(i)],ax)
    ax.set_ylim(-7,7)
    ax.tick_params(labelsize = ticksize)
    ax.set_xlabel('LMASS',fontsize = axisfont)
    ax.set_ylabel(lab,fontsize = axisfont)
    
    
fig.tight_layout()
#plt.show()
fig.savefig(figout + 'Av_hist.pdf')
plt.close(fig)
