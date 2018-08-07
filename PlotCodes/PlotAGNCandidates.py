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

#Where the calibrated spectra are stored
caldatapath = '/Users/blorenz/COSMOS/COSMOSData/flxFitsFileOut/'

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
axisfont = 18
ticksize = 16
titlefont = 24
legendfont = 16
textfont = 16

#Division function
def divz(X,Y):
        return X/np.where(Y,Y,Y+1)*np.not_equal(Y,0)

'''
#BPT
lines = ['6563','4861','6583','5007']
flagidx = (fluxdata['6563_flag']==0)
flagidx = np.logical_and(flagidx,(fluxdata['4861_flag']==0))
flagidx = np.logical_and(flagidx,(fluxdata['5007_flag']==0))
flagidx = np.logical_and(flagidx,(fluxdata['6583_flag']==0))
badidx = np.logical_not(flagidx)
badflux = fluxdata[badidx]
Hblow = np.logical_not((fluxdata['4861_flux'] > mad_df['4861_mad'][0]*3))
O3low = np.logical_not((fluxdata['5007_flux'] > mad_df['5007_mad'][0]*3))
Halow = np.logical_not((fluxdata['6563_flux'] > mad_df['6563_mad'][0]*3))
N2low = np.logical_not((fluxdata['6583_flux'] > mad_df['6583_mad'][0]*3))
flagidx = np.logical_and(flagidx,np.logical_not(Hblow))
flagidx = np.logical_and(flagidx,np.logical_not(O3low))
flagidx = np.logical_and(flagidx,np.logical_not(Halow))
flagidx = np.logical_and(flagidx,np.logical_not(N2low))
Hblow = np.logical_and(Hblow,np.logical_not(badidx))
O3low = np.logical_and(O3low,np.logical_not(badidx))
Halow = np.logical_and(Halow,np.logical_not(badidx))
N2low = np.logical_and(N2low,np.logical_not(badidx))
alllow = np.logical_and(np.logical_and(Halow,Hblow),np.logical_and(N2low,O3low))
lowoneline = np.logical_and(np.logical_not(flagidx),np.logical_not(badidx))
lowflux = fluxdata[np.logical_and(lowoneline,np.logical_not(alllow))]
goodflux = fluxdata[flagidx]
alllowflux = fluxdata[alllow]
print len(fluxdata)
print len(badflux),len(lowflux),len(goodflux)



#Takes the dataframe and the four lines to combine into the bpt
def getbpt(pd_df,Ha=lines[0],Hb=lines[1],N2=lines[2],O3=lines[3]):
    calHa = divz(pd_df[Ha+'_flux'],pd_df[Ha+'_scale'])
    calHb = divz(pd_df[Hb+'_flux'],pd_df[Hb+'_scale'])
    calN2 = divz(pd_df[N2+'_flux'],pd_df[N2+'_scale'])
    calO3 = divz(pd_df[O3+'_flux'],pd_df[O3+'_scale'])
    Harat = np.log10(divz(calN2,calHa))
    Hbrat = np.log10(divz(calO3,calHb))
    ratidx = (0 > 0.61/(Harat-0.5)+1.3-Hbrat)
    return Harat,Hbrat,ratidx
'''

#BPT
#Strings of all of the lines needed for bpt
lines = ['6563_fix','6583_fix','4861','5007']
Ha = lines[0]
N2 = lines[1]
Hb = lines[2]
O3 = lines[3]


#fig2,axarr2 = plt.subplots(2,2,figsize=(15,12))
#ax1,ax2,ax3,ax4 = axarr2[0,0],axarr2[0,1],axarr2[1,0],axarr2[1,1]

#Takes the dataframe and the four lines to combine into the bpt
def getbpt(pd_df,err_df,Ha,Hb,N2,O3):
    errHa = err_df[Ha]
    errN2 = err_df[N2]
    errHb = err_df[Hb]
    errO3 = err_df[O3]
    #Divide by the scale to calibrate the flux
    calHa = divz(pd_df[Ha+'_flux'],pd_df[Ha+'_scale'])
    calHb = divz(pd_df[Hb+'_flux'],pd_df[Hb+'_scale'])
    calN2 = divz(pd_df[N2+'_flux'],pd_df[N2+'_scale'])
    calO3 = divz(pd_df[O3+'_flux'],pd_df[O3+'_scale'])
    #Find the ratios
    Harat = np.log10(divz(calN2,calHa))
    Hbrat = np.log10(divz(calO3,calHb))
    #Find the errors
    eHarat = (1/np.log(10))*divz(calHa,calN2)*np.sqrt((divz(1,calHa) * errN2)**2 + (divz(-calN2,(calHa**2)) * errHa)**2)
    eHbrat = (1/np.log(10))*divz(calHb,calO3)*np.sqrt((divz(1,calHb) * errO3)**2 + (divz(-calO3,(calHb**2)) * errHb)**2)
    ratidx = (0 > 0.61/(Harat-0.5)+1.3-Hbrat)
    return (Harat,Hbrat,eHarat,eHbrat,ratidx)


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

msb,msl,msg=(5,5,5)
mark='o'

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
    plotframe[colname+'ex'] = bpts[2]
    plotframe[colname+'ey'] = bpts[3]
    plotframe[colname+'idx'] = bpts[4]

alllowx = plotframe['badx']
alllowy = plotframe['bady']
alllowidx = plotframe['badidx']
goodx = plotframe['goodx']
goody = plotframe['goody']
goodidx = plotframe['goodidx']
lowx = plotframe['lowx']
lowy = plotframe['lowy']
lowidx = plotframe['lowidx']

goodobjs = fluxdata[goodidx]
lowobjs = fluxdata[lowidx]
alllowobjs = fluxdata[alllowidx]


for typestr in ['good','low','all_low']:
    #typestr = 'good'
    #typestr = 'low'
    #typestr = 'all_low'

    colors = ['C0','C1','C2','C3','C4','C5','C6','b','C8','C9','navy','lightsalmon','lime','C0','C1','C2','C3','C4','C5','C6','b','C8','C9','navy','lightsalmon','lime']

    gridsize = 4
    fig,axarr = plt.subplots(gridsize,gridsize,figsize=(30,25))
    axarr = np.reshape(axarr,gridsize**2)
    axbpt = axarr[0]

    if typestr == 'good':
        data = goodobjs
        xtype = goodx
        ytype = goody
        idxtype = goodidx
    elif typestr == 'low':
        data = lowobjs
        xtype = lowx
        ytype = lowy
        idxtype = lowidx
    elif typestr == 'all_low':
        data = alllowobjs
        xtype = alllowx
        ytype = alllowy
        idxtype = alllowidx
        
    axbpt.plot(alllowx[alllowidx],alllowy[alllowidx],color='grey',ls='None',ms=msb,marker=mark,label='Low in all lines (' + str(len(fluxdata[baddata])) + ')')
    axbpt.plot(lowx[lowidx],lowy[lowidx],color='grey',ls='None',ms=msl,marker=mark,label='Low in one line ('+ str(len(fluxdata[somelow])) + ')')
    axbpt.plot(goodx[goodidx],goody[goodidx],color='grey',ls='None',ms=msl,marker=mark,label='Low in one line ('+ str(len(fluxdata[allgood])) + ')')
    
    
    counter = 0
    for j in range (0,len(data)):
        obj = data.iloc[j]
        ax = axarr[counter+1]
        #Calculate the wavelength range for the data
        flxdata = fits.open(caldatapath + obj.fluxfile)[0].data
        flxhead = fits.open(caldatapath + obj.fluxfile)[0].header
        #Read in the spectrum and model
        spec = flxdata[0]
        noise = flxdata[1] 
        model = flxdata[3]
        #Calculate the wavelength range for the data
        crval1 = flxhead["crval1"]
        crpix1 = flxhead["crpix1"]
        cdelt1 = flxhead["cdelt1"]
        naxis1 = flxhead["naxis1"]
        dcflag = flxhead["dc-flag"]
        exptime = flxhead['exptime']
        wavelength = (1.0+np.arange(naxis1)-crpix1)*cdelt1 + crval1
        ax.plot(wavelength,spec,color='grey')
        ax.set_xlim(4880,9000)
        ax.set_ylim(-0.2,3*np.median(spec))
        ax.set_title('Spectrum of ' + str(int(obj.OBJID)) + '_' + obj.Mask, fontsize = titlefont, color=colors[counter])
        ax.set_xlabel('Wavelegnth $\AA$')
        ax.set_xlabel('Wavelength ($\AA$)',fontsize = axisfont)
        ax.set_ylabel('Flux ($10^{-17}$ erg/s/${cm}^2/\AA$)',fontsize = axisfont)
        axbpt.plot(xtype[idxtype].iloc[counter],ytype[idxtype].iloc[counter],color=colors[counter],ls='None',ms=msg,marker=mark)
        counter = counter+1
     
     

    xline = np.arange(-3.0,0.49,0.001)
    yline = 0.61/(xline-0.5)+1.3 #Groves and Kewley
    
    
    
    axbpt.plot(xline,yline,color='black',lw=1)
    
    #Titles, axes, legends
    axbpt.set_ylabel('Log(O[III] 5007 / H$\\beta$)',fontsize = axisfont)
    axbpt.set_title('BPT Diagram (' + typestr + ' colored)' ,fontsize = titlefont)
    axbpt.set_xlabel('Log(N[II] 6583 / H$\\alpha$)',fontsize = axisfont)
    axbpt.tick_params(labelsize = ticksize)
    #ax.set_xlim(0.001,10)
    #ax.set_ylim(0.01,25)
    axbpt.set_xlim(-2,1)
    axbpt.set_ylim(-1.2,1.5)
    #ax.set_xscale("log")
    #ax.set_yscale("log")
    #axbpt.legend(fontsize = axisfont-2)
    fig.tight_layout()
    #fig2.tight_layout()
    #plt.show()
    fig.savefig(figout + 'agn_' + typestr + '_cand.pdf')
    plt.close(fig)
