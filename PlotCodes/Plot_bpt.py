#Plots a lot of the various quantities that cna be extracted from emission lines
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

#File for the MAD of the difference in flux of duplicates in each line
maddatapath = '/Users/blorenz/COSMOS/COSMOSData/linemad.txt'

#Read in the mad of the lines 
mad_df = ascii.read(maddatapath).to_pandas()


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


#Get the strings of each line
lines = ['4861','4959','5007','6563']#'6548','6583']






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
lowflux = fluxdata[np.logical_and(np.logical_not(flagidx),np.logical_not(badidx))]
goodflux = fluxdata[flagidx]
print len(fluxdata)
print len(badflux),len(lowflux),len(goodflux)

errHa = mad_df['6563_fix_mad'][0]
errN2 = mad_df['6583_fix_mad'][0]
errHb = mad_df['4861_mad'][0]
errO3 = mad_df['5007_mad'][0]

fig,ax = plt.subplots(figsize=(13,7))
fig2,axarr = plt.subplots(2,2,figsize=(15,12))
ax1,ax2,ax3,ax4 = axarr[0,0],axarr[0,1],axarr[1,0],axarr[1,1]

#Takes the dataframe and the four lines to combine into the bpt
def getbpt(pd_df,Ha=lines[0],Hb=lines[1],N2=lines[2],O3=lines[3]):
    calHa = divz(pd_df[Ha+'_flux'],pd_df[Ha+'_scale'])
    calHb = divz(pd_df[Hb+'_flux'],pd_df[Hb+'_scale'])
    calN2 = divz(pd_df[N2+'_flux'],pd_df[N2+'_scale'])
    calO3 = divz(pd_df[O3+'_flux'],pd_df[O3+'_scale'])
    Harat = np.log10(divz(calN2,calHa))
    Hbrat = np.log10(divz(calO3,calHb))
    eHarat = np.sqrt((divz(1,calHa) * errN2)**2 + (divz(-calN2,(calHa**2)) * errHa)**2)
    eHbrat = np.sqrt((divz(1,calHb) * errO3)**2 + (divz(-calO3,(calHb**2)) * errHb)**2)
    return Harat,Hbrat,eHarat,eHbrat

msb,msl,msg=(3,3,3)
lw=0.5
mark='o'

goodx,goody,egoodx,egoody = getbpt(goodflux)
lowx,lowy,elowx,elowy = getbpt(lowflux)
badx,bady,ebadx,ebady = getbpt(badflux)

xline = np.arange(-3.0,0.49,0.001)
yline = 0.61/(xline-0.5)+1.3 #Groves and Kewley

#ax.plot(badx,bady,color='red',ls='None',ms=msb,marker=mark,label='Bad in one line (' + str(len(badflux)) + ')')
ax.errorbar(lowx,lowy,elowx,elowy,color='orange',ls='None',ms=msl,lw=lw,marker=mark,label='Low in one line ('+ str(len(lowflux)) + ')')
ax.errorbar(goodx,goody,egoodx,egoody,color='cornflowerblue',ls='None',ms=msg,marker=mark,label='Good in all lines ('+ str(len(goodflux)) + ')',lw=lw)

alllowx,alllowy,ealllowx,ealllowy = getbpt(fluxdata[alllow])
counter = 0
for elmt in [ax1,ax2,ax3,ax4]:
    elmt.plot(xline,yline,color='black',lw=1)
    elmt.set_ylabel('Log(O[III] 5007 / H$\\beta$)',fontsize = axisfont)
    elmt.set_xlabel('Log(N[II] 6583 / H$\\alpha$)',fontsize = axisfont)
    elmt.tick_params(labelsize = ticksize)
    elmt.set_xlim(-2,1)
    elmt.set_ylim(-1.2,1.5)
    if counter == 0:
        titlestring = 'H$\\alpha$'
        lowlinex,lowliney,elowlinex,elowliney = getbpt(fluxdata[Halow])
    if counter == 1:
        titlestring = 'H$\\beta$'
        lowlinex,lowliney,elowlinex,elowliney = getbpt(fluxdata[Hblow])
    if counter == 2:
        titlestring = 'N[II]'
        lowlinex,lowliney,elowlinex,elowliney = getbpt(fluxdata[N2low])
    if counter == 3:
        titlestring = 'O[III}'
        lowlinex,lowliney,elowlinex,elowliney = getbpt(fluxdata[O3low])
    
    elmt.set_title('BPT Diagram, low ' + titlestring,fontsize = titlefont)
    elmt.errorbar(goodx,goody,egoodx,egoody,color='cornflowerblue',ls='None',ms=msg,marker=mark,label='Good in all lines ('+ str(len(goodflux)) + ')',lw=lw)
    elmt.errorbar(lowlinex,lowliney,elowlinex,elowliney,color='orange',ls='None',ms=msg,marker=mark,label='Low '  + titlestring + ' ('+ str(len(lowlinex)-len(alllowx)) + ')',lw=lw)
    counter = counter+1
    elmt.errorbar(alllowx,alllowy,ealllowx,ealllowy,color='red',ls='None',ms=msg,marker=mark,label='Low in all lines' + ' ('+ str(len(alllowx)) + ')',lw=lw)
    elmt.legend(fontsize = axisfont-6)



ax.plot(xline,yline,color='black',lw=1)

#Titles, axes, legends
ax.set_ylabel('Log(O[III] 5007 / H$\\beta$)',fontsize = axisfont)
ax.set_title('BPT Diagram',fontsize = titlefont)
ax.set_xlabel('Log(N[II] 6583 / H$\\alpha$)',fontsize = axisfont)
ax.tick_params(labelsize = ticksize)
#ax.set_xlim(0.001,10)
#ax.set_ylim(0.01,25)
ax.set_xlim(-2,1)
ax.set_ylim(-1.2,1.5)
#ax.set_xscale("log")
#ax.set_yscale("log")
ax.legend(fontsize = axisfont-2)
fig.tight_layout()
fig2.tight_layout()
plt.show()
fig.savefig(figout + 'bpt_err.pdf')
fig2.savefig(figout + 'bpt_low_err.pdf')
plt.show()
plt.close(fig)
plt.close(fig2)
