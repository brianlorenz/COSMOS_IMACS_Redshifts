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

#Merge our data with the UVISTA catalog
#The location of the muzzin et al data:
mdatapath = '/Users/blorenz/COSMOS/muzzin_data/UVISTA_final_colors_sfrs_v4.1.dat'

#Read in the muzzin data
mdata = ascii.read(mdatapath).to_pandas()
zidx = np.logical_and(mdata['z_peak'] > 0.3, mdata['z_peak'] < 0.4)
midx = np.logical_and(mdata['LMASS'] > 9, mdata['LMASS'] < 10)
idx = np.logical_and(zidx,midx)
total = mdata[idx]

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



#Ha, Hb, Hg
lines = ['6563','4861','4340']
flagHa = (fluxdata['6563_fix_flag']==0)
flagHb = (fluxdata['4861_flag']==0)
flagHg = (fluxdata['4340_flag']==0)
badHaHb = np.logical_not(np.logical_and(flagHa,flagHb))
badHbHg = np.logical_not(np.logical_and(flagHb,flagHg))
Halow = np.logical_not((fluxdata['6563_fix_flux'] > mad_df['6563_fix_mad'][0]*3))
Hblow = np.logical_not((fluxdata['4861_flux'] > mad_df['4861_mad'][0]*3))
Hglow = np.logical_not((fluxdata['4340_flux'] > mad_df['4340_mad'][0]*3))
HaHblow = np.logical_or(Halow,Hblow)
HbHglow = np.logical_or(Hblow,Hglow)
HaHblow = np.logical_and(HaHblow,np.logical_not(badHaHb))
HbHglow = np.logical_and(HbHglow,np.logical_not(badHbHg))
HaHbgood = fluxdata[np.logical_and(np.logical_not(HaHblow),np.logical_not(badHaHb))]
HbHggood = fluxdata[np.logical_and(np.logical_not(HbHglow),np.logical_not(badHbHg))]



fig,axarr = plt.subplots(2,2,figsize=(20,15))
axarr = np.reshape(axarr,4)
ax1,ax2,ax3,ax4 = axarr[0],axarr[1],axarr[2],axarr[3]

def getrat2line(pd_df,L1,L2):
    errL1 = mad_df[L1+'_mad'][0]
    errL2 = mad_df[L2+'_mad'][0]
    calL1 = divz(pd_df[L1+'_flux'],pd_df[L1+'_scale'])
    calL2 = divz(pd_df[L2+'_flux'],pd_df[L2+'_scale'])  
    rat = divz(calL1,calL2)
    erat = np.sqrt((divz(1,calL2) * errL1)**2 + (divz(-calL1,(calL2**2)) * errL2)**2)
    return rat,erat

ms = 4
mark='o'

HaHb,HaHberr = getrat2line(HaHbgood,'6563_fix','4861')
HbHg,HbHgerr = getrat2line(HaHbgood,'4861','4340')

#ax.errorbar(xdata,ydata,xerr,yerr,color='cornflowerblue',ls='None',ms=ms,marker=mark,lw=0.5)
#ax2.errorbar(np.log10(xdata2/CaseBHaHb),np.log10(ydata2/CaseBHgHb),0.434*divz((xerr2/CaseBHaHb),xdata2),0.434*divz((yerr2/CaseBHgHb),ydata2),color='cornflowerblue',ls='None',ms=ms,marker=mark,lw=0.5)

Rv = 4.05 #+-0.80

def Calzetti_k(wave):
    waveum = wave*.0001
    if ((waveum >= 0.63) and (waveum <= 2.2)):
        k = 2.659*(-1.857+divz(1.040,waveum))+Rv
    if ((waveum < 0.63) and (waveum >= 0.12)):
        k = 2.659*(-2.156+divz(1.509,waveum)-divz(0.198,waveum)**2+divz(0.011,waveum)**3)+Rv
    return k


Ha_k = Calzetti_k(6563)
Hb_k = Calzetti_k(4861)
Hg_k = Calzetti_k(4340)

AvHaHb = divz(np.log10(HaHb/2.86),(0.4*((Hb_k/Rv)-(Ha_k/Rv))))
AvHbHg = divz(np.log10(HbHg/2.86),(0.4*((Hg_k/Rv)-(Hb_k/Rv))))
AvHaHberr = divz(0.434*divz((HaHberr/2.86),HaHb),(0.4*((Hb_k/Rv)-(Ha_k/Rv))))
AvHbHgerr = divz(0.434*divz((HbHgerr/2.86),HbHg),(0.4*((Hg_k/Rv)-(Hb_k/Rv))))

mHaHb = np.median(AvHaHberr)
mHbHg = np.median(AvHbHgerr)

bins = 30
ax1.hist(AvHaHb,bins=bins,color='grey')
ax2.hist(AvHbHg[AvHbHg>-100],bins=bins,color='grey')
ax3.errorbar(AvHaHb,AvHbHg,AvHaHberr,AvHbHgerr,color='cornflowerblue',marker='o',ms=4,lw=0.5,ls='None')

#Titles, axes, legends
ax2.set_xlabel('(H$\\beta$ / H$\\gamma$) Av',fontsize = axisfont)
ax1.set_xlabel('(H$\\alpha$ / H$\\beta$) Av',fontsize = axisfont)
ax3.set_xlabel('(H$\\alpha$ / H$\\beta$) Av',fontsize = axisfont)
ax3.set_ylabel('(H$\\beta$ / H$\\gamma$) Av',fontsize = axisfont)
ax1.text(0.02,0.9,'Median(error): ' + str(np.round(mHaHb,4)),fontsize = axisfont, transform=ax1.transAxes)
ax2.text(0.02,0.9,'Median(error): ' + str(np.round(mHbHg,4)),fontsize = axisfont, transform=ax2.transAxes)
ax3.set_title('Comparison of Av',fontsize = titlefont)
for ax in [ax1,ax2]:
    ax.set_title('Extinction Histogram',fontsize = titlefont)
    ax.tick_params(labelsize = ticksize)
    ax.set_ylabel('Counts',fontsize = axisfont)
#ax.set_xlim(0.001,10)
ax3.set_ylim(-5,5)
ax3.set_xlim(-5,5)
#ax.set_xscale("log")
#ax.set_yscale("log")
#ax.legend(fontsize = axisfont-2)
#ax2.legend(fontsize = axisfont-2)
fig.tight_layout()
#plt.show()
fig.savefig(figout + 'Av_hist.pdf')
plt.close(fig)

#Plot flux vs sumflux for each line
#Plot (flux - sumflux)/formal error possibly by signal bins
#Ha/Hb v Hb/Hg
#Look at AGN spectra, plot them next to each other
