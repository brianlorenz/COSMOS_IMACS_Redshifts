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


#Ha, Hb, Hg
lines = ['6563_fix','4861','4340']
idxarr =  [(fluxdata[line + '_flag']==0) for line in lines]
sigidx0 =  (divz(fluxdata['4861_flux'],fluxdata['4861_scale']) > 3*mad_df['4861_mad'][0])
sigidx1 =  (divz(fluxdata['6563_fix_flux'],fluxdata['6563_fix_scale']) > 3*mad_df['6563_fix_mad'][0])
sigidx2 =  (divz(fluxdata['4340_flux'],fluxdata['4340_scale']) > 3*mad_df['4340_mad'][0])
noflagidx = np.logical_and(np.logical_and(idxarr[0],idxarr[1]),idxarr[2])
highsigidx = np.logical_and(np.logical_and(sigidx0,sigidx1),sigidx2)
allidx = np.logical_and(noflagidx,highsigidx)
highsn = fluxdata[allidx]

errHa = mad_df['6563_fix_mad'][0]
errHb = mad_df['4861_mad'][0]
errHg = mad_df['4340_mad'][0]

fig,ax = plt.subplots(figsize=(13,7))
fig2,ax2 = plt.subplots(figsize=(13,7))

CaseBHaHb = 2.86
CaseBHbHg = divz(1,0.468)
CaseBHgHb = 0.468


def getrats(pd_df,Ha=lines[0],Hb=lines[1],Hg=lines[2]):
    calHa = divz(pd_df[Ha+'_flux'],pd_df[Ha+'_scale'])
    calHb = divz(pd_df[Hb+'_flux'],pd_df[Hb+'_scale'])
    calHg = divz(pd_df[Hg+'_flux'],pd_df[Hg+'_scale'])
    rat1 = divz(calHa,calHb)
    rat2 = divz(calHb,calHg)
    erat1 = np.sqrt((divz(1,calHb) * errHa)**2 + (divz(-calHa,(calHb**2)) * errHb)**2)
    erat2 = np.sqrt((divz(1,calHg) * errHb)**2 + (divz(-calHb,(calHg**2)) * errHg)**2)
    return rat1,rat2,erat1,erat2


def getrats2(pd_df,Ha=lines[0],Hb=lines[1],Hg=lines[2]):
    calHa = divz(pd_df[Ha+'_flux'],pd_df[Ha+'_scale'])
    calHb = divz(pd_df[Hb+'_flux'],pd_df[Hb+'_scale'])
    calHg = divz(pd_df[Hg+'_flux'],pd_df[Hg+'_scale'])
    rat1 = divz(calHa,calHb)
    rat2 = divz(calHg,calHb)
    erat1 = np.sqrt((divz(1,calHb) * errHa)**2 + (divz(-calHa,(calHb**2)) * errHb)**2)
    erat2 = np.sqrt((divz(1,calHb) * errHg)**2 + (divz(-calHg,(calHb**2)) * errHb)**2)
    return rat1,rat2,erat1,erat2

ms = 4
mark='o'

ydata,xdata,yerr,xerr = getrats(highsn)
xdata2,ydata2,xerr2,yerr2 = getrats2(highsn)

ax.errorbar(xdata,ydata,xerr,yerr,color='cornflowerblue',ls='None',ms=ms,marker=mark,lw=0.5)
ax2.errorbar(np.log10(xdata2/CaseBHaHb),np.log10(ydata2/CaseBHgHb),0.434*divz((xerr2/CaseBHaHb),xdata2),0.434*divz((yerr2/CaseBHgHb),ydata2),color='cornflowerblue',ls='None',ms=ms,marker=mark,lw=0.5)

def Calzetti(wave,Av):
    waveum = wave*.0001
    Rv = 4.05 #+-0.80
    if ((waveum >= 0.63) and (waveum <= 2.2)):
        k = 2.659*(-1.857+divz(1.040,waveum))+Rv
    if ((waveum < 0.63) and (waveum >= 0.12)):
        k = 2.659*(-2.156+divz(1.509,waveum)-divz(0.198,waveum)**2+divz(0.011,waveum)**3)+Rv
    F_rat = 10**(divz(-0.4*k*Av,Rv))
    return F_rat


ax.plot(CaseBHbHg,CaseBHaHb,color='red',marker='x',ms=10)
ax2.plot(0,0,color='red',marker='x',ms=10)

Avs = np.arange(0,3,0.001)

Ha_att = Calzetti(6563,Avs)
Hb_att = Calzetti(4861,Avs)
Hg_att = Calzetti(4340,Avs)

HaHb_att = divz(Ha_att,Hb_att)
HbHg_att = divz(Hb_att,Hg_att)
HgHb_att = divz(Hg_att,Hb_att)

ax.plot(CaseBHbHg*HbHg_att,CaseBHaHb*HaHb_att,color='red',ls='-',label='Case B with extinction')
ax2.plot(np.log10(HaHb_att),np.log10(HgHb_att),color='red',ls='-',label='Case B with extinction')


#Titles, axes, legends
ax.set_ylabel('H$\\alpha$ / H$\\beta$',fontsize = axisfont)
ax2.set_ylabel('H$\\gamma$ / H$\\beta$',fontsize = axisfont)
ax.set_title('Balmer Decrements',fontsize = titlefont)
ax2.set_title('Balmer Decrements',fontsize = titlefont)
ax.set_xlabel('H$\\beta$ / H$\\gamma$',fontsize = axisfont)
ax2.set_xlabel('H$\\alpha$ / H$\\beta$',fontsize = axisfont)
ax.tick_params(labelsize = ticksize)
#ax.set_xlim(0.001,10)
#ax.set_ylim(0.01,25)
#ax.set_xlim(-2,1)
ax.set_ylim(0,10)
ax2.set_ylim(-0.82,0.75)
ax2.set_xlim(-0.4,0.6)
#ax.set_xscale("log")
#ax.set_yscale("log")
ax.legend(fontsize = axisfont-2)
ax2.legend(fontsize = axisfont-2)
fig.tight_layout()
fig2.tight_layout()
#plt.show()
fig.savefig(figout + 'HaHb_HbHg.pdf')
fig2.savefig(figout + 'HgHb_HaHb.pdf')
plt.close(fig)
plt.close(fig2)
