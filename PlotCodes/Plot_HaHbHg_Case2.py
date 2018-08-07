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

#The location to store the scale and its stddev of each line
qualdatapath = '/Users/blorenz/COSMOS/COSMOSData/dataqual.txt'
#Read in the scale of the lines 
dataqual = ascii.read(qualdatapath).to_pandas()

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



#Ha, Hb, Hg
lines = ['6563_fix','4861','4340']

d = {'True': True, 'False': False}

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

highsn=fluxdata[allgood]

errHa = err_df['6563_fix']
errHb = err_df['4861']
errHg = err_df['4340']

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

ax.errorbar(xdata,ydata,xerr=xerr.dropna(),yerr=yerr.dropna(),color='blue',ls='None',ms=ms,marker=mark,lw=0.5,zorder=1)
ax2.errorbar(np.log10(xdata2/CaseBHaHb),np.log10(ydata2/CaseBHgHb),xerr=0.434*divz((xerr2.dropna()/CaseBHaHb),xdata2),yerr=0.434*divz((yerr2.dropna()/CaseBHgHb),ydata2),color='blue',ls='None',ms=ms,marker=mark,lw=0.5,zorder=1)

def Calzetti(wave,Av):
    waveum = wave*.0001
    Rv = 4.05 #+-0.80
    if ((waveum >= 0.63) and (waveum <= 2.2)):
        k = 2.659*(-1.857+divz(1.040,waveum))+Rv
    if ((waveum < 0.63) and (waveum >= 0.12)):
        k = 2.659*(-2.156+divz(1.509,waveum)-divz(0.198,waveum)**2+divz(0.011,waveum)**3)+Rv
    F_rat = 10**(divz(-0.4*k*Av,Rv))
    return F_rat


ax.plot(CaseBHbHg,CaseBHaHb,color='red',marker='x',ms=15,zorder=2)
ax2.plot(0,0,color='red',marker='x',ms=15,zorder=2)

Avs = np.arange(0,3,0.001)

Ha_att = Calzetti(6563,Avs)
Hb_att = Calzetti(4861,Avs)
Hg_att = Calzetti(4340,Avs)

HaHb_att = divz(Ha_att,Hb_att)
HbHg_att = divz(Hb_att,Hg_att)
HgHb_att = divz(Hg_att,Hb_att)

ax.plot(CaseBHbHg*HbHg_att,CaseBHaHb*HaHb_att,color='red',ls='-',label='Case B with extinction',lw=2,zorder=3)
ax2.plot(np.log10(HaHb_att),np.log10(HgHb_att),color='red',ls='-',label='Case B with extinction',lw=2,zorder=3)


#Titles, axes, legends
ax.set_ylabel('H$\\alpha$ / H$\\beta$',fontsize = axisfont)
ax2.set_ylabel('log((H$\\gamma$ / H$\\beta$)/(H$\\gamma$ / H$\\beta$),)',fontsize = axisfont)
ax.set_xlabel('H$\\beta$ / H$\\gamma$',fontsize = axisfont)
ax2.set_xlabel('log((H$\\alpha$ / H$\\beta$)/(H$\\alpha$ / H$\\beta$),)',fontsize = axisfont)
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
