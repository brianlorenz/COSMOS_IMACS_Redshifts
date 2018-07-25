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





'''
#fluxdiff
c = 0
lines = ['4861','4959','5007','6563','4340','4102','6548','6583','3727','6563_fix','6548_fix','6583_fix']
lines = np.sort(lines)
#We want good Ha, HB, and OIII
idxarr = [(fluxdata[line + '_flag']==0) for line in lines]

goodfluxarr = [fluxdata[idx] for idx in idxarr]
badfluxarr = [fluxdata[np.logical_not(idx)] for idx in idxarr]

fig,axarr = plt.subplots(3,4,figsize=(40,20))
axarr = np.reshape(axarr,12)

fig2,axarr2 = plt.subplots(3,4,figsize=(40,20))
axarr2 = np.reshape(axarr2,12)

fig3,axarr3 = plt.subplots(3,4,figsize=(40,20))
axarr3 = np.reshape(axarr3,12)

mark = 'o'
marksize = 4


for ax in axarr:
    xdata = divz(goodfluxarr[c][lines[c] + '_flux'],goodfluxarr[c][lines[c] + '_scale'])
    ydata = divz(goodfluxarr[c][lines[c] + '_sumflux'],goodfluxarr[c][lines[c] + '_scale'])
    ax.plot(xdata,ydata,color='cornflowerblue',marker=mark,ms=marksize,ls='None')
    ax.set_ylabel(lines[c] + ' Unweighted Flux',fontsize = axisfont)
    ax.set_title('Flux Comparison for ' + lines[c],fontsize = titlefont)
    ax.set_xlabel(lines[c] + ' Gaussian Flux',fontsize = axisfont)
    ax.plot((0,1000),(0,1000),color='black',ls='--')
    mad = mad_df[lines[c] + '_mad']
    ax.plot((mad,0),(mad,mad),color='orange',ls='-',lw=2)
    ax.plot((mad,mad),(mad,0),color='orange',ls='-',lw=2)
    xmin,xmax,ymin,ymax = (0.001,np.max(xdata)*1.1,0.001,np.max(ydata)*1.1)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)
    ax.tick_params(labelsize = ticksize)
    c = c+1

c = 0
for ax in axarr2:
    xdata = divz(goodfluxarr[c][lines[c] + '_flux'],goodfluxarr[c][lines[c] + '_scale'])
    ydata = xdata-divz(goodfluxarr[c][lines[c] + '_sumflux'],goodfluxarr[c][lines[c] + '_scale'])
    ax.plot(xdata,ydata,color='cornflowerblue',marker=mark,ms=marksize,ls='None')
    ax.set_ylabel(lines[c] + ' Gaussian Flux - Unweighted Flux',fontsize = axisfont)
    ax.set_title('Flux Difference for ' + lines[c],fontsize = titlefont)
    ax.set_xlabel(lines[c] + ' Gaussian Flux',fontsize = axisfont)
    #ax.plot((0,1000),(0,1000),color='black',ls='--')
    mad = mad_df[lines[c] + '_mad']
    ax.plot((0,100000),(mad,mad),color='orange',ls='-',lw=2)
    #ax.plot((mad,mad),(mad,0),color='orange',ls='-',lw=2)
    xmin,xmax,ymin,ymax = (0.001,np.max(xdata)*1.1,np.min(ydata)*0.9,np.max(ydata)*1.1)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)
    ax.tick_params(labelsize = ticksize)
    c = c+1

c = 0
#Separate by 3* mad?
sepsn = True
#usig or wsig?
usig = False
if usig: sigstr = 'u'
else: sigstr = 'w'

bins = (np.arange(30)-15)*2

for ax in axarr3:
    snidx = divz(goodfluxarr[c][lines[c] + '_flux'],goodfluxarr[c][lines[c] + '_scale']) > 3*mad_df[lines[c] + '_mad'][0]
    highsn = goodfluxarr[c][snidx]
    lowsn = goodfluxarr[c][np.logical_not(snidx)]
    if sepsn:
        highsnflux = divz(highsn[lines[c] + '_flux'],highsn[lines[c] + '_scale'])
        highsnsumflux = divz(highsn[lines[c] + '_sumflux'],highsn[lines[c] + '_scale'])
        lowsnflux = divz(lowsn[lines[c] + '_flux'],lowsn[lines[c] + '_scale'])
        lowsnsumflux = divz(lowsn[lines[c] + '_sumflux'],lowsn[lines[c] + '_scale'])
        xdata = divz(highsnflux-highsnsumflux,highsn[lines[c] + '_' + sigstr + 'sig'])
        xdatalow = divz(lowsnflux-lowsnsumflux,lowsn[lines[c] + '_' + sigstr + 'sig'])
        sepstr = '_sep'
    else:
        xdata = divz(goodfluxarr[c][lines[c] + '_flux']-goodfluxarr[c][lines[c] + '_sumflux'],goodfluxarr[c][lines[c] + '_' + sigstr + 'sig'])
        sepstr = ''
    xdata = xdata[np.abs(xdata) < 30]
    
    if sepsn:
        ax.hist(xdatalow,color='C1',bins=bins,alpha=0.5,label='Low S/N')
        ax.hist(xdata,color='C0',bins=bins,alpha=0.5,label='High S/N')
        ax.legend()
    else: ax.hist(xdata,color='grey',bins=bins)
    ax.set_xlabel(lines[c] + ' (Gaussian Flux - Unweighted Flux)/' + sigstr + 'sig',fontsize = axisfont)
    ax.set_title('Error Histogram for ' + lines[c],fontsize = titlefont)
    ax.set_ylabel('Counts',fontsize = axisfont)
    #ax.plot((0,1000),(0,1000),color='black',ls='--')
    #mad = mad_df[lines[c] + '_mad']
    #ax.plot((0,100000),(mad,mad),color='orange',ls='-',lw=2)
    #ax.plot((mad,mad),(mad,0),color='orange',ls='-',lw=2)
    xmin,xmax = (-30,30)
    #ax.set_xscale('log')
    ax.set_xlim(xmin,xmax)
    #ax.set_ylim(ymin,ymax)
    ax.tick_params(labelsize = ticksize)
    c = c+1


#Titles, axes, legends

fig.tight_layout()
fig2.tight_layout()
fig3.tight_layout()
fig.savefig(figout + 'flux_sumflux.pdf')
fig2.savefig(figout + 'fluxdiff_flux.pdf')
fig3.savefig(figout + 'fluxdiff_hist' + sepstr + '_' + sigstr + 'sig.pdf')
plt.close(fig)
plt.close(fig2)
plt.close(fig3)
'''
