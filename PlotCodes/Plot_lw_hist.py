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


#Histogram of Line widths

#Filter Ha for good fits
fluxarr = []
for line in lines:
    flagidx = np.logical_and((fluxdata[line + '_flag']==0), (fluxdata[line + '_stddev']>1))
    goodflux = fluxdata[flagidx]
    badflux = fluxdata[np.logical_not(flagidx)]
    fluxarr.append(goodflux)

bins = 30
h1 = fluxarr[0][lines[0]+'_stddev']
h2 = fluxarr[1][lines[1]+'_stddev']
h3 = fluxarr[2][lines[2]+'_stddev']
h4 = fluxarr[3][lines[3]+'_stddev']
tot = pd.concat([h1,h2,h3,h4])
                                        
#ax.hist([h1,h2,h3,h4],color=['darkblue','mediumslateblue','dodgerblue','lightblue'],bins=bins,label=[lines[0],lines[1],lines[2],lines[3]])
fig, axarr = plt.subplots(2,2,figsize=(13,7))#,sharex=True)
fig2, ax4 = plt.subplots(figsize=(13,7))
ax0,ax1,ax2,ax3 = axarr[0,0],axarr[1,0],axarr[0,1],axarr[1,1]

ax0.hist(h1,color='cornflowerblue',bins=bins,label=lines[0])
ax1.hist(h2,color='cornflowerblue',bins=bins,label=lines[1])
ax2.hist(h3,color='cornflowerblue',bins=bins,label=lines[2])
ax3.hist(h4,color='cornflowerblue',bins=bins,label=lines[3])
ax4.hist(tot,color='gray',bins=bins,label='Total')

ax4.set_ylabel('Count',fontsize = axisfont)
ax4.set_title('Line widths for all Lines',fontsize = titlefont)
ax4.set_xlabel('Line width $(\AA)$',fontsize = axisfont)
ax4.tick_params(labelsize = ticksize)
ax4.set_xlim(1.5,7.5)
    

    
#Titles, axes, legends
linecount = 0
for ax in [ax0,ax1,ax2,ax3]:
    ax.set_ylabel('Count',fontsize = axisfont)
    ax.set_title('Line widths for ' + lines[linecount],fontsize = titlefont)
    ax.set_xlabel('Line width $(\AA)$',fontsize = axisfont)
    ax.tick_params(labelsize = ticksize)
    ax.set_xlim(1.5,7.5)
    linecount = linecount + 1

fig.tight_layout()
plt.show()
fig.savefig(figout + 'linewidth_hist.pdf')
fig2.savefig(figout + 'linewidth_hist_tot.pdf')


#Line widths against each other
#We want good Ha, HB, and OIII
flagidx = np.logical_and((fluxdata['6563_flag']==0), (fluxdata['6563_stddev']>1))
flagidx = np.logical_and(flagidx,(fluxdata['4861_flag']==0))
flagidx = np.logical_and(flagidx,(fluxdata['5007_flag']==0))
#flagidx = np.logical_and(flagidx,(fluxdata['4959_flag']==0))
#flagidx = np.logical_and(flagidx,(fluxdata['4959_stddev']>1))
flagidx = np.logical_and(flagidx,(fluxdata['4861_stddev']>1))
flagidx = np.logical_and(flagidx,(fluxdata['5007_stddev']>1))
goodflux = fluxdata[flagidx]
badflux = fluxdata[np.logical_not(flagidx)]

                                        
#ax.hist([h1,h2,h3,h4],color=['darkblue','mediumslateblue','dodgerblue','lightblue'],bins=bins,label=[lines[0],lines[1],lines[2],lines[3]])
fig,ax = plt.subplots(figsize=(13,7))

colors = goodflux['5007_stddev']
sizes = goodflux['4959_stddev']
figure = ax.scatter(goodflux['6563_stddev'],goodflux['4861_stddev'],c=colors)

#cax = plt.axes([0.85, 0.1, 0.075, 0.8])
cb = fig.colorbar(figure)
ax.text(1.145,0.7,'O[III] Line Width $(\AA)$',fontsize = axisfont, transform=ax.transAxes,rotation=90)      

#Titles, axes, legends
ax.set_ylabel('H$\\beta$ Line Width $(\AA)$',fontsize = axisfont)
ax.set_title('Line Width Comparison',fontsize = titlefont)
ax.set_xlabel('H$\\alpha$ Line width $(\AA)$',fontsize = axisfont)
ax.tick_params(labelsize = ticksize)
cb.ax.tick_params(labelsize = ticksize)
ax.set_xlim(1.5,7.5)
ax.set_ylim(1.5,7.5)
fig.tight_layout()
plt.show()
fig.savefig(figout + 'linewidth_by_gal.pdf')
plt.close(fig)
