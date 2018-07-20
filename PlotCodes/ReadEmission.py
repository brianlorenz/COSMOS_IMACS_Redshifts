#Plots a lot of the various quantities that cna be extracted from emission lines
'''

'''

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
#Plot of Haflux vs Hbflux

#Filter to get the good and back liens
flagidx = np.logical_and((fluxdata['6563_flag']==0), (fluxdata['4861_flag']==0))
goodflux = fluxdata[flagidx]
badflux = fluxdata[np.logical_not(flagidx)]

lw=0.25
mark='.'
erHa = mad_df.iloc[0]['6563_mad']
erHb = mad_df.iloc[0]['4861_mad']
fig,ax = plt.subplots(figsize=(13,7))
ax.errorbar(goodflux['6563_flux'],goodflux['4861_flux'],erHa,erHb,color='black',marker=mark,label='Good Fits',ls='none',lw=lw)
ax.errorbar(badflux['6563_flux'],badflux['4861_flux'],erHa,erHb,color='cornflowerblue',marker=mark,label='Flagged Data',ls='none',alpha=0.2,lw=lw)


#Titles, axes, legends
ax.set_title('H$\\beta$ Flux vs H$\\alpha$ Flux',fontsize = titlefont)
ax.legend(fontsize = legendfont)
ax.set_xlabel('H$\\alpha$ Flux ($10^{-17}$ erg/s/$cm^2$)',fontsize = axisfont)
ax.set_ylabel('H$\\beta$ Flux ($10^-{17}$ erg/s/$cm^2$)',fontsize = axisfont)
ax.set_xscale("log")
ax.set_yscale("log")
ax.tick_params(labelsize = ticksize)
plt.show()
fig.savefig(figout + 'HaHB_flux.pdf')
plt.close(fig)
'''



'''
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
plt.shoow()
fig.savefig(figout + 'linewidth_hist.pdf')
fig2.savefig(figout + 'linewidth_hist_tot.pdf')
'''

'''
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
'''

'''
#Line offsets against each other
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



colors = goodflux['5007_mean']-5007*(1+goodflux['zcc'])
p1 = goodflux['6563_mean']-6562.80*(1+goodflux['zcc'])
p2 = goodflux['4861_mean']-4861.3*(1+goodflux['zcc'])
figure = ax.scatter(p1,p2,c=colors)

cb = fig.colorbar(figure)
ax.text(1.145,0.7,'O[III] 4959 Offset $(\AA)$',fontsize = axisfont, transform=ax.transAxes,rotation=90)      

#Titles, axes, legends
ax.set_ylabel('H$\\beta$ Offset $(\AA)$',fontsize = axisfont)
ax.set_title('Offset Comparison',fontsize = titlefont)
ax.set_xlabel('H$\\alpha$ Offset $(\AA)$',fontsize = axisfont)
ax.tick_params(labelsize = ticksize)
cb.ax.tick_params(labelsize = ticksize)
ax.set_xlim(-5,5)
ax.set_ylim(-5,5)
fig.tight_layout()
plt.show()
fig.savefig(figout + 'offset_by_gal.pdf')
plt.close(fig)
##Fit slope to this line and see if it is ratio of the wavelengths
'''

'''
#Plot of Hbflux vs Hbflux of duplicates

minlim= 0.005
maxlim= 100

#Filter to get the good and back liens
flagidx =(fluxdata['4861_flag']==0)
goodflux = fluxdata[flagidx]
badflux = fluxdata[np.logical_not(flagidx)]
dupids = [item for item, count in collections.Counter(goodflux['OBJID']).items() if count > 1]
dupobjsarr = []
hbflux1 = []
hbflux2 = []
for obj in dupids:
    dupobjsarr.append(fluxdata[fluxdata.OBJID == obj])
for i in range(0,len(dupobjsarr)):
    hbflux1.append(divz(dupobjsarr[i].iloc[0]['4861_flux'],dupobjsarr[i].iloc[0]['4861_scale']))
    hbflux2.append(divz(dupobjsarr[i].iloc[1]['4861_flux'],dupobjsarr[i].iloc[1]['4861_scale']))

dupidsb = [item for item, count in collections.Counter(badflux['OBJID']).items() if count > 1]
dupobjsarrb = []
hbflux1b = []
hbflux2b = []
for obj in dupidsb:
    dupobjsarrb.append(fluxdata[fluxdata.OBJID == obj])
for i in range(0,len(dupobjsarrb)):
    hbflux1b.append(divz(dupobjsarrb[i].iloc[0]['4861_flux'],dupobjsarrb[i].iloc[0]['4861_scale']))
    hbflux2b.append(divz(dupobjsarrb[i].iloc[1]['4861_flux'],dupobjsarrb[i].iloc[1]['4861_scale']))

lw=0.25
mark='.'

fig,ax = plt.subplots(figsize=(8,7))
ax.scatter(hbflux1,hbflux2,color='blue',marker=mark,label='Galaxy with good fit')
ax.scatter(hbflux1b,hbflux2b,color='red',marker=mark,label='Galaxy with bad fit')
ax.plot((0,1000),(0,1000),color='black',ls='--')

#Titles, axes, legends
ax.set_title('H$\\beta$ Flux in Duplicate Measurements',fontsize = titlefont)
ax.legend(fontsize = legendfont)
ax.set_xlabel('H$\\beta$ Flux ($10^{-17}$ erg/s/$cm^2$)',fontsize = axisfont)
ax.set_ylabel('H$\\beta$ Flux ($10^-{17}$ erg/s/$cm^2$)',fontsize = axisfont)
ax.text(0.02,0.79,'1.49*MAD:     ' + str(round(mad_df['4861_mad'],2)),fontsize = textfont, transform=ax.transAxes)      
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(minlim,maxlim)
ax.set_ylim(minlim,maxlim)
ax.tick_params(labelsize = ticksize)
plt.show()
fig.savefig(figout + 'HB_dup_flux.pdf')
plt.close(fig)
'''





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
'''



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

#Plot flux vs sumflux for each line
#Plot (flux - sumflux)/formal error possibly by signal bins
#Ha/Hb v Hb/Hg
#Look at AGN spectra, plot them next to each other

