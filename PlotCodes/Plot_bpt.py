#Creates a BPT diagram for all objects, and a second figure that shows objects for which single lines are low
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import sys, os, string
import pandas as pd
from astropy.io import fits
import collections
from cycler import cycler
from matplotlib import cm


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
    return (Harat,Hbrat,eHarat,eHbrat)

#Plotting parameters
ms = 3
lw=0.5
mark='o'

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


#Find the duplicates
dupids = [item for item, count in collections.Counter(fluxdata[allgood]['OBJID']).items() if count > 1]
dupobjsarr = []
for obj in dupids:
    dupobjsarr.append(fluxdata[fluxdata.OBJID == obj].OBJID)

    
#Line that defines the agn cutoff
xline = np.arange(-3.0,0.469,0.001)
yline = 0.61/(xline-0.47)+1.19 #Kewley (2001)
xlineemp = np.arange(-3.0,0.049,0.001)
ylineemp = 0.61/(xlineemp-0.05)+1.3 #Kauffman (2003)





#If duplicates, set to 1, otherwise 0
dup = 0




#Make the figures
if dup:
    fig,dupax = plt.subplots(1,1,figsize=(8.5,7),sharex=True,sharey=True)
    axarr = [1,2]
else: fig,axarr = plt.subplots(1,3,figsize=(24,7),sharex=True,sharey=True)
fig2,axarr2 = plt.subplots(1,3,figsize=(24,7),sharex=True,sharey=True)

#Color wheel if looking at duplicates
prop_iter = iter(plt.rcParams['axes.prop_cycle'])

colorsmask=['r','b','g','y','orange','purple','cyan','lime','gold','brown','pink']
masks=['a6','b6','d6','e6','f6','g6','h6','i6','b7','e7','j7']



#Plot the data with error bars
#Counter
c = 0
for ax in axarr:
    if ax==2:
        break
    if c==0:
        col = 'good'
        color = 'blue'
    elif c==1:
        col = 'low'
        color = 'orange'
    else:
        col = 'bad'
        color = 'red'
    p = plotframe
    if (c==0 and dup):
        ax = dupax
        '''
        for j in range(0,len(colorsmask)):
            ax.plot(-100,-100,color=colorsmask[j],ms=ms,marker=mark,label=masks[j])
        ax.legend()
        '''
    if dup: ax.errorbar(p[col+'x'],p[col+'y'],xerr=p[col+'ex'],yerr=p[col+'ey'],color='grey',ls='None',ms=ms,marker=mark,lw=lw,label=None)
    else: ax.errorbar(p[col+'x'],p[col+'y'],xerr=p[col+'ex'],yerr=p[col+'ey'],color=color,ls='None',ms=ms,marker=mark,lw=lw,label=None)
    #Titles, axes, legends
    if c==0:
        ax.set_ylabel('log(O[III] 5007 / H$\\beta$)',fontsize = axisfont)
    ax.set_xlabel('log(N[II] 6583 / H$\\alpha$)',fontsize = axisfont)
    ax.tick_params(labelsize = ticksize, size=ticks)
    ax.set_xlim(-2,1)
    ax.set_ylim(-1.2,1.5)
    #Plot the bpt line
    ax.plot(xline,yline,color='black',lw=1,ls='--',label=None)
    ax.plot(xlineemp,ylineemp,color='black',lw=1,ls='-',label=None)
    #Draw lines between duplicates:
    if (c==0 and dup):
            for i in range(0,len(dupobjsarr)):
                pdup = p.iloc[dupobjsarr[i].index]
                #colo = [colorsmask[masks.index(pdup.iloc[0]['Mask'])],colorsmask[masks.index(pdup.iloc[1]['Mask'])]]
                #ax.errorbar(pdup[col+'x'],pdup[col+'y'],xerr=pdup[col+'ex'],yerr=pdup[col+'ey'],color=cm.jet(1.*i/len(dupobjsarr)),ls='-',ms=ms,marker=mark,lw=1,zorder=2)
                
                ax.plot(pdup[col+'x'],pdup[col+'y'],color='red',ls='-',ms=ms,marker=mark,lw=1,zorder=2,label=None)
                #ax.errorbar(pdup[col+'x'].iloc[0],pdup[col+'y'].iloc[0],xerr=pdup[col+'ex'].iloc[0],yerr=pdup[col+'ey'].iloc[0],color=colorsmask[masks.index(pdup.iloc[0]['Mask'])],ls='-',ms=ms,marker=mark,lw=1,zorder=2,label=None)
                #ax.errorbar(pdup[col+'x'].iloc[1],pdup[col+'y'].iloc[1],xerr=pdup[col+'ex'].iloc[1],yerr=pdup[col+'ey'].iloc[1],color=colorsmask[masks.index(pdup.iloc[1]['Mask'])],ls='-',ms=ms,marker=mark,lw=1,zorder=2,label=None)
                ax.errorbar(pdup[col+'x'].iloc[0],pdup[col+'y'].iloc[0],xerr=pdup[col+'ex'].iloc[0],yerr=pdup[col+'ey'].iloc[0],color='blue',ls='-',ms=ms,marker=mark,lw=1,zorder=2,label=None)
                ax.errorbar(pdup[col+'x'].iloc[1],pdup[col+'y'].iloc[1],xerr=pdup[col+'ex'].iloc[1],yerr=pdup[col+'ey'].iloc[1],color='blue',ls='-',ms=ms,marker=mark,lw=1,zorder=2,label=None)
    c = c+1
fig.tight_layout()
if dup: fig.savefig(figout + 'bpt_dup.pdf')
else: fig.savefig(figout + 'bpt.pdf')
plt.close(fig)


#Plot the data with error bars
c = 0
for ax in axarr2:
    if c==0:
        col = 'good'
        filt = allgood
    elif c==1:
        col = 'low'
        filt = somelow
    else:
        col = 'bad'
        filt = baddata
    colors = fluxdata[filt]['LMASS']
    p = plotframe
    ax.errorbar(p[filt][col+'x'],p[filt][col+'y'],xerr=p[filt][col+'ex'],yerr=p[filt][col+'ey'],c='black',ls='None',ms=ms,marker=mark,lw=lw,zorder=1)
    colordata = ax.scatter(p[filt][col+'x'],p[filt][col+'y'],c=colors,marker=mark,zorder=3)
    #Titles, axes, legends
    if c==0:
        ax.set_ylabel('log(O[III] 5007 / H$\\beta$)',fontsize = axisfont)
    ax.set_xlabel('log(N[II] 6583 / H$\\alpha$)',fontsize = axisfont)
    ax.tick_params(labelsize = ticksize, size=ticks)
    ax.set_xlim(-2,1)
    ax.set_ylim(-1.2,1.5)
    #Plot the bpt line
    ax.plot(xline,yline,color='black',lw=1,ls='--')
    ax.plot(xlineemp,ylineemp,color='black',lw=1,ls='-')
    c = c+1
cb = fig2.colorbar(colordata)
cb.ax.tick_params(labelsize = ticksize)
fig2.text(0.98,0.76,'log(Stellar Mass) ($M_{sun}$)',fontsize = axisfont,rotation=90)      
fig2.tight_layout()
fig2.savefig(figout + 'bpt_mass.pdf')
plt.close(fig2)
    
'''
alllowerr = geterrarr(fluxdata[alllow],lines,alllowew)
Halowerr = geterrarr(fluxdata[Halow],lines,Halowew)
Hblowerr = geterrarr(fluxdata[Hblow],lines,Hblowew)
N2lowerr = geterrarr(fluxdata[N2low],lines,N2lowew)
O3lowerr = geterrarr(fluxdata[O3low],lines,O3lowew)

alllowx,alllowy,ealllowx,ealllowy = getbpt(fluxdata[alllow],alllowerr)
#Second figure:
counter = 0
for elmt in [ax1,ax2,ax3,ax4]:
    #Plot the bpt line
    elmt.plot(xline,yline,color='black',lw=1)
    elmt.set_ylabel('Log(O[III] 5007 / H$\\beta$)',fontsize = axisfont)
    elmt.set_xlabel('Log(N[II] 6583 / H$\\alpha$)',fontsize = axisfont)
    elmt.tick_params(labelsize = ticksize)
    elmt.set_xlim(-2,1)
    elmt.set_ylim(-1.2,1.5)
    #Depending on which loop it is, set the x and y data to a particular line
    if counter == 0:
        titlestring = 'H$\\alpha$'
        lowlinex,lowliney,elowlinex,elowliney = getbpt(fluxdata[Halow],Halowerr)
    if counter == 1:
        titlestring = 'H$\\beta$'
        lowlinex,lowliney,elowlinex,elowliney = getbpt(fluxdata[Hblow],Hblowerr)
    if counter == 2:
        titlestring = 'N[II]'
        lowlinex,lowliney,elowlinex,elowliney = getbpt(fluxdata[N2low],N2lowerr)
    if counter == 3:
        titlestring = 'O[III}'
        lowlinex,lowliney,elowlinex,elowliney = getbpt(fluxdata[O3low],O3lowerr)
    
    elmt.set_title('BPT Diagram, low ' + titlestring,fontsize = titlefont)
    #Plot the good data
    elmt.errorbar(goodx,goody,xerr=egoodx,yerr=egoody,color='cornflowerblue',ls='None',ms=msg,marker=mark,label='Good in all lines ('+ str(len(goodflux)) + ')',lw=lw)
    #Plot the ones that are low in this line
    elmt.errorbar(lowlinex,lowliney,xerr=elowlinex,yerr=elowliney,color='orange',ls='None',ms=msg,marker=mark,label='Low '  + titlestring + ' ('+ str(len(lowlinex)-len(alllowx)) + ')',lw=lw)
    counter = counter+1
    #highlight the ones that are low in all lines, so they probably should not be trusted
    elmt.errorbar(alllowx,alllowy,xerr=ealllowx,yerr=ealllowy,color='red',ls='None',ms=msg,marker=mark,label='Low in all lines' + ' ('+ str(len(alllowx)) + ')',lw=lw)
    elmt.legend(fontsize = axisfont-6)





fig2.tight_layout()
#plt.show()

fig2.savefig(figout + 'bpt_low_err.pdf')
#plt.show()

plt.close(fig2)
'''
