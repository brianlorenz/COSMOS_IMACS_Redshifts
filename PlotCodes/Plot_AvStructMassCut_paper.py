#Creates a BPT diagram for all objects, and a second figure that shows objects for which single lines are low
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import sys, os, string
import pandas as pd
from astropy.io import fits
import collections
from astropy.cosmology import WMAP9 as cosmo
from astropy.stats import biweight_midvariance

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

#File with the structural properties
spropdatapath = '/Users/blorenz/COSMOS/COSMOSData/struct_prop.txt'
#Read in the scale of the lines 
sprop_df = ascii.read(spropdatapath).to_pandas()
sprop_df = sprop_df.rename(columns={'id':'OBJID'})
fluxdata = pd.merge(fluxdata,sprop_df)

#Read in the sfr file
sfdata = '/Users/blorenz/COSMOS/COSMOSData/sfrs.txt'
sfr_df = ascii.read(sfdata).to_pandas()
fluxdata = pd.merge(fluxdata,sfr_df,on='fluxfile')

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

lines=['6563_fix','4861']

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


combinemass = 1

if not combinemass:
    fig,axarr = plt.subplots(2,3,figsize=(24,15),sharex=False,sharey=False)
    axarr = np.reshape(axarr,6)
else:
    fig,axarr = plt.subplots(2,2,figsize=(16,15),sharex=False,sharey=False)
    axarr = np.reshape(axarr,4)
#Gets rid of objects with bad ellipticities
filtar = fluxdata['ar']>0


#Plot the data with error bars
#Counter
c = 0
plotdata = 'ar'
ylabel = 'b/a'
savename = 'AxisRatio'

#fluxdata['n']=np.log10(fluxdata['n'])
#fluxdata['SMD']=np.log10(divz(fluxdata['LMASS'],(4*np.pi*fluxdata['re_kpc']**2)))



ms=12
lwbw=2

notbad = np.logical_not(baddata)

#colormap = np.log10(fluxdata['sSFR'])
colormap = fluxdata['re_kpc']
#cut = 2.7
cut = 3.83
colorcut = 1
colorcut1 = 'blue'
colormed1 = 'darkblue'
colorcut2 = 'red'
colormed2 = 'maroon'
propname = 're'

for ax in axarr:
    color1='dodgerblue'
    color3='darkblue'
    color2= 'blue'
    color4='black'
    if c in [0,1,2]:
        massfilt = fluxdata['LMASS']<9.5
    else:
        massfilt = fluxdata['LMASS']>=9.5
    if c in [0,2,3,5]:
        col = 'good'
        filt = notbad
        if combinemass: filt = allgood
        color='blue'
    elif c in [1,4]:
        col = 'low'
        filt = notbad
        if combinemass: filt = (fluxdata.OBJID < 0)
        color='orange'
    else:
        col = 'bad'
        filt = baddata
        color='red'
    #ax.errorbar(fluxdata[filt][filtar]['av'],fluxdata[filt][filtar]['ar'],xerr=fluxdata[filt][filtar]['dav1'],color=color,marker='o',ms=4,lw=0.5,ls='None')
    #ax2.errorbar(fluxdata[filt][filtar]['av'],fluxdata[filt][filtar]['re_kpc'],xerr=fluxdata[filt][filtar]['dav1'],color=color,marker='o',ms=4,lw=0.5,ls='None')
    #Titles, axes, legends
    acount = 0
    filttype = (fluxdata[plotdata]>-98.9)
    if c==0:
        ax.set_ylabel(ylabel+', LMASS < 9.5',fontsize = axisfont)
    if c==3:
        ax.set_ylabel(ylabel+', LMASS >= 9.5',fontsize = axisfont)
    ax.set_xlabel('Av (mag)',fontsize = axisfont)
    ax.tick_params(labelsize = ticksize, size=ticks)
    filters = np.logical_and(filt,massfilt)
    filters = np.logical_and(filters,filttype)
    if c in [0,2,3,5]:
        loc1 = np.sin(22.5/180*np.pi)
        loc2 = np.sin(45.0/180*np.pi)
        loc3 = np.sin(67.5/180*np.pi)
        mr1 = (fluxdata[filters]['ar']<loc1)
        mr2 = np.logical_and(fluxdata[filters]['ar']>=loc1,fluxdata[filters]['ar']<loc2)
        mr3 = np.logical_and(fluxdata[filters]['ar']>=loc2,fluxdata[filters]['ar']<loc3)
        mr4 = (fluxdata[filters]['ar']>=loc3)
        med1 = np.median(fluxdata[filters][mr1].av)
        med2 = np.median(fluxdata[filters][mr2].av)
        med3 = np.median(fluxdata[filters][mr3].av)
        med4 = np.median(fluxdata[filters][mr4].av)
        med751 = np.percentile(fluxdata[filters][mr1].av,75)
        med752 = np.percentile(fluxdata[filters][mr2].av,75)
        med753 = np.percentile(fluxdata[filters][mr3].av,75)
        med754 = np.percentile(fluxdata[filters][mr4].av,75)
        emed1 = np.sqrt(biweight_midvariance(fluxdata[filters][mr1].av))/len(fluxdata[filters][mr1])
        emed2 = np.sqrt(biweight_midvariance(fluxdata[filters][mr2].av))/len(fluxdata[filters][mr2])
        emed3 = np.sqrt(biweight_midvariance(fluxdata[filters][mr3].av))/len(fluxdata[filters][mr3])
        emed4 = np.sqrt(biweight_midvariance(fluxdata[filters][mr4].av))/len(fluxdata[filters][mr4])
        s1 = np.median(fluxdata[filters][mr1]['ar'])
        s2 = np.median(fluxdata[filters][mr2]['ar'])
        s3 = np.median(fluxdata[filters][mr3]['ar'])
        s4 = np.median(fluxdata[filters][mr4]['ar'])
        if c in [0,3]:
            ax.errorbar(fluxdata[filters][mr1]['av'],fluxdata[filters][mr1][plotdata],xerr=fluxdata[filters][mr1]['dav1'],color=color1,marker='o',ms=4,lw=0.5,ls='None',label=None)
            ax.errorbar(fluxdata[filters][mr2]['av'],fluxdata[filters][mr2][plotdata],xerr=fluxdata[filters][mr2]['dav1'],color=color2,marker='o',ms=4,lw=0.5,ls='None',label=None)
            ax.errorbar(fluxdata[filters][mr3]['av'],fluxdata[filters][mr3][plotdata],xerr=fluxdata[filters][mr3]['dav1'],color=color3,marker='o',ms=4,lw=0.5,ls='None',label=None)
            ax.errorbar(fluxdata[filters][mr4]['av'],fluxdata[filters][mr4][plotdata],xerr=fluxdata[filters][mr4]['dav1'],color=color4,marker='o',ms=4,lw=0.5,ls='None',label=None)
            if colorcut:
                #Cut so that we only have SF galaxies
                sf = np.log10(fluxdata.sSFR)>-10.5
                above = (colormap[filters]>cut)
                below = (colormap[filters]<=cut)
                above = np.logical_and(above,sf[filters])
                below = np.logical_and(below,sf[filters])
                ax.errorbar(fluxdata[filters][above]['av'],fluxdata[filters][above][plotdata],xerr=fluxdata[filters][above]['dav1'],color=colorcut1,marker='o',ms=4,lw=0.5,ls='None',label=propname+'>'+str(cut))
                ax.errorbar(fluxdata[filters][below]['av'],fluxdata[filters][below][plotdata],xerr=fluxdata[filters][below]['dav1'],color=colorcut2,marker='o',ms=4,lw=0.5,ls='None',label=propname+'<'+str(cut))
                medsabove = [np.median(fluxdata[filters][np.logical_and(above,g)].av) for g in [mr1,mr2,mr3,mr4]]
                medsbelow = [np.median(fluxdata[filters][np.logical_and(below,g)].av) for g in [mr1,mr2,mr3,mr4]]
                #emedsabove = [np.std(fluxdata[filters][np.logical_and(above,g)].av) for g in [mr1,mr2,mr3,mr4]]
                emedsabove = 1.49*np.array([np.median(np.abs(fluxdata[filters][np.logical_and(above,g)].av-np.median(fluxdata[filters][np.logical_and(above,g)].av))) for g in [mr1,mr2,mr3,mr4]])
                emedsabove = 1.49*np.array([np.median(np.abs(fluxdata[filters][np.logical_and(below,g)].av-np.median(fluxdata[filters][np.logical_and(below,g)].av))) for g in [mr1,mr2,mr3,mr4]])
                #emedsbelow = [np.std(fluxdata[filters][np.logical_and(below,g)].av) for g in [mr1,mr2,mr3,mr4]]
                ax.legend(fontsize=axisfont-6,loc=4)
                s = 12
                ax.errorbar(medsabove,[s1,s2,s3,s4],xerr=emedsabove,label='Median ' + propname + ' > ' + str(cut),ms=s,ls='None',marker='x',zorder=10**10, markerfacecolor='None', markeredgecolor=colormed1,mew=4)
                ax.errorbar(medsbelow,[s1,s2,s3,s4],xerr=emedsbelow,label='Median ' + propname + ' < ' + str(cut),ms=s,ls='None',marker='o',zorder=10**10, markerfacecolor='None', markeredgecolor=colormed2,mew=4)
               
            else:
                ax.errorbar(med1,s1,xerr=emed1,color='red',marker='o',ms=ms,lw=lwbw,ls='None',label=None)
                ax.errorbar(med2,s2,xerr=emed2,color='red',marker='o',ms=ms,lw=lwbw,ls='None',label=None)
                ax.errorbar(med3,s3,xerr=emed3,color='red',marker='o',ms=ms,lw=lwbw,ls='None',label=None)
                ax.errorbar(med4,s4,xerr=emed4,color='red',marker='o',ms=ms,lw=lwbw,ls='None',label='Median in bin')
                ax.text(0.685,0.02,'Median in bin',fontsize = axisfont-2, transform=ax.transAxes,color='red')
            #ax.errorbar(med751,loc1/2,xerr=emed1,color='red',marker='o',ms=ms,lw=lwbw,ls='None')
            #ax.errorbar(med752,(loc1+loc2)/2,xerr=emed2,color='red',marker='o',ms=ms,lw=lwbw,ls='None')
            #ax.errorbar(med753,(loc2+loc3)/2,xerr=emed3,color='red',marker='o',ms=ms,lw=lwbw,ls='None')
            #ax.errorbar(med754,(1+loc3)/2,xerr=emed4,color='red',marker='o',ms=ms,lw=lwbw,ls='None')
            ax.plot((-100,100),(loc1,loc1),color='black',ls='--',label=None)
            ax.plot((-100,100),(loc2,loc2),color='black',ls='--',label=None)
            ax.plot((-100,100),(loc3,loc3),color='black',ls='--',label=None)
           
        ydist1 = np.arange(len(fluxdata[filters][mr1]['av']))/float(len(fluxdata[filters][mr1]['av']))
        xdist1 = np.sort(fluxdata[filters][mr1]['av'])
        ydist2 = np.arange(len(fluxdata[filters][mr2]['av']))/float(len(fluxdata[filters][mr2]['av']))
        xdist2 = np.sort(fluxdata[filters][mr2]['av'])
        ydist3 = np.arange(len(fluxdata[filters][mr3]['av']))/float(len(fluxdata[filters][mr3]['av']))
        xdist3 = np.sort(fluxdata[filters][mr3]['av'])
        ydist4 = np.arange(len(fluxdata[filters][mr4]['av']))/float(len(fluxdata[filters][mr4]['av']))
        xdist4 = np.sort(fluxdata[filters][mr4]['av'])
        if c in [2,5]:
            ax.plot(xdist1,ydist1,color=color1)
            ax.plot(xdist2,ydist2,color=color2)
            ax.plot(xdist3,ydist3,color=color3)
            ax.plot(xdist4,ydist4,color=color4)
            ax.set_ylabel('Cumulative Distribution',fontsize=axisfont)
        
    ax.set_xlim(-0.1,3)
    ax.set_ylim(0,1)
    c = c+1
    if (combinemass and (c in [1,4])): c = c+1


fig.tight_layout()
if colorcut: fig.savefig(figout + 'Av_'+savename+'_'+propname+'_cut.pdf')
elif combinemass: fig.savefig(figout + 'Av_'+savename+'_combmass.pdf')
else:fig.savefig(figout + 'Av_'+savename+'_mass.pdf')
plt.close(fig)




#Color for BoT > 0.1 or 0.2
#Red for BoT <0.2 and r<2.7...
#Remove bulge galaxies since we only want to look at disks
#thinka bout whether to use stddev or se in the mean (stddev/sqrt(n))

#<Av> vs i, arccos(b/a)
