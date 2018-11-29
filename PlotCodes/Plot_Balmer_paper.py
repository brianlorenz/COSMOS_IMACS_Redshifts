import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import sys, os, string
import pandas as pd
from astropy.io import fits
import collections
from matplotlib.lines import Line2D

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
axisfont = 18
ticksize = 16
titlefont = 24
legendfont = 14
textfont = 16

#Division function
def divz(X,Y):
        return X/np.where(Y,Y,Y+1)*np.not_equal(Y,0)

lines = ['6563_fix','4861','4340']

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




step=0



def Calzetti(wave,Av,Rv):
    waveum = wave*.0001
    if ((waveum >= 0.63) and (waveum <= 2.2)):
        k = 2.659*(-1.857+divz(1.040,waveum))+Rv
    if ((waveum < 0.63) and (waveum >= 0.12)):
        k = 2.659*(-2.156+divz(1.509,waveum)-divz(0.198,waveum**2)+divz(0.011,waveum**3))+Rv
    F_rat = 10**(divz(-0.4*k*Av,Rv))
    return F_rat

def Cardelli(wave,Av,Rv):
    waveum = wave*.0001
    waveum = 1/waveum
    if ((waveum >= 0.3) and (waveum <= 1.1)):
        a = 0.574*waveum**1.61
        b = -0.527*waveum**1.61
        k = Rv*a+b
    if ((waveum >= 1.1) and (waveum <= 3.3)):
        y = waveum - 1.82
        a = 1+0.17699*y-0.50447*y**2-0.02427*y**3+0.72085*y**4+0.01979*y**5-0.77530*y**6+0.32999*y**7
        b = 1.41338*y+2.28305*y**2+1.07233*y**3-5.38434*y**4-0.62251*y**5+5.30260*y**6-2.09002*y**7
        k = Rv*a+b
    F_rat = 10**(divz(-0.4*k*Av,Rv))
    return F_rat

allgood = np.logical_and(allgood,fluxdata['4102_flux']>0.5)
fluxspec = fluxdata[allgood].iloc[0]
fitsfile = '/Users/blorenz/COSMOS/COSMOSData/flxFitsFileOut/' + fluxspec['fluxfile']
zcc = fluxspec['zcc']

data = fits.open(fitsfile)[0].data
head = fits.open(fitsfile)[0].header

crval1 = head["crval1"]
crpix1 = head["crpix1"]
cdelt1 = head["cdelt1"]
naxis1 = head["naxis1"]
dcflag = head["dc-flag"]
exptime = head['exptime']
wavelength = (1.0+np.arange(naxis1)-crpix1)*cdelt1 + crval1
spec = data[0]

fig,axarr = plt.subplots(2,1,figsize=(8,8),sharex=False)




   
lw=0.25
mark='.'

legcolor = 'grey'
#legend_elements = [Line2D([0], [0],ls=lss[0],color=legcolor, label='Rv = 3.25'),Line2D([0], [0],color=legcolor,ls=lss[1], label='Rv = 4.05'),Line2D([0], [0],ls=lss[2],color=legcolor, label='Rv = 4.85'),Line2D([0], [0],color=colors[0], label='Calzetti'),Line2D([0], [0],color=colors[1], label='Cardelli'),]

for ax in axarr:
 #ax.legend(handles=legend_elements, fontsize=legendfont, loc=2,framealpha=1)
 if step==0:
     #ax.plot((6562.8,6562.8),(-100,100),color='black',ls='--')
     #ax.plot((4861.3,4861.3),(-100,100),color='black',ls='--')
     #ax.plot((4340.5,4340.5),(-100,100),color='black',ls='--')
     #ax.plot((4101.7,4101.7),(-100,100),color='black',ls='--')
     lams = [6562.8,4861.3,4340.5,4101.7]
     restwave = wavelength/(1+zcc)
     ax.plot(restwave*(1+zcc),spec,color='black')
     Hazone = np.logical_and((restwave<lams[0]+10),(restwave>lams[0]-10))
     Hbzone = np.logical_and((restwave<lams[1]+10),(restwave>lams[1]-10))
     Hgzone = np.logical_and((restwave<lams[2]+10),(restwave>lams[2]-10))
     Hdzone = np.logical_and((restwave<lams[3]+10),(restwave>lams[3]-10))
     lamcolor = 'cornflowerblue'
     ax.plot(restwave[Hazone]*(1+zcc),spec[Hazone],color=lamcolor)
     ax.plot(restwave[Hbzone]*(1+zcc),spec[Hbzone],color=lamcolor)
     ax.plot(restwave[Hgzone]*(1+zcc),spec[Hgzone],color=lamcolor)
     ax.plot(restwave[Hdzone]*(1+zcc),spec[Hdzone],color=lamcolor)
     ax.set_ylim(-0.1,4)
     ax.set_ylabel('Flux ($10^{-17}$ erg/s/${cm}^2/\AA$)',fontsize = axisfont)
 if step>0:
     flxs = [fluxspec['6563_fix_flux'],fluxspec['4861_flux'],fluxspec['4340_flux'],fluxspec['4102_flux']]
     
     err_df = err_df[allgood].iloc[0]
     errs = [err_df['6563_fix'],err_df['4861'],err_df['4340'],err_df['4102']]
     ax.errorbar(lams,flxs,yerr=errs,marker='o',markersize=12,color='cornflowerblue',zorder=10,ls='None')
     ax.set_ylim(-0.25)
     ax.set_ylabel('Line Flux ($10^{-17}$ erg/s/${cm}^2$)',fontsize = axisfont)


 Avs=(0,0.5,1)
 wavelengths = np.arange(3800,7002,2)
 colors = ('red','blue','black')
 Rv = 4.05
 lss = ('-','--','-.')
 count = 0
 cAv=0
 if step==3: cAv=1
 if step==4:
     cAv=fluxspec.av
     dAv1=fluxspec.dav1
     dAv2=fluxspec.dav2
 for Av in [cAv]:
     cal = [Calzetti(i,Av,Rv) for i in wavelengths]
     car = [Cardelli(i,Av,Rv) for i in wavelengths]
     cal = np.array(cal)
     if (step>1) and (Av<1.5):
         calHa = Calzetti(lams[0],Av,Rv)
         calHb = Calzetti(lams[1],Av,Rv)
         calHg = Calzetti(lams[2],Av,Rv)
         calHd = Calzetti(lams[3],Av,Rv)
         Hbpred = (calHb/calHa)*flxs[0]*(1/2.86)
         Hgpred = (calHg/calHa)*flxs[0]*(1/(2.86*2.137))
         Hdpred = (calHd/calHa)*flxs[0]*(1/(2.86/0.259))
         preds = [flxs[0],Hbpred,Hgpred,Hdpred]
         if step==4:
             calHad = Calzetti(lams[0],Av-dAv1,Rv)
             calHbd = Calzetti(lams[1],Av-dAv1,Rv)
             calHgd = Calzetti(lams[2],Av-dAv1,Rv)
             calHdd = Calzetti(lams[3],Av-dAv1,Rv)
             Hbpredd = (calHbd/calHad)*flxs[0]*(1/2.86)
             Hgpredd = (calHgd/calHad)*flxs[0]*(1/(2.86*2.137))
             Hdpredd = (calHdd/calHad)*flxs[0]*(1/(2.86/0.259))
             predsd = [flxs[0],Hbpredd,Hgpredd,Hdpredd]
             errsd = np.abs(np.array(preds)-np.array(predsd))
         wavelengths = np.arange(3800,7002,2)
         cal = [Calzetti(i,Av,4.05) for i in wavelengths]
         shape = (np.array(wavelengths)*np.array(cal))/6562.8

         #ax.plot(wavelengths,shape*flxs[0],color='red',ls='-',marker=None)

             #ax.errorbar(lams,preds,yerr=errsd,color='red',marker='o',ms=6,ls='None',label='Calzetti, Rv ' + str(Rv),zorder=4)
         ax.text(0.02,0.92,'Av = '+str(np.round(cAv,3)) + ' (mag)',fontsize = axisfont, transform=ax.transAxes,color='Red')
         ax.plot(lams,preds,color='red',marker='o',ms=6,ls='None',label='Calzetti, Rv ' + str(Rv),zorder=200)
     #ax.plot(wavelengths,car,color=colors[1],ls=lss[count],label='Cardelli, Rv ' + str(Rv))
     count = count+1

 #ax.text(0.865,0.02,'Av = 1',fontsize = axisfont-2, transform=ax.transAxes)
 ax.set_xlim(3800,7000)
 if step==0: ax.set_xlim(3800*(1+zcc),7000*(1+zcc))
 




 #Titles, axes, legends
 ax.set_xlabel('Observed Wavelength ($\AA$)',fontsize = axisfont)
 #ax.set_title('Reddening Laws, Av=' + str(Av),fontsize = titlefont)
 if step==1: ax.set_xlabel('Rest Wavelength ($\AA$)',fontsize = axisfont)
 step=1
 ax.tick_params(labelsize = ticksize)
fig.tight_layout()
fig.savefig(figout + 'Balmer_examp_'+str(step)+'.pdf')
plt.close(fig)
