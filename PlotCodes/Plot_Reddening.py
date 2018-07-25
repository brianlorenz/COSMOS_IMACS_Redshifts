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
legendfont = 14
textfont = 16

#Division function
def divz(X,Y):
        return X/np.where(Y,Y,Y+1)*np.not_equal(Y,0)

def Calzetti(wave,Av,Rv):
    waveum = wave*.0001
    if ((waveum >= 0.63) and (waveum <= 2.2)):
        k = 2.659*(-1.857+divz(1.040,waveum))+Rv
    if ((waveum < 0.63) and (waveum >= 0.12)):
        k = 2.659*(-2.156+divz(1.509,waveum)-divz(0.198,waveum**2)+divz(0.011,waveum**3))+Rv
    F_rat = 10**(divz(-0.4*k*Av,Rv))
    return np.log10(F_rat)

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
    return np.log10(F_rat)


fig,ax = plt.subplots(figsize=(8,7))

Av=1
wavelengths = np.arange(3800,7002,2)
colors = ('red','blue','black')
Rvs = (3.25,4.05,4.85)
lss = ('--','-','-.')
count = 0
for Rv in Rvs:
    cal = [Calzetti(i,Av,Rv) for i in wavelengths]
    car = [Cardelli(i,Av,Rv) for i in wavelengths]
    ax.plot(wavelengths,cal,color=colors[0],ls=lss[count],label='Calzetti, Rv ' + str(Rv))
    ax.plot(wavelengths,car,color=colors[1],ls=lss[count],label='Cardelli, Rv ' + str(Rv))
    count = count+1

   
lw=0.25
mark='.'


legcolor = 'grey'
legend_elements = [Line2D([0], [0],ls=lss[0],color=legcolor, label='Rv = 3.25'),Line2D([0], [0],color=legcolor,ls=lss[1], label='Rv = 4.05'),Line2D([0], [0],ls=lss[2],color=legcolor, label='Rv = 4.85'),Line2D([0], [0],color=colors[0], label='Calzetti'),Line2D([0], [0],color=colors[1], label='Cardelli'),]

ax.legend(handles=legend_elements, fontsize=legendfont)


#Titles, axes, legends
ax.set_title('Reddening Laws, Av=' + str(Av),fontsize = titlefont)
ax.set_xlabel('Wavelength ($\AA$)',fontsize = axisfont)
ax.set_ylabel('log10(Flux Transmitted)',fontsize = axisfont)
ax.tick_params(labelsize = ticksize)
plt.show()
fig.savefig(figout + 'Redlaws.pdf')
plt.close(fig)
