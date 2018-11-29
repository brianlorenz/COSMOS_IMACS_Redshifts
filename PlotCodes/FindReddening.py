#Find which objects are bad and low based on various cuts through the data. Output this as a dataframe containing True False for every object in every line
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import sys, os, string
import pandas as pd
from astropy.io import fits
import collections
from astropy.stats import biweight_midvariance
from matplotlib.font_manager import FontProperties

#Folder to save the figures
figout = '/Users/blorenz/COSMOS/Reports/2018/Images/'

#The location with the file for all of our data
qualout = '/Users/blorenz/COSMOS/COSMOSData/dataqual.txt'

#The location with the file for all of our data
fluxdatapath = '/Users/blorenz/COSMOS/COSMOSData/lineflux_zero.txt'

#The location with the file for all of our data
outpath = '/Users/blorenz/COSMOS/COSMOSData/lineflux_red.txt'

#Path for outdated avs, computed using balmer ratio
avpath = '/Users/blorenz/COSMOS/COSMOSData/balmer_avs.txt'
av_df = ascii.read(avpath).to_pandas()


#The location to store the scale and its stddev of each line
scaledata = '/Users/blorenz/COSMOS/COSMOSData/scales.txt'
#Read in the scale of the lines 
#scale_df = ascii.read(scaledata).to_pandas()

#Location of the equivalent width data
ewdata = '/Users/blorenz/COSMOS/COSMOSData/lineew.txt'
#Read in the ew of the lines 
ew_df = ascii.read(ewdata).to_pandas()

#File with the error array
errdatapath = '/Users/blorenz/COSMOS/COSMOSData/errs.txt'
#Read in the scale of the lines 
err_df = ascii.read(errdatapath,data_start=1,header_start=0,format='csv').to_pandas()

#The location to store the scale and its stddev of each line
qualdatapath = '/Users/blorenz/COSMOS/COSMOSData/dataqual.txt'
#Read in the scale of the lines 
dataqual = ascii.read(qualdatapath).to_pandas()
d = {'True': True, 'False': False}

#File with the error array
errreddatapath = '/Users/blorenz/COSMOS/COSMOSData/errs_red.txt'

#Read the datafile:
fluxdata = ascii.read(fluxdatapath).to_pandas()

#Division function
def divz(X,Y):
        return X/np.where(Y,Y,Y+1)*np.not_equal(Y,0)


#Check if bpt correlates with stellar mass
#The location of the muzzin et al data:
mdatapath = '/Users/blorenz/COSMOS/muzzin_data/UVISTA_final_colors_sfrs_v4.1.dat'
#Read in the muzzin data
mdata = ascii.read(mdatapath).to_pandas()
mdata = mdata.rename(columns={'ID':'OBJID'})

fluxdata = pd.merge(fluxdata,mdata)


#Location of the reddening data
reddata = '/Users/blorenz/COSMOS/COSMOSData/reddenings.dat'
#Read in the ew of the lines 
red_df = ascii.read(reddata).to_pandas()
red_df = red_df.drop('zcc',axis=1)
fluxdata = pd.merge(fluxdata,red_df,on='fluxfile')
#Read in the medians by mass
avmeddata = '/Users/blorenz/COSMOS/COSMOSData/av_med_df.txt'
av_med_df = ascii.read(avmeddata).to_pandas()
for i in range(0,len(fluxdata)):
    if (fluxdata.iloc[i].av < 0):
        fluxdata.at[i,'dav1'] = 2
        fluxdata.at[i,'dav2'] = 2
        if (fluxdata.iloc[i].LMASS < 9.25):
            fluxdata.at[i,'av'] = av_med_df.mr1[0]
        elif (fluxdata.iloc[i].LMASS < 9.5):
            fluxdata.at[i,'av'] = av_med_df.mr2[0]
        elif (fluxdata.iloc[i].LMASS < 9.75):
            fluxdata.at[i,'av'] = av_med_df.mr3[0]
        else:
            fluxdata.at[i,'av'] = av_med_df.mr4[0]


#Define the reddening law
def Calzetti(wave,Av,Rv):
    waveum = wave*.0001
    if ((waveum >= 0.63) and (waveum <= 2.2)):
        k = 2.659*(-1.857+divz(1.040,waveum))+Rv
    if ((waveum < 0.63) and (waveum >= 0.12)):
        k = 2.659*(-2.156+divz(1.509,waveum)-divz(0.198,waveum**2)+divz(0.011,waveum**3))+Rv
    F_rat = 10**(divz(-0.4*k*Av,Rv))
    return F_rat


err_dfred=pd.DataFrame()
err_dfred['fluxfile'] = fluxdata['fluxfile']
alllines = ['3727','4102','4340','4861','4959','5007','6548','6548_fix','6563','6563_fix','6583','6583_fix']
#Set low objects to an upper limit
for line in alllines:
    for i in range(0,len(fluxdata)):
        if (fluxdata.iloc[i][line+'_flux'] == 0) and (dataqual.iloc[i][line+'_low']=='True'):
            fluxdata.at[i,line+'_flux'] = 5*err_df.iloc[i][line]
            #print ('New flux ' + str(fluxdata.iloc[i][line+'_flux']))

for line in alllines:
    if len(line) == 4: lineint = int(line)
    else: lineint = int(line[0:4])
    Av = fluxdata.av
    dAv1 = fluxdata.dav1
    dAv2 = fluxdata.dav2
    F_rat = Calzetti(lineint,Av,4.05)
    dFratu = Calzetti(lineint,Av+np.abs(dAv2),4.05)
    dFratd = Calzetti(lineint,Av-np.abs(dAv1),4.05)
    erru = (fluxdata[line+'_flux']/dFratu) - (fluxdata[line+'_flux']/F_rat)
    errd = (fluxdata[line+'_flux']/F_rat) - (fluxdata[line+'_flux']/dFratd)
    fluxdata[line+'_flux_red'] = fluxdata[line+'_flux']/F_rat
    err_dfred[line+'_u'] = np.sqrt(erru**2+(err_df[line]/F_rat)**2)
    err_dfred[line+'_d'] = np.sqrt(errd**2+(err_df[line]/F_rat)**2)


err_dfred.to_csv(errreddatapath,index=False)
fluxdata.to_csv(outpath,index=False)
