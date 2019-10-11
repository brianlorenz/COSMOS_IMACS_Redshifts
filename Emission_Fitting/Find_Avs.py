#Plots the Av magnitude due to the balmer decerment
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
#Read the datafile:
fluxdata = ascii.read(fluxdatapath).to_pandas()

#The location to store the scale and its stddev of each line
qualdatapath = '/Users/blorenz/COSMOS/COSMOSData/dataqual.txt'
#Read in the scale of the lines 
dataqual = ascii.read(qualdatapath).to_pandas()


#USIG ERRORS
#File with the error array
errdatapath = '/Users/blorenz/COSMOS/COSMOSData/errs.txt'
#Read in the scale of the lines 
err_df = ascii.read(errdatapath,data_start=1,header_start=0,format='csv').to_pandas()

'''
#BIWT DUP ERRORS
#File with the error array
errdatapath = '/Users/blorenz/COSMOS/COSMOSData/biwt_dup.txt'
#Read in the scale of the lines 
err_df = ascii.read(errdatapath,data_start=1,header_start=0,format='csv').to_pandas()
'''
#File to write the Av array to
dataout = '/Users/blorenz/COSMOS/COSMOSData/balmer_avs.txt'

#The location of the muzzin et al data:
mdatapath = '/Users/blorenz/COSMOS/muzzin_data/UVISTA_final_colors_sfrs_v4.1.dat'
#Read in the muzzin data
mdata = ascii.read(mdatapath).to_pandas()
mdata = mdata.rename(columns={'ID':'OBJID'})

fluxdata = pd.merge(fluxdata,mdata)

#Fontsizes for plotting
axisfont = 18
ticksize = 16
titlefont = 24
legendfont = 16
textfont = 16

#Division function
def divz(X,Y):
        return X/np.where(Y,Y,Y+1)*np.not_equal(Y,0)

lines = ['6563_fix','4861','4340']

#Fontsizes for plotting
axisfont = 18
ticksize = 16
titlefont = 24
legendfont = 16
textfont = 16

#Division function
def divz(X,Y):
        return X/np.where(Y,Y,Y+1)*np.not_equal(Y,0)


#Create the figure
fig,axarr = plt.subplots(3,3,figsize=(25,22))

#Plotting parameters
ms = 3
lw=0.5
mark='o'


#Set the Rv value - this one is taken for Calzetti
Rv = 4.05 #+-0.80

#Reddening law from Calzetti et al (2000)
def Calzetti_k(wave):
    waveum = wave*.0001
    if ((waveum >= 0.63) and (waveum <= 2.2)):
        k = 2.659*(-1.857+divz(1.040,waveum))+Rv
    if ((waveum < 0.63) and (waveum >= 0.12)):
        k = 2.659*(-2.156+divz(1.509,waveum)-divz(0.198,waveum**2)+divz(0.011,waveum**3))+Rv
    return k

#Finds the ratio and errors in that ratio of any two lines
def getAv(pd_df,err_df,L1,L2,bdec):
    #Calibrate the fluxes by dividing by scale
    calL1 = divz(pd_df[L1+'_flux'],pd_df[L1+'_scale'])
    calL2 = divz(pd_df[L2+'_flux'],pd_df[L2+'_scale'])
    #Find the ratio
    rat = divz(calL1,calL2)
    #Find the error in the ratio
    erat = np.sqrt((divz(1,calL2) * err_df[L1])**2 + (divz(-calL1,(calL2**2)) * err_df[L2])**2)
    #Get the integer of the line
    if len(L1)==8: iL1=int(L1[0:4])
    else: iL1 = int(L1)
    if len(L2)==8: iL2=int(L2[0:4])
    else: iL2 = int(L2)
    #Get the k value for each line
    L1k = Calzetti_k(iL1)
    L2k = Calzetti_k(iL2)
    #Compute the Av
    Av = divz(np.log10(rat/bdec),(0.4*((L2k/Rv)-(L1k/Rv))))
    #And its error
    eAv = divz((1/np.log(10))*divz((erat/bdec),rat),(0.4*((L2k/Rv)-(L1k/Rv))))
    return Av,eAv


d = {'True': True, 'False': False}

lines0 = ['6563_fix','4861']
lines1 = ['4861','4340']
lines2 = ['6563_fix','4340']

c = 0
Av_df = pd.DataFrame()
#Add the fluxfile so that this can later be merged with the main frame
Av_df['fluxfile'] = fluxdata['fluxfile']
Av_df['LMASS'] = fluxdata['LMASS']
for lines in [lines0,lines1,lines2]:
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
    if c==0:
        bdec=2.86
    elif c==1:
        bdec=2.137
    else:
        bdec=2.86*2.137
    Av,eAv = getAv(fluxdata,err_df,lines[0],lines[1],bdec)
    if c==0:
        Av_df['AvHaHb'] = Av
        Av_df['AvHaHberr'] = eAv
    elif c==1:
        Av_df['AvHbHg'] = Av
        Av_df['AvHbHgerr'] = eAv
    elif c==2:
        Av_df['AvHaHg'] = Av
        Av_df['AvHaHgerr'] = eAv
    c=c+1


#Get the average between the two good Avs and its error
Av_df['AvHa_avg'] = (Av_df['AvHaHb']+Av_df['AvHaHg'])/2
Av_df['AvHa_avgerr'] = (Av_df['AvHaHberr']+Av_df['AvHaHgerr'])/2


#Find the mass_weighted medians
mr1 = (fluxdata['LMASS']<9.25)
mr2 = np.logical_and(fluxdata['LMASS']>=9.25,fluxdata['LMASS']<9.5)
mr3 = np.logical_and(fluxdata['LMASS']>=9.5,fluxdata['LMASS']<9.75)
mr4 = (fluxdata['LMASS']>=9.75)
med1 = np.median(Av_df[np.logical_and(allgood,mr1)]['AvHaHb'])
med2 = np.median(Av_df[np.logical_and(allgood,mr2)]['AvHaHb'])
med3 = np.median(Av_df[np.logical_and(allgood,mr3)]['AvHaHb'])
med4 = np.median(Av_df[np.logical_and(allgood,mr4)]['AvHaHb'])
'''
#Linear fit for the medians
coeff = np.polyfit(fluxdata[goodidx]['LMASS'],AvHaHg[goodidx],1)
'''

Av_df = Av_df.replace(-np.inf,-99.99999999999)

d = {'True': True, 'False': False}

ulim=4 #Upper Av limit to consider good
#Number of stddevs away form median to be good
sig = 2

for i in range(0,len(fluxdata)):
    '''
    use 0 - 'AvHaHb'
    use 1 - 'AvHbHg'
    use 2 - 'AvHaHg'
    use 3 - 'AvMedian'
    use 4 - 'AvHa_avf'
    '''

    
    row = fluxdata.iloc[i]

    #Mass-weighted medians
    if (row['LMASS'] < 9.25): Av_df.at[i,'AvMedian'] = med1
    elif np.logical_and(row['LMASS'] >= 9.25,row['LMASS'] < 9.5,): Av_df.at[i,'AvMedian'] = med2
    elif np.logical_and(row['LMASS'] >= 9.5,row['LMASS'] < 9.75): Av_df.at[i,'AvMedian'] = med3
    elif (row['LMASS'] >= 9.75): Av_df.at[i,'AvMedian'] = med4

    '''
    #Linear fit for medians
    Av_df.at[i,'AvMedian'] = coeff[0]*row['LMASS']+coeff[1]
    '''
    Avrow = Av_df.iloc[i]

    if np.logical_or((Avrow['AvHaHb'] < 0),((Avrow['AvHaHb'] > ulim))): cHaHb = 10**80
    else: cHaHb = Avrow['AvHaHb']
    if np.logical_or((Avrow['AvHbHg'] < 0),((Avrow['AvHbHg'] > ulim))): cHbHg = 10**90
    else: cHbHg = Avrow['AvHbHg']
    if np.logical_or((Avrow['AvHaHg'] < 0),((Avrow['AvHaHg'] > ulim))): cHaHg = 10**100
    else: cHaHg = Avrow['AvHaHg']
        
    use = 3
    #Find out which lines are good. (Ha,Hb,Hg)
    l1g = dataqual['6563_fix_good'].map(d).iloc[i]
    l2g = dataqual['4861_good'].map(d).iloc[i]
    l3g = dataqual['4340_good'].map(d).iloc[i]
    
    goodlines = (l1g,l2g,l3g)
    #If only Ha and Hb are good, check those
    if goodlines == (1,1,0):
        Av_df.at[i,'AvHaHbok'] = Avrow['AvHaHb']
        if (cHaHb < 10): use=0
    #If only Ha and Hg are good, check those
    if goodlines == (1,0,1):
        if (cHaHg < 10): use=1
    #If all lines are good, 
    if goodlines == (1,1,1):
        #Compare HaHb and HaHg. if they are within each other's errors, average them
        diff = np.abs(cHaHb-cHaHg)
        err = Avrow['AvHaHberr']+Avrow['AvHaHgerr']
        #If they are within each other's errors, check if the error bars are large on one of the measurements
        if (diff < err):
            if  (divz(Avrow['AvHaHberr'],Avrow['AvHaHgerr']) < 0.5):
                if (cHaHb < 10): use=0
            elif (divz(Avrow['AvHaHberr'],Avrow['AvHaHgerr']) > 2):
                if (cHaHg < 10): use=1
            else:
                if (Avrow['AvHa_avg'] > 0): use = 4
        #If they are not close, use whichever is closest to the median
        else:
            diffHaHb = np.abs(cHaHb-Avrow['AvMedian'])
            diffHaHg = np.abs(cHaHg-Avrow['AvMedian'])
            arr = np.array([diffHaHb,diffHaHg])
            if (5 > arr[np.argmin(arr)]):
                use = np.argmin(arr)
        
                
        
            
            
            
    if use == 0: usestr = 'AvHaHb'
    elif use == 1: usestr = 'AvHaHg'
    elif use == 3: usestr = 'AvMedian'
    elif use == 4: usestr = 'AvHa_avg'
    Av_df.at[i,'useAv'] = usestr
     

#Write to csv
Av_df = Av_df.reindex(sorted(Av_df.columns), axis=1)
Av_df.to_csv(dataout,index=False)
