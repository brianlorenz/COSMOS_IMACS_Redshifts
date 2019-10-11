#Converts the verb files to .txt (eg j7_verb.txt to j7.txt)

###Usage - run verb_to_txt.py 'a6'
#this will convert verb_a6.txt to a6.txt

import numpy as np
from astropy.io import ascii
import sys, os, string
import pandas as pd

letnum = sys.argv[1]

#Location of verb_xx.txt
verbloc = '/Users/blorenz/COSMOS/COSMOSData/corFitsFileOut/verb_' + letnum + '.txt'

ourdata = ascii.read(verbloc, delimiter=' ').to_pandas()

#Search each row for the temp number, find the redshift associated with that number and append it to a list
z_cc = []
dzhi = []
dzlo = []
ccmax = []
chi2 = []
rchi2 = []
for i in range(0,len(ourdata)):
    temp =  ourdata.iloc[i].temp
    if temp == 23:
        z_cc.append(ourdata.iloc[i].z23)
        dzhi.append(ourdata.iloc[i].dzhi23)
        dzlo.append(ourdata.iloc[i].dzlo23)
        ccmax.append(ourdata.iloc[i].ccmax23)
        chi2.append(ourdata.iloc[i].chi223)
        rchi2.append(ourdata.iloc[i].rchi223)
    elif temp == 24:
        z_cc.append(ourdata.iloc[i].z24)
        dzhi.append(ourdata.iloc[i].dzhi24)
        dzlo.append(ourdata.iloc[i].dzlo24)
        ccmax.append(ourdata.iloc[i].ccmax24)
        chi2.append(ourdata.iloc[i].chi224)
        rchi2.append(ourdata.iloc[i].rchi224)
    elif temp == 25:
        z_cc.append(ourdata.iloc[i].z25)
        dzhi.append(ourdata.iloc[i].dzhi25)
        dzlo.append(ourdata.iloc[i].dzlo25)
        ccmax.append(ourdata.iloc[i].ccmax25)
        chi2.append(ourdata.iloc[i].chi225)
        rchi2.append(ourdata.iloc[i].rchi225)
    elif temp == 26:
        z_cc.append(ourdata.iloc[i].z26)
        dzhi.append(ourdata.iloc[i].dzhi26)
        dzlo.append(ourdata.iloc[i].dzlo26)
        ccmax.append(ourdata.iloc[i].ccmax26)
        chi2.append(ourdata.iloc[i].chi226)
        rchi2.append(ourdata.iloc[i].rchi226)
    elif temp == 27:
        z_cc.append(ourdata.iloc[i].z27)
        dzhi.append(ourdata.iloc[i].dzhi27)
        dzlo.append(ourdata.iloc[i].dzlo27)
        ccmax.append(ourdata.iloc[i].ccmax27)
        chi2.append(ourdata.iloc[i].chi227)
        rchi2.append(ourdata.iloc[i].rchi227)
    else:
        z_cc.append(None)
        dzhi.append(None)
        dzlo.append(None)
        ccmax.append(None)
        chi2.append(None)
        rchi2.append(None)

#Turn the list into a df and join it with the data
z_ccdf = pd.DataFrame({'z_cc':z_cc})
dzhidf = pd.DataFrame({'dzhi':dzhi})
dzlodf = pd.DataFrame({'dzlo':dzlo})
ccmaxdf = pd.DataFrame({'ccmax':ccmax})
chi2df = pd.DataFrame({'chi2':chi2})
rchi2df = pd.DataFrame({'rchi2':rchi2})
ourdata = ourdata.join(z_ccdf)
ourdata = ourdata.join(dzhidf)
ourdata = ourdata.join(dzlodf)
ourdata = ourdata.join(ccmaxdf)
ourdata = ourdata.join(chi2df)
ourdata = ourdata.join(rchi2df)



##IF getting errors, change row.Revisit, row.Note, and row.Unusable to row.Flag1,row.Flag2, and row.Flag3 respectively
f = open(verbloc.replace('verb_' + letnum + '.txt', letnum + '.txt'),'w+')
f.write('#OBJID  temp    z           dzhi        dzlo        ccmax     chi2    rchi2  eclip     S/N      Star Bad  Unsure    ImageName              Flag1   Flag2    Flag3      Confidence\n')
for idx,row in ourdata.iterrows():
    f.write(('%06d   %d    %.6f    %.6f   %.6f    %7.2f   %7.2f   %2.2f     %d      %2.2f      %d    %d     %d ' % (row.OBJID,row.temp,row.z_cc,row.dzhi,row.dzlo, row.ccmax, row.chi2,row.rchi2,row.eclip,row['S/N'],row.Star,row.Bad,row.Unsure)) + row.ImageName + ('      %d      %d      %d      %d' % (row.Revisit, row.Note, row.Unusable, row.Confidence)) + '\n')
f.close()


