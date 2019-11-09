# Computes the ssfr and mass of all target objects

'''

'''

import sys
import os
import string
import numpy as np
import pandas as pd
from astropy.io import ascii

# Read the redshift data
red_data_loc = '/Users/galaxies-air/COSMOS/COSMOSData/all_c_hasinger.txt'
red_data = ascii.read(red_data_loc).to_pandas()

# Read the muzzin data
muzzin_data_loc = '/Users/galaxies-air/COSMOS/Muzzin_Data/UVISTA_final_colors_sfrs_v4.1.dat'
muzzin_data = ascii.read(muzzin_data_loc).to_pandas()

# Preserving order and number of objects in red_data (our measurements)
merged_data = red_data.merge(
    muzzin_data, how='left', left_on='id', right_on='ID')

''' Duplicate Filtering '''


def filter_df(df):
    # Input a pandas dataframe (df), gives back indicies of rows that have good measurements
    return np.logical_and.reduce((df['Bad'] == 0, df['Unsure'] == 0, df['Star'] == 0))


# First we need to isolate the duplicates into a dataframe
duplicate_indicies = merged_data.duplicated(subset='id', keep=False)
duplicate_df = merged_data[duplicate_indicies]

# Then, find all objects that have at least one good measurement, and remove any duplicates from that
duplicates_df_filt = duplicate_df[filter_df(duplicate_df)]
# This keeps only the first measurement of each object - possibly chagne to be the one with lower errors later
duplicates_df_filt = duplicates_df_filt[np.logical_not(duplicates_df_filt.duplicated(
    subset='id'))]

# Find which single objects have good measurements
single_objs_df = merged_data[np.logical_not(duplicate_indicies)]
single_objs_df_filt = single_objs_df[filter_df(single_objs_df)]


# We are now left with a df of object IDs that were duplicates, but now have at least one good measurement
# Concatenate this df with the objects that were NOT duplicates, to get a full count of objects that have good redshifts
good_objects_indices = pd.concat([single_objs_df_filt,
                                  duplicates_df_filt]).index

# Count the total number of objects
no_stars = merged_data[merged_data['Star'] == 0]
unique_objs = np.logical_not(no_stars.duplicated(subset='id'))
unique_indicies = no_stars[unique_objs].index
total_measured = len(good_objects_indices)
total_objects = len(unique_indicies)
print('Measured objects: ' + str(total_measured))
print('Total objects: ' + str(total_objects))
print('Fraction measured: ' + str(total_measured/total_objects))


muzzin_mass = merged_data['LMASS']
muzzin_ssfr = np.log10(merged_data['SFR_tot'] / (10**muzzin_mass))

# Putting only the relevant data into a frame
all_objs = pd.DataFrame(
    np.transpose([merged_data['OBJID'], muzzin_mass, muzzin_ssfr, np.zeros(len(merged_data))]), columns=['OBJID', 'LMASS', 'sSFR', 'Measured'])

# Setting the objects that had a measurement to 1 before we cahnge the indices
all_objs['Measured'].iloc[good_objects_indices] = 1.

# Selecting only the unique objects
unique_objs = all_objs.iloc[unique_indicies]
