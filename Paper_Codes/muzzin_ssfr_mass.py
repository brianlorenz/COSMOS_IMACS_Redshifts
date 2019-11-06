# Plots the ssfr vs mass of the galaxies using UltraVISTA data.

import sys
import os
import string
import numpy as np
import pandas as pd
from astropy.io import ascii
from tabulate import tabulate
import matplotlib.pyplot as plt

# Location to store figure output
figout = '/Users/galaxies-air/COSMOS/COSMOS_IMACS_Redshifts/Paper_Images/'

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


# Fontsizes for plotting
axisfont = 24
ticksize = 18
ticks = 8
titlefont = 24
legendfont = 16
textfont = 16

# Figure setup
fig, ax = plt.subplots(figsize=(8, 7))

mark = '.'
xlab = 'log(Stellar Mass) (M$_\odot$)'
ylab = 'sSFR (yr$^{-1}$)'
xlim = (8.95, 10.05)
ylim = (-12, -8)

# Plot all the objects in red
plt.scatter(muzzin_mass.iloc[unique_indicies], muzzin_ssfr.iloc[unique_indicies], color='orange',
            marker=mark, label='Not measured')
# Overplot the objects that are good, which will cover the bad ones
plt.scatter(muzzin_mass.iloc[good_objects_indices], muzzin_ssfr.iloc[good_objects_indices], color='cornflowerblue',
            marker=mark, label='Measured')

# Set the axis labels
ax.set_xlabel(xlab, fontsize=axisfont)
ax.set_ylabel(ylab, fontsize=axisfont)

# Set the limits
ax.set_xlim(xlim)
ax.set_ylim(ylim)

# Set the tick size
ax.tick_params(labelsize=ticksize, size=ticks)
ax.legend(fontsize=axisfont-4)

# Save the figure
fig.savefig(figout+'ssfr_mass.pdf')
plt.close('all')
