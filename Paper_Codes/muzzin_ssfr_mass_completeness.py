# Plots the ssfr vs mass of the galaxies using UltraVISTA data.

'''
Usage: run muzzing_ssfr_mass_completeness dmass dssfr
e.g. run muzzing_ssfr_mass_completeness 0.2 0.5
dmass - size of log(mass) bin - default 0.1
dssfr - size of ssfr bin - default 0.5
'''

# muzzin_mass, muzzin_ssfr
from find_muzzin_mass_ssfr import unique_objs
import sys
import os
import string
import numpy as np
import pandas as pd
from astropy.io import ascii
from itertools import product
from tabulate import tabulate
import matplotlib.pyplot as plt

'''
unique_obs = dataframe of unique objects in our dataset, already filtered for stars. Has columns:
OBJID - Object ID
LMASS - Log Stellar Mass
sSFR - Specific Star Formation Rate
Measured - 1 if we have a measured redshift, 0 if not
'''

try:
    dmass = float(sys.argv[1])
    dssfr = float(sys.argv[2])
except:
    print('No input detected - using default dmass = 0.1, dssfr = 0.5')
    dmass, dssfr = (0.1, 0.5)

# Location to store figure output
figout = '/Users/galaxies-air/COSMOS/Paper_Images/'


def get_bins(dmass, dssfr):
    # Makes bins in mass and ssfr of size dmass and dssfr
    mass_bins = np.arange(9, 10+dmass, dmass)
    ssfr_bins = np.arange(-12, -8+dssfr, dssfr)
    return mass_bins, ssfr_bins

# WORKING HERE - NEED TO SEPARATE ALL AND UNIQUE WHEN BINNING, THEN DIVIDE FOR COMPLETENESS


def compute_completeness(mass_bins, ssfr_bins, unique_objs, measured_objs):
    # Computes the completeness in each box
    unique_binned = np.histogram2d(unique_objs['LMASS'], unique_objs['sSFR'], bins=[
        mass_bins, ssfr_bins])
    measured_binned = np.histogram2d(measured_objs['LMASS'], measured_objs['sSFR'], bins=[
        mass_bins, ssfr_bins])
    # Need to transpose it to get the axes to line up
    completeness_arr = np.transpose(measured_binned[0]/unique_binned[0])
    # The flip to account for ssfr being backwards. Compared by eye to the ssfr_mass plot
    completeness_arr = np.flip(completeness_arr, axis=0)
    return completeness_arr


measured_objs = unique_objs[unique_objs['Measured'] == 1]

mass_bins, ssfr_bins = get_bins(dmass, dssfr)
completeness_arr = compute_completeness(
    mass_bins, ssfr_bins, unique_objs, measured_objs)

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
xlim = (mass_bins[0], mass_bins[-1])
ylim = (ssfr_bins[0], ssfr_bins[-1])


# Extent sets where to do the imshow in data coords
extent = [mass_bins[0], mass_bins[-1], ssfr_bins[0], ssfr_bins[-1]]
# aspect stretches the pixes to make them boxes
image = ax.imshow(completeness_arr, extent=extent,
                  interpolation='nearest', aspect='auto')
plt.colorbar(image)

'''
# Plot all the objects in red
plt.scatter(unique_objs['LMASS'], unique_objs['sSFR'], color='black',
            marker=mark, label='Not measured')
# Overplot the objects that are good, which will cover the bad ones
plt.scatter(measured_objs['LMASS'], measured_objs['sSFR'], color='white',
            marker=mark, label='Measured')
'''


for (x, y) in product(np.arange(len(mass_bins)-1), np.arange(len(ssfr_bins)-1)):
    try:
        ax.text(mass_bins[x]+dmass/2, ssfr_bins[::-1][y+1]+dssfr/2, f'{int(round(completeness_arr[y,x],2)*100)}%',
                horizontalalignment='center', verticalalignment='center')
    except:
        pass


# Set the axis labels
ax.set_xlabel(xlab, fontsize=axisfont)
ax.set_ylabel(ylab, fontsize=axisfont)

# Set the limits
ax.set_xlim(xlim)
ax.set_ylim(ylim)

# Set the tick size
ax.tick_params(labelsize=ticksize, size=ticks)
# ax.legend(fontsize=axisfont-6)

plt.tight_layout()

# Save the figure
fig.savefig(figout+'ssfr_mass_density.pdf')
plt.close('all')
