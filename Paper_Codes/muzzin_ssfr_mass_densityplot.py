# Plots the ssfr vs mass of the galaxies using UltraVISTA data.

# muzzin_mass, muzzin_ssfr
from find_muzzin_mass_ssfr import muzzin_mass, muzzin_ssfr, unique_indicies, good_objects_indices
import sys
import os
import string
import numpy as np
import pandas as pd
from astropy.io import ascii
from tabulate import tabulate
import matplotlib.pyplot as plt

# Location to store figure output
figout = '/Users/galaxies-air/COSMOS/Paper_Images/'


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
ax.legend(fontsize=axisfont-6)

plt.tight_layout()

# Save the figure
fig.savefig(figout+'ssfr_mass.pdf')
plt.close('all')
