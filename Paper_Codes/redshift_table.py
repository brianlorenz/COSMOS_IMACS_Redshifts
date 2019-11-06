# Generates a table of objids, redshifts, uncertainties, and any flags

import sys
import os
import string
import numpy as np
import pandas as pd
from astropy.io import ascii
from tabulate import tabulate

# Location to save outputs
figout = '/Users/galaxies-air/COSMOS/COSMOS_IMACS_Redshifts/Paper_Images/'

# Read the data
red_data_loc = '/Users/galaxies-air/COSMOS/COSMOSData/all_c_hasinger_dz.txt'
red_data = ascii.read(red_data_loc).to_pandas()

# Filter stars
red_data = red_data[red_data['Star'] == 0]

# Clean up the OBJIDS to be readable
# Turning to int to get rid of .0, then adding leading zeros with str.zfill
objids = red_data.OBJID.astype('int')
objids_zero = [str(obj).zfill(6) for obj in objids]


# Set the columns and their order for the table
objids_data = [objids_zero, 'OBJID']
zcc_data = [red_data.z_cc, 'z']
dzhi_data = [red_data.dzhi, 'Uncertainty_Above']
dzlo_data = [red_data.dzlo, 'Uncertainty_Below']
unsure_data = [red_data.Unsure.astype(int), 'Flag_Unsure']
bad_data = [red_data.Bad.astype(int), 'Flag_Bad']
table_data = [objids_data, zcc_data, dzhi_data,
              dzlo_data, unsure_data, bad_data]

# Extract the data and corresponding names in order
columns = np.transpose([table_data[i][0] for i in range(len(table_data))])
column_names = [table_data[i][1] for i in range(len(table_data))]

# Generates the output table
out_table = pd.DataFrame(columns, columns=column_names)

# Sort on OBJID
out_table.sort_values('OBJID', inplace=True)


# Write the table in LaTeX format
f = open(figout+'redshifts_latex.txt', "w+")
f.write(tabulate(out_table, headers=(out_table.columns),
                 tablefmt="latex", showindex=False, floatfmt=".5f"))
f.close()

# Write as a csv
out_table.to_csv(figout+'redshifts.df', index=False)
