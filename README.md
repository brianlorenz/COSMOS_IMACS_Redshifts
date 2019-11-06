# COSMOS_IMACS_Redshifts
A pipeline to take processed big.fits files from the [IMACS spectrograph](http://www.lco.cl/telescopes-information/magellan/instruments/imacs) on the Baade telescope at Magellan and compute redshifts to galaxies


## Redshift_Measurement
Files located in Redshift_Measurement
Intended to run in this order:

ExtractSpectra.py - takes a 2-D big.fits file and outputs 1-D spectra files

FitStars.py - generates a telluric correction file from a given hot star in the mask

TCorr.py - applies the correction file to every object in the mask through cross correlation

CrossCor.py - GUI that finds the redshift to each object

## Dependencies
SDSS DR**2** Spectral Templates 23 - 27, available at [this URL](http://classic.sdss.org/dr2/algorithms/spectemplates/) in .FIT format

[Numpy](http://www.numpy.org/)

[SciPy](https://www.scipy.org/)

[Astropy](http://www.astropy.org/)

## Additional Processing
Files located in Data_Conversion
Pipeline for coverting data from redshifts into a more usable form and combining it with UVISTA and deimos data

run MergeData.py - this combines our data with the ULTRAVista data, selecting only those rows that we have observed. Outputs all_data_UVISTA.txt

run AddZ.py - this adds two missing columns to the above, the final redshift 'z_cc' and the object identifier 'OBJID'. outputs all_data_UVISTA_z.txt

run Combine_Hassinger.py - Finds the nearest object to each of ours in the Hassigner catalog and combines with ours, outputs all_data_hasinger.txt

run Add_Zerrs.py - Adds the missing error columsn from redshift measurement, outputs all_data_hasinger_dz.txt

## Emission fitting and line measurements
Files located in Emission_Fitting
Fits common emission lines in each galaxies, provides a confidence level for them, and computes reddening and star formation rates

run FitEmissionOnesig.py - Fits emission lines in all of the data and outputs quality graphs. Creates lineflux_tot.txt

run Drop_Bad_z.py - Drops the "galaxies" that are outside of the target redshift range from the dataframe (these are likely continanimants). Outputs lineflux.txt

run Find_MAD.py - computes the median absolute deviation in each line. Outputs linemad.txt

run Find_bad_low.py - Computes the cutoff for signal-to-noise for each line, classifying it as 'good', 'low', or 'bad.' Outputs dataqual.text

run Find_EW_errdf.py - Computes the equivvalent width of the Balmer absorption lines as well as their uncertainties. Outputs lineew.txt and errs.txt

Recieved the Avs of each of the lines from Dan. Cam also compute with Find_Avs.py, but it's lest Robust

run FindReddening.py - Computes the reddened flux of the measured emisison lines. Outputs lineflux_red.txt and err_df_red.txt

run Find_sfr.py - computes the star formation rate for every galaxy. Outputs sfrs.txt


