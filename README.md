# COSMOS_IMACS_Redshifts
A pipeline to take processed big.fits files from the IMACS spectrograph on the Baade telescope at Magellan and compute redshifts to galaxies


Intended to run in this order:

ExtractSpectra.py - takes a 2-D big.fits file and outputs 1-D spectra files

FitStars.py - generates a telluric correction file from a given hot star in the mask

TCorr.py - applies the correction file to every object in the mask through corss correlation

CrossCor.py - GUI that finds the redshift to each object

            - Requires spectral templates 23-27 from Sloan Digital Sky Survey DR2
