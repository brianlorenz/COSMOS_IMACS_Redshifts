# COSMOS_IMACS_Redshifts
A pipeline to take processed big.fits files from the [IMACS spectrograph](http://www.lco.cl/telescopes-information/magellan/instruments/imacs) on the Baade telescope at Magellan and compute redshifts to galaxies


Intended to run in this order:

ExtractSpectra.py - takes a 2-D big.fits file and outputs 1-D spectra files

FitStars.py - generates a telluric correction file from a given hot star in the mask

TCorr.py - applies the correction file to every object in the mask through corss correlation

CrossCor.py - GUI that finds the redshift to each object


## Dependencies

SDSS DR**2** Spectral Templates 23 - 27, available at [this URL](http://classic.sdss.org/dr2/algorithms/spectemplates/) in .FIT format

[Numpy](http://www.numpy.org/)

[SciPy](https://www.scipy.org/)

[Astropy](http://www.astropy.org/)
