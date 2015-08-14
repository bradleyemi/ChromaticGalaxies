# ChromaticGalaxies

This is the code I used to generate a basis set of parametric galaxies from the AEGIS data for ChromaticRealGalaxy in Galsim (https://github.com/GalSim-developers/GalSim) for issue #551 (https://github.com/GalSim-developers/GalSim/issues/551). 

Requirements:

Python 2.7
Standard NumPy, SciPy distributions
PyPlot (optional)
SExtractor (http://www.astromatic.net/software/sextractor)
Astro Ascii Data (http://www.stecf.org/software/PYTHONtools/astroasciidata/)
Pyfits (http://www.stsci.edu/institute/software_hardware/pyfits)
TinyTim starfields corresponding to focus positions -10um to +5um saved as .fits files
A single file containing the centroids of the TinyTim stars (assumed to be the same for all starfields)- text file with x, y columns

Your data should include:

A list of CTI-corrected, multi-drizzled HST image .fits files with a WCS, and their corresponding inverse variance files. These should have drz.fits and wht.fits extensions respectively. 

The workflow is as follows:

1. Generate catalogs.

Run HST_sextractor.py on the data, initiate one GalaxyCatalog for each image/weight file. Run make_catalog on each GalaxyCatalog to make the SExtractor catalog and save the catalog to disk. Assemble these into an instance of GalaxyCatalogList, and then run add_focus to add the focus positions to the catalogs (you will need this to get the PSFs). Then run overlap.py on all your catalog files (with WCS) to delete duplicate objects. A manual masking file is optional, and the format is described within HST_sextractor.py. 

2. Get postage stamps.

Run postage_stamps.py on the list of catalogs with focus, images, and also input the filter in nm. This is accomplished by get_postage_stamps_all. You will need to change the path of the TinyTim files at the top of postage_stamps.py to correspond to where the TinyTim starfields are located, and also the path of the list of starfield centroids in the PSF function. You will also need to set the path in which the postage stamps are saved (this process will be made easier in the future). Run generate_masks to get mask files. 

3. Get parametric fits.

Right now, we are dependent on Claire Lackner's code to get parametric models. This requires IDL and the library IDLUtils. See Claire's page for more information (https://github.com/clackner2007/bdfit_lite).
