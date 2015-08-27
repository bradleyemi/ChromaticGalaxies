# ChromaticGalaxies

This is the code I used to generate a basis set of parametric galaxies from the AEGIS 
data for ChromaticRealGalaxy in Galsim (https://github.com/GalSim-developers/GalSim) for 
issue #551 (https://github.com/GalSim-developers/GalSim/issues/551). 

I can be contacted at bemi@stanford.edu to discuss the code or fix issues.

Software Requirements:

Python 2.7,

Standard NumPy, SciPy distributions,

PyPlot,

SExtractor (http://www.astromatic.net/software/sextractor),

Astro Ascii Data (http://www.stecf.org/software/PYTHONtools/astroasciidata/),

Pyfits (http://www.stsci.edu/institute/software_hardware/pyfits),

GalSim (https://github.com/GalSim-developers/GalSim),

For f606w and f814w filters: 

TinyTim starfield folders. This should be downloaded and placed into your path. 
For now you can download them off my Google Drive:
https://drive.google.com/folderview?id=0B5k2MfMTwAYdfndkREhuUG9aLWN2N3pQWUdxcTdtZlh3SXgxLUMyT2FuaVZvRWlNQmFINjQ&usp=sharing

They have to be unzipped before you can use them. You'll need to put in the correct
path in run.py to these folders.


Your data should include:

A list of CTI-corrected, multi-drizzled HST image .fits files with a WCS, and their 
corresponding inverse variance files. These should have drz.fits and wht.fits extensions 
respectively. 
The names of these filters IN ORDER should be saved in two separate text files. 
An optional file with coordinates of manual masks. See HST_Sextractor_new.py for details.

Description:

Makes a clean sample of postage stamps of galaxies and stores them in folders
compatible with Claire Lackner's parametric fitting code. 
Specifically it:

Runs SExtractor on the data
Cleans for blended objects using the hot-cold method (Leauthaud et. al. 2007)
Classifies stars and galaxies
Deletes objects on the noisy border
Automatic detection of star diffraction spikes, removal of objects affected
Checks for overlapping objects not caught by the hot-cold method
Finds the PSF by estimating the focus position from TT star fields
Generates postage-stamps and PSFs from simulated stars from the TT fields

Usage:

You need to have all of the files in this repository in the same directory. But all you
need to modify is the top of the script run.py. 

Then all you need to do is:
python run.py

For developers, here is the dependency tree of the modules:

run.py: requires assoc_catalogs, HST_sextractor_new, overlap, and postage_stamps. 
(will later require generate_masks, for now no masks are generated). 

HST_sextractor_new.py: requires focus_positions and cleanutils. 






