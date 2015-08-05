'''
Script Name: focus_positions.py
Author: Bradley Emi
Date: July 2015

v. 1.0 (7/21/2015) script finds the mean focus positions of images given a set of SExtractor catalogs, the corresponding 
image files, and set of TT starfields (as GalSim images) from -10 to +5 um. It outputs them to a text file.

########### Description ##########

This script determines the mean focus positions of the HST given TinyTim starfields for 16 focus positions from -10 to 5 um. 
Stars must be checked manually for input. 

########## Input ##########

-list of images (.fits files)
-list of SExtractor catalogs (.txt files) with the following parameters:
NUMBER
X_IMAGE
Y_IMAGE
FLUX_RADIUS
IS_STAR (MU_CLASS following convention in Leauthaud et. al. 2007, this can be modified in the code below. Just need a binary star/galaxy classifier with stars=1 and non-stars=0)
SNR
-list of TT GalSim images. Can be imported from fits files as shown below, just change the filenames below. 
-one text file that gives the centroids of the TT starfields in two columns, x and y. Can be an ascii table, SExtractor catalog, or regular text file.

*Images and catalogs must be in the same order, i.e., image 1 has to correspond to catalog 1


########## Dependencies ##########

AstroAsciiData
PyFits/AstroPy
NumPy
MatPlotLib.PyPlot
GalSim
SciPy
Ds9

########## Usage ##########

If you don't have the TT images as galsim images:
Change the filenames below, the script will automatically convert them to galsim images in the list "tt_galsim_images".

Then run:

focus(catalogs, filenames, tt_galsim_images, tt_star_file, out_name, match_dist = 200., stamp_size = 6., nstars=20)

The script will write out the files to a text file called <out_name>. 

Adjustable parameters:
match_dist -- approximately 1/2 of the distance between TinyTim stars. The script will print "no star found" if this is too small, if this is the case, make it larger.
stamp_size -- in Gaussian standard deviations, the size of the postage stamp of each star.
nstars --- the number of input stars to review

When the script prompts you to accept or reject input stars, open a new terminal window, go to the directory where the script is running and type:

ds9 Star*.fits

to review the input stars. Type y or n to accept or reject stars in order. Once you are done, it will display a histogram of the focus positions, and
write out the focus positions to a text file. The stars used will be written out to files called <catalog>.stars. 
'''

import asciidata
import subprocess
import pyfits
import numpy as np
import matplotlib.pyplot as plt
import time
import galsim
from scipy.optimize import curve_fit

'''
Steps: 
1) Import all the TinyTim fields
2) Feed an image into GalSim
3) Pick the best S/N stars in the image
4) Get the corresponding TinyTim stars
5) Run pattern match, find focus position
6) Write focuses into a text file
'''

###### Modify the filenames below to get TT fits files into GalSim image format. ######

def get_tt_files(filter):
    if filter == 606:
         focusDict = {-1 : "/Users/bemi/JPL/F606W_TT/TinyTim_f-1.fits",
                   -2 : "/Users/bemi/JPL/F606W_TT/TinyTim_f-2.fits",
                   -3 : "/Users/bemi/JPL/F606W_TT/TinyTim_f-3.fits",
                   -4 : "/Users/bemi/JPL/F606W_TT/TinyTim_f-4.fits",
                   -5 : "/Users/bemi/JPL/F606W_TT/TinyTim_f-5.fits",
                   -6 : "/Users/bemi/JPL/F606W_TT/TinyTim_f-6.fits",
                   -7 : "/Users/bemi/JPL/F606W_TT/TinyTim_f-7.fits",
                   -8 : "/Users/bemi/JPL/F606W_TT/TinyTim_f-8.fits",
                   -9 : "/Users/bemi/JPL/F606W_TT/TinyTim_f-9.fits",
                   -10: "/Users/bemi/JPL/F606W_TT/TinyTim_f-10.fits",
                   0  : "/Users/bemi/JPL/F606W_TT/TinyTim_f0.fits",
                   1  : "/Users/bemi/JPL/F606W_TT/TinyTim_f1.fits",
                   2  : "/Users/bemi/JPL/F606W_TT/TinyTim_f2.fits",
                   3  : "/Users/bemi/JPL/F606W_TT/TinyTim_f3.fits",
                   4  : "/Users/bemi/JPL/F606W_TT/TinyTim_f4.fits",
                   5  : "/Users/bemi/JPL/F606W_TT/TinyTim_f5.fits"}
                   
    elif filter == 814:
        focusDict = {-1 : "/Users/bemi/JPL/F814W_TT/TinyTim_f-1.fits",
                     -2 : "/Users/bemi/JPL/F814W_TT/TinyTim_f-2.fits",
                     -3 : "/Users/bemi/JPL/F814W_TT/TinyTim_f-3.fits",
                     -4 : "/Users/bemi/JPL/F814W_TT/TinyTim_f-4.fits",
                     -5 : "/Users/bemi/JPL/F814W_TT/TinyTim_f-5.fits",
                     -6 : "/Users/bemi/JPL/F814W_TT/TinyTim_f-6.fits",
                     -7 : "/Users/bemi/JPL/F814W_TT/TinyTim_f-7.fits",
                     -8 : "/Users/bemi/JPL/F814W_TT/TinyTim_f-8.fits",
                     -9 : "/Users/bemi/JPL/F814W_TT/TinyTim_f-9.fits",
                    -10: "/Users/bemi/JPL/F814W_TT/TinyTim_f-10.fits",
                      0  : "/Users/bemi/JPL/F814W_TT/TinyTim_f0.fits",
                      1  : "/Users/bemi/JPL/F814W_TT/TinyTim_f1.fits",
                      2  : "/Users/bemi/JPL/F814W_TT/TinyTim_f2.fits",
                      3  : "/Users/bemi/JPL/F814W_TT/TinyTim_f3.fits",
                      4  : "/Users/bemi/JPL/F814W_TT/TinyTim_f4.fits",
                      5  : "/Users/bemi/JPL/F814W_TT/TinyTim_f5.fits"}
    else:
        raise ValueError("No data for input filter.")
    tt_list = []
    for i in range(-10,6):
        tt_list.append(focus_dict[i])
    i = 0
    tt_galsim_images = []
    for image in tt_list:
        print "importing image", i
        f = pyfits.open(image)
        image_data = f[0].data 
        img = galsim.Image(image_data)
        tt_galsim_images.append(img)
        f.close()
        i += 1
    return tt_galsim_images
    
def get_star_file(filter)
    if filter == 606:
       return "Users/bemi/JPL/606_stars.txt"
    if filter == 814:
       return "Users/bemi/JPL/F814W_TT/TinyTim_f-1.stars.dat"

###### End .fits to GalSim import ####

output_params = ["X_IMAGE",
"Y_IMAGE"]

def run_sextractor_tt(file,output_params,out_name,clean=True): 
    #Create params_file and write out to a file
    param_ascii = asciidata.create(1,len(output_params))
    row_counter = 0
    for param in output_params:
        param_ascii[0][row_counter] = param 
        row_counter += 1
    param_fname = out_name + ".param"
    param_ascii.writeto(param_fname)
    
    #Create config newfiles[i] and write out to a file
    config_ascii = asciidata.create(2,5)
        
    #File-Specific Configurations
    config_ascii[0][0] = 'CATALOG_NAME'
    config_ascii[1][0] = out_name
    config_ascii[0][1] = 'PARAMETERS_NAME'
    config_ascii[1][1] = out_name + ".param"
    config_ascii[0][2] = 'BACK_TYPE'
    config_ascii[1][2] = 'MANUAL'
    config_ascii[0][3] = 'BACK_VALUE'
    config_ascii[1][3] = 1e-8
    config_ascii[0][4] = 'DETECT_THRESH'
    config_ascii[1][4] = 1e-10
    config_fname = out_name + ".config"
    config_ascii.writeto(config_fname)
        
    #Run sextractor and get the catalog
    subprocess.call(["sex", file , "-c", config_fname])
    
    #Optional Clean
    subprocess.call(["rm", config_fname])
    subprocess.call(["rm", param_fname])

'''
for key in tt_606:
    run_sextractor_tt(tt_606[key], output_params, "TinyTim_f" + str(key) + ".stars.dat", clean=True)
'''

def renumber(catalog):
    cat = asciidata.open(catalog)
    for i in range(cat.nrows):
        cat['NUMBER'][i] = i
    cat.writeto(catalog)
    
def select_good_stars(catalog, out_name, nstars=10):
    renumber(catalog)
    cat = asciidata.open(catalog)
    numbers = []
    snrs = []
    for i in range(cat.nrows):
        if cat['IS_STAR'][i] == 1: # and cat['MAG_AUTO'][i]+25 < 25.0:
            numbers.append(cat['NUMBER'][i])
            snrs.append(cat['SNR'][i])
    snr_arr = np.asarray(snrs)
    nstars = len(snrs)
    max_indexes = np.argsort(snr_arr)[-1*nstars:]
    keep = []
    for i in range(nstars):
        keep.append(numbers[max_indexes[i]])
    star_table = asciidata.create(3,nstars)
    for i in range(nstars):
        x = cat['X_IMAGE'][keep[i]]
        y = cat['Y_IMAGE'][keep[i]]
        r = cat['FLUX_RADIUS'][keep[i]]
        star_table[0][i] = x
        star_table[1][i] = y
        star_table[2][i] = r
    star_table.writeto(out_name)

def get_subImages(image_star_table, image, stamp_size = 6.):
    subImages = []
    table = asciidata.open(image_star_table)
    #out_table = asciidata.create(4,table.nrows)
    f = pyfits.open(image)
    image_data = f[0].data
    img = galsim.Image(image_data)
    for i in range(table.nrows):
        x0 = table[0][i]
        y0 = table[1][i]
        r = table[2][i]
        L = stamp_size * r
        b = galsim.BoundsI(int(x0-0.5*L), int(x0+0.5*L), int(y0-0.5*L), int(y0+0.5*L))
        sub = img.subImage(b)
        subImages.append(sub)
    return subImages

#Gets the centroid of a close TT star
def match_to_tt(image_star_table, tt_star_data_file, dist=200.):
    centroid_table = asciidata.open(tt_star_data_file)
    star_table = asciidata.open(image_star_table)
    out_table = asciidata.create(6, star_table.nrows)
    tt_centroids = []
    for i in range(star_table.nrows):
        x0 = star_table[0][i]
        y0 = star_table[1][i]
        r = star_table[2][i]
        for j in range(centroid_table.nrows):
            x = centroid_table[0][j]
            y = centroid_table[1][j]
            if abs(x0-x) < dist and abs(y0-y) < dist:
                tt_centroids.append((x,y,r))
                break
            if j == centroid_table.nrows-1:
                print "no star found"
    return tt_centroids

def get_tt_subImages(centroid, tt_images, stamp_size = 6.0):
    tt_subImages = []
    for img in tt_images:
        x0 = centroid[0]
        y0 = centroid[1]
        r = centroid[2]
        L = stamp_size * r
        b = galsim.BoundsI(int(x0-0.5*L), int(x0+0.5*L), int(y0-0.5*L), int(y0+0.5*L))
        sub = img.subImage(b)
        tt_subImages.append(sub)
    return tt_subImages

def get_cost(subImage_moments, tt_moment_lists):
    #Given a list of star moments, and a list of 15 possible TT moments corresponding to each one, find focus position.
    #Calculate a cost function for each focus position.
    costFunction = np.zeros(len(tt_moment_lists[0]))
    for i in range(len(costFunction)):
        for j in range(len(subImage_moments)):
            e1 = subImage_moments[j].observed_shape.getE1()
            e2 = subImage_moments[j].observed_shape.getE2()
            moment = (tt_moment_lists[j])[i]
            costFunction[i] += (moment.observed_shape.getE1() - e1)**2 + (moment.observed_shape.getE2() - e2)**2
    return costFunction
    #return np.argmin(costFunction)

def poly2(x, a, b, c):
    return a*x**2 + b*x + c

def find_focus_position(focus_list, cost_list, plot=True):
    popt, pcov = curve_fit(poly2, focus_list, cost_list)
    b = popt[1]
    a = popt[0]
    perr = np.sqrt(np.diag(pcov))
    sigma_a = perr[0]
    sigma_b = perr[1]
    if plot == True:
        plt.xlabel("Focus position (um)")
        plt.ylabel("Cost")
        plt.title("Error associated with TinyTim focus positions for 10 high S/N stars")
        plt.scatter(focus_list, cost_list)
        plt.plot(focus_list, poly2(focus_list, a, b, popt[2]))
        plt.show()
    return -b/(2*a), (1/(2*a**2))*np.sqrt(a**2*sigma_b**2 + b**2*sigma_a**2)   
    
def getMoments(image_star_table, image, tt_galsim_images, tt_star_file, match_dist = 200., stamp_size = 6., plot=True):  
    #One subimage for each star
    subImages = get_subImages(image_star_table, image, stamp_size = stamp_size)
    #One centroid for each star
    tt_centroids = match_to_tt(image_star_table, tt_star_file, dist=match_dist)
    #One list of tt_subimages for each star
    tt_subImages = []
    #One moment for each star
    subImage_moments = []
    keep = []
    i = 10
    for im in subImages:
        im.write("Star" + str(i) + ".fits")
        i += 1
    i = 10
    for im in subImages:
        print "View image", "Star" + str(i) + ".fits"
        accept = ""
        accept = raw_input("Accept? (y/n) -->")
        if accept == "y":
            print "Input accepted."
            try:
                subImage_moments.append(galsim.hsm.FindAdaptiveMom(im))
                keep.append(i)
            except:
                print "Adaptive moment did not converge, skipping..."
        else:
            print "Input rejected."
        i += 1
    j = 10
    print "Keeping stars", keep
    for centroid in tt_centroids:
        if j in keep:
            tt_subImages.append(get_tt_subImages(centroid, tt_galsim_images, stamp_size = stamp_size))
        j += 1
    #One moment for each tt_subimage
    tt_moment_lists = []
    for tt_focus_position in tt_subImages:
        tt_moments = []
        for tt_subImage in tt_focus_position:
            tt_moments.append(galsim.hsm.FindAdaptiveMom(tt_subImage))
        tt_moment_lists.append(tt_moments)
    cost = np.asarray(get_cost(subImage_moments, tt_moment_lists))
    focus, focus_err = find_focus_position(np.asarray(range(-10,6)), cost, plot=plot)
    subprocess.call(["rm", "Users/bemi/JPL/Star*.fits"])
    return focus, len(keep)

def focus(catalogs, filenames, tt_galsim_images, tt_star_file, out_name, match_dist = 200., stamp_size = 6., nstars=20, plot=False, generate_new_star_files=True, histogram = True):
    if generate_new_star_files:
        n = 0
        for catalog in catalogs:
            select_good_stars(catalog, catalog+".stars", nstars=nstars)
            print "generating star file", n
            n += 1
    foci = []
    err = []
    for i in range(len(filenames)):
        focus, focus_nstars = getMoments(catalogs[i]+".stars", filenames[i], tt_galsim_images, tt_star_file, match_dist=match_dist, stamp_size=stamp_size, plot=plot)
        print "Focus is", focus, "using", focus_nstars, "stars for calibration."
        out = open(out_name, "a")
        out.write(filenames[i] + " ")
        out.write(str(focus) + " ")
        out.write(str(focus_nstars) + "\n")
        out.close()
        foci.append(focus)
        err.append(focus_nstars)
    if histogram:
        plt.hist(foci, bins=20)
        plt.xlabel("focus position (um)")
        plt.ylabel("frequency")
        plt.title("Focus positions")
        plt.show()

def label_catalogs(focus_text_file):
    f = open(focus_text_file)
    for line in f.readlines():
         split = line.split()
         filename = split[0]
         focus = np.float32(split[1])
         catalog = asciidata.open(filename + ".cat")
         for i in range(catalog.nrows):
              catalog['FILENAME'][i] = filename
              catalog['FOCUS'][i] = focus
         catalog['FILENAME'].set_colcomment('Original name of image file for object')
         catalog['FOCUS'].set_colcomment('Focus position in um')
         catalog.writeto(filename + ".focus.cat")

#label_catalogs("606_focus_positions.txt")
#label_catalogs("814_focus_positions.txt")
