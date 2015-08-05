'''
We used methods from this script to get an estimate on the focus uncertainty. It is mostly identical to focus_positions.py
'''

import asciidata
import subprocess
import pyfits
import numpy as np
import matplotlib.pyplot as plt
import time
import random
import galsim
from scipy.optimize import curve_fit

###### Modify the filenames below to get TT fits files into GalSim image format. ######

tt_606 = {-1 : "/Users/bemi/JPL/F606W_TT/TinyTim_f-1.fits",
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

tt_606_list = []
for i in range(-10,6):
    tt_606_list.append(tt_606[i])

tt_814 = {-1 : "/Users/bemi/JPL/F814W_TT/TinyTim_f-1.fits",
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

tt_814_list = []
for i in range(-10,6):
    tt_814_list.append(tt_814[i])

i = 0
tt_galsim_images = []
for image in tt_814_list:
    print "importing image", i
    f = pyfits.open(image)
    image_data = f[0].data 
    img = galsim.Image(image_data)
    tt_galsim_images.append(img)
    f.close()
    i += 1


output_params = ["NUMBER",
"X_IMAGE",
"Y_IMAGE",
"A_IMAGE",
"B_IMAGE",
"ALPHA_SKY",
"DELTA_SKY",
"XMIN_IMAGE",
"XMAX_IMAGE",
"YMIN_IMAGE",
"YMAX_IMAGE",
"FLAGS",
"MU_MAX",
"MAG_AUTO",
"CLASS_STAR",
"FLUX_RADIUS",
"FLUX_AUTO",
"FLUXERR_AUTO"]

#Set the global config files

bright_config_dict = { 'DETECT_MINAREA' : 140 ,
'DETECT_THRESH' : 2.2 ,
'DEBLEND_NTHRESH' : 64 ,
'DEBLEND_MINCONT' : 0.04 ,
'CLEAN_PARAM' : 1.0 ,
'BACK_SIZE' : 400 ,
'BACK_FILTERSIZE' : 5 ,
'BACKPHOTO_TYPE' : "LOCAL" ,
'BACKPHOTO_THICK' : 200,
'PIXEL_SCALE' : 0.03}

def run_sextractor(file,use_dict,output_params,out_name,clean=True): 
    #Make sure the file opens
    try:
        hdulist = pyfits.open(file)
        hdulist.close()
    except IOError:
        print "IO Error for filename: ", file
    #Create params_file and write out to a file
    param_ascii = asciidata.create(1,len(output_params))
    row_counter = 0
    for param in output_params:
        param_ascii[0][row_counter] = param
        row_counter += 1
    param_fname = out_name + ".param"
    param_ascii.writeto(param_fname)
    
    #Create config newfiles[i] and write out to a file
    config_ascii = asciidata.create(2,2+len(use_dict))
        
    #File-Specific Configurations
    config_ascii[0][0] = 'CATALOG_NAME'
    config_ascii[1][0] = out_name + ".cat"
    config_ascii[0][1] = 'PARAMETERS_NAME'
    config_ascii[1][1] = out_name + ".param"
    row_counter = 2
    for key, value in use_dict.iteritems():
        config_ascii[0][row_counter] = key
        config_ascii[1][row_counter] = value
        row_counter += 1
    config_fname = out_name + ".config"
    config_ascii.writeto(config_fname)                
    #Run sextractor and get the catalog
    subprocess.call(["sex", file , "-c", config_fname])
    #Optional Clean
    if clean:
        subprocess.call(["rm", config_fname])
        subprocess.call(["rm", param_fname])
    
#run_sextractor("j6ll01ycq_drz.fits", bright_config_dict, output_params, "47Tuc_f606w", clean=True)
#run_sextractor("j6ll01ykq_drz.fits", bright_config_dict, output_params, "47Tuc_f814w", clean=True)

def mu_mag(catalog):
    mu = []
    mag = []
    cat = asciidata.open(catalog)
    for i in range(cat.nrows):
        mu.append(cat['MU_MAX'][i])
        mag.append(cat['MAG_AUTO'][i])
    mu_array = np.asarray(mu)
    mag_array = np.asarray(mag)
    plt.scatter(mag_array, mu_array)
    plt.show()

#mu_mag("47Tuc_f814w.cat")

def delete_mu_cutoff(catalog, out_name, mu_cutoff=-14.6):
    new_cat = open(out_name, "w")
    cat = open(catalog, "r")
    for line in cat.readlines():
        split = line.split()
        if split[0] == "#" or np.float32(split[12]) > mu_cutoff:
            new_cat.write(line)
    cat.close()
    new_cat.close()

#delete_mu_cutoff("47Tuc_f814w.cat","47Tuc_f814w_unsat.cat")

def delete_flagged(catalog, out_name):
    new_cat = open(out_name, "w")
    cat = open(catalog, "r")
    for line in cat.readlines():
        split = line.split()
        if split[0] == "#" or np.float32(split[11]) == 0:
            new_cat.write(line)
    cat.close()
    new_cat.close()

#delete_flagged("47Tuc_f814w_unsat.cat", "47Tuc_f814w_filter.cat")

def add_used(catalog):
    cat = asciidata.open(catalog)
    for i in range(cat.nrows):
        cat['NUMBER'][i] = i
        cat['USED'][i] = 0
    cat.writeto(catalog)

#add_used("47Tuc_f814w_filter.cat")

def select_good_stars(catalog, out_name, nstars=10):
    cat = asciidata.open(catalog)
    numbers = []
    for i in range(cat.nrows):
        numbers.append(cat['NUMBER'][i])
    random.shuffle(numbers)
    keep = []
    for i in range(nstars):
        keep.append(numbers[i])
    star_table = asciidata.create(3,nstars)
    for i in range(nstars):
        x = cat['X_IMAGE'][keep[i]]
        y = cat['Y_IMAGE'][keep[i]]
        r = cat['FLUX_RADIUS'][keep[i]]
        star_table[0][i] = x
        star_table[1][i] = y
        star_table[2][i] = r
    star_table.writeto(out_name)

#select_good_stars("47Tuc_f814w_filter.cat", "47Tuc_f814w_filter.cat.stars", nstars=60)

def fix(n):
    if n < 1:
        return 1
    else:
        return n

def get_subImages(image_star_table, image, stamp_size = 6.):
    subImages = []
    table = asciidata.open(image_star_table)
    #out_table = asciidata.create(4,table.nrows)
    f = pyfits.open(image)
    image_data = f[1].data
    img = galsim.Image(image_data)
    for i in range(table.nrows):
        x0 = table[0][i]
        y0 = table[1][i]
        r = table[2][i]
        L = stamp_size * r
        b = galsim.BoundsI(fix(int(x0-0.5*L)), fix(int(x0+0.5*L)), fix(int(y0-0.5*L)), fix(int(y0+0.5*L)))
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
        x0 = star_table[0][i]*(7500./4210.)
        y0 = star_table[1][i]*(7500./4242.)
        r = star_table[2][i]
        for j in range(centroid_table.nrows):
            x = centroid_table[0][j]
            y = centroid_table[1][j]
            if abs(x0-x) < dist and abs(y0-y) < dist:
                tt_centroids.append((x,y,r))
                break
            if j == centroid_table.nrows-1:
                print "no star found at location", x0, y0
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
    
def get_cost_multiple(subImage_moments, tt_moment_lists, star_n_list):
    out_costs = []
    for nstars in star_n_list:
        selected_moments = random.sample(subImage_moments, nstars)
        out_costs.append(get_cost(selected_moments, tt_moment_lists))
    return out_costs

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
    return -b/(2*a) 
    
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
    costs = np.asarray(get_cost_multiple(subImage_moments, tt_moment_lists, range(1,21)))
    focus = [find_focus_position(np.asarray(range(-10,6)), costs[i], plot=plot) for i in range(len(costs))]
    return focus

def focus(catalogs, filenames, tt_galsim_images, tt_star_file, out_name, match_dist = 200., stamp_size = 6., nstars=60, plot=False, generate_new_star_files=True, histogram = True):
    if generate_new_star_files:
        n = 0
        for catalog in catalogs:
            select_good_stars(catalog, catalog+".stars", nstars=nstars)
            print "generating star file", n
            n += 1
    foci = []
    err = []
    for i in range(len(filenames)):
        focus = getMoments(catalogs[i]+".stars", filenames[i], tt_galsim_images, tt_star_file, match_dist=match_dist, stamp_size=stamp_size, plot=plot)
        print "Focus is", focus
        out = open(out_name, "a")
        out.write(filenames[i] + " ")
        for item in focus:
             out.write(str(item) + " ")
        out.write("\n")
        out.close()

for i in range(20):
    focus(["47Tuc_f814w_filter.cat"],["j6ll01ykq_drz.fits"],tt_galsim_images,"/Users/bemi/JPL/F814W_TT/TinyTim_f-1.stars.dat", "focus_uncertainty_814.txt", match_dist=300., stamp_size=6.0)

def uncertainty_boxplot(focus_file, n):
    f = open(focus_file)
    focus_for_n = []
    for i in range(20):
        focus_for_n.append([])
    for line in f.readlines():
        split = line.split()
        for j in range(1,20):
            focus_for_n[j].append(np.float32(split[j]))
    plt.xlabel("Number of stars used for focus calibration")
    plt.ylabel("Measured focus positions")
    plt.title("Distribution of measured focus positions for random samples of n stars for HST ACS image of 47Tuc")
    plt.boxplot(focus_for_n)
    plt.show()
    sample_sd = [np.std(np.asarray(x)) for x in focus_for_n[1:]]
    mean_sd = [x/np.sqrt(n) for x in sample_sd]
    plt.xlabel("Number of stars used for focus calibration")
    plt.ylabel("Standard deviation of sample means of distribution of measured focus positions")
    plt.title("Standard deviation of measured focus positions for random samples of n stars for HST ACS image of 47Tuc")
    plt.scatter(np.asarray(range(2,21)), np.asarray(mean_sd))
    plt.show()
    
#uncertainty_boxplot("focus_uncertainty_606.txt", 20)
#uncertainty_boxplot("focus_uncertainty_814.txt", 12)
            