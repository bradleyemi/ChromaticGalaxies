#Finds the focus positions of each COSMOS image

import asciidata
import subprocess
import pyfits
import numpy as np
import matplotlib.pyplot as plt
import time
import galsim
from scipy.optimize import curve_fit
import random
import pandas


'''
Steps: 
1) Import all the TinyTim fields
2) Feed an image into GalSim
3) Pick the best S/N stars in the image
4) Get the corresponding TinyTim stars
5) Run pattern match, find focus position
6) Write focus position into catalog
'''

test_image = "EGS_10134_01_acs_wfc_f606w_30mas_unrot_drz.fits"
test_catalog = "EGS_10134_01_acs_wfc_f606w_30mas_unrot_drz.fits.cat"

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

'''     
i = 0
tt_galsim_images = []
for image in tt_606_list:
    print "importing image", i
    f = pyfits.open(image)
    image_data = f[0].data 
    img = galsim.Image(image_data)
    tt_galsim_images.append(img)
    f.close()
    i += 1  
'''
output_params = ["X_IMAGE",
"Y_IMAGE"]

def run_sextractor(file,output_params,out_name,clean=True): 
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
    run_sextractor(tt_606[key], output_params, "TinyTim_f" + str(key) + ".stars.dat", clean=True)
'''
    
def select_good_stars(catalog, out_name, nstars=10):
    cat = asciidata.open(catalog)
    numbers = []
    snrs = []
    for i in range(cat.nrows):
        if cat['IS_STAR'][i] == 1 and cat['MAG_AUTO'][i]+25 < 25.0:
            numbers.append(cat['NUMBER'][i])
            snrs.append(cat['SNR'][i])
    snr_arr = np.asarray(snrs)
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
def match_to_tt(image_star_table, tt_star_data_file, dist=300.):
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
    
def getMoments(image_star_table, image, tt_galsim_images, tt_star_file, stamp_size = 6., plot=True):  
    #One subimage for each star
    subImages = get_subImages(image_star_table, image)
    #One centroid for each star
    tt_centroids = match_to_tt(image_star_table, tt_star_file)
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
        print "accept is", accept
        if accept == "y":
            try:
                subImage_moments.append(galsim.hsm.FindAdaptiveMom(im))
                keep.append(i)
            except:
                print "Adaptive moment did not converge, skipping..."
        i += 1
    j = 10
    print "Keeping stars", keep
    for centroid in tt_centroids:
        if j in keep:
            tt_subImages.append(get_tt_subImages(centroid, tt_galsim_images))
        j += 1
    #One moment for each tt_subimage
    tt_moment_lists = []
    for tt_focus_position in tt_subImages:
        tt_moments = []
        for tt_subImage in tt_focus_position:
            tt_moments.append(galsim.hsm.FindAdaptiveMom(tt_subImage))
        tt_moment_lists.append(tt_moments)
    print len(subImage_moments), len(tt_moment_lists), len(tt_moment_lists[0])
    cost = np.asarray(get_cost(subImage_moments, tt_moment_lists))
    focus, focus_err = find_focus_position(np.asarray(range(-10,6)), cost, plot=plot)
    return focus, focus_err
    
f = open("f606w_catalogs.txt")   
f606w_catalogs = []
for line in f.readlines():
    f606w_catalogs.append(line.strip())
f.close()

g = open("f606w_filenames.txt")
f606w_files = []
for line in g.readlines():
    f606w_files.append(line.strip())
g.close()

'''
n = 0
for catalog in f606w_catalogs:
    select_good_stars(catalog, catalog + ".stars", nstars=20)
    print "generating star file", n
    n += 1
'''
'''    
foci = []
err = []
n = 0
for file in f606w_files:
    print "working on file", n
    n += 1
    try:
        focus, focus_err = getMoments(file + ".cat.stars", file, tt_galsim_images, "TinyTim_f-1.stars.dat", plot=True)
        print "Focus is", focus, "plus/minus", focus_err 
        foci.append(focus)
        err.append(focus_err)
    except:
        raise
        foci.append(None)
        err.append(None)
        

plt.hist(foci, bins=20)
plt.xlabel("focus position (um)")
plt.ylabel("frequency")
plt.title("Focus positions of f606w COSMOS images")
plt.show()       
    
out = open("606_focus_positions.txt", "w")
for i in range(len(foci)):
    out.write(f606w_files[i] + " ")
    out.write(str(foci[i]) + " ")
    out.write(str(err[i]) + "\n")
'''

f = open("606_focus_positions.txt")
times = []
focuses = []
for line in f.readlines():
    split = line.split()
    image = split[0]
    focus = np.float32(split[1])
    hdulist = pyfits.open(image)
    time = np.float32(hdulist[0].header['EXPSTART'])
    hdulist.close()
    if time < 53225:
       times.append(time)
       focuses.append(focus)
f.close()

#print times
#print focuses
plt.scatter(times,focuses)
plt.title("Interpolated focus position of camera vs. time")
plt.xlabel("Modified Julian Day Number")
plt.ylabel("Focus Position (um)")
plt.show()
    
    
