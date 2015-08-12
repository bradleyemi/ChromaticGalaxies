'''Basic script to use SExtractor to get a mask for each postage stamp.'''

import galsim
import numpy as np
import pyfits
import asciidata
import matplotlib.pyplot as plt
import time
import subprocess
import os

faint_config_dict = { 'DETECT_MINAREA' : 18 ,
    'DETECT_THRESH' : 1.0 ,
    'DEBLEND_NTHRESH' : 64 ,
    'DEBLEND_MINCONT' : 0.065 ,
    'CLEAN_PARAM' : 1.0 ,
    'BACK_SIZE' : 100 ,
    'BACK_FILTERSIZE' : 3 ,
    'BACKPHOTO_TYPE' : "LOCAL" ,
    'BACKPHOTO_THICK' : 200,
    'PIXEL_SCALE' : 0.03}

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

#Get an initial segmentation map for the postage stamps
def run_sextractor(file,weight_file,use_dict=faint_config_dict,output_params=output_params,clean=True): 
    #Create params_file and write out to a file
    param_ascii = asciidata.create(1,len(output_params))
    row_counter = 0
    for param in output_params:
        param_ascii[0][row_counter] = param
        row_counter += 1
    param_fname = file + ".param"
    param_ascii.writeto(param_fname)
    
    #Create config newfiles[i] and write out to a file
    config_ascii = asciidata.create(2,6+len(use_dict))
        
    #File-Specific Configurations
    config_ascii[0][0] = 'CATALOG_NAME'
    config_ascii[1][0] = file + ".cat"
    config_ascii[0][1] = 'PARAMETERS_NAME'
    config_ascii[1][1] = file + ".param"
    config_ascii[0][2] = 'WEIGHT_TYPE'
    config_ascii[1][2] = 'MAP_WEIGHT'
    config_ascii[0][3] = 'WEIGHT_IMAGE'
    config_ascii[1][3] = weight_file
    config_ascii[0][4] = 'CHECKIMAGE_TYPE'
    config_ascii[1][4] = 'SEGMENTATION'
    config_ascii[0][5] = 'CHECKIMAGE_NAME'
    config_ascii[1][5] = file[:len(file)-5] + ".mask.fits"
    row_counter = 6
    for key, value in use_dict.iteritems():
        config_ascii[0][row_counter] = key
        config_ascii[1][row_counter] = value
        row_counter += 1
    config_fname = file + ".config"
    config_ascii.writeto(config_fname)
                
    #Run sextractor and get the catalog
    subprocess.call(["sex", file , "-c", config_fname])
    if clean:
        #Optional Clean
        subprocess.call(["rm", config_fname])
        subprocess.call(["rm", param_fname])

#Make the "real" segmentation map by identifying the object detection in the center of the stamp.

def make_seg_map(init_map):
    #identify centroid
    hdulist = pyfits.open(init_map)
    x_dim = int(hdulist[0].header['NAXIS1'])
    y_dim = int(hdulist[0].header['NAXIS2'])
    x_c = x_dim/2
    y_c = y_dim/2
    data = hdulist[0].data
    mainObjNumber = data[y_c, x_c]
    #Make an empty numpy array of zeros in those dimensions
    segmentation_map_array = np.zeros((x_dim, y_dim))
    for x in range(x_dim):
        for y in range(y_dim):
            if data[y,x] == mainObjNumber:
                segmentation_map_array[y,x] = 1
    #Write out to a fits file
	hdu_out = pyfits.PrimaryHDU(segmentation_map_array)
    hdu_out.writeto(init_map,clobber=True)

#run_sextractor("/Users/bemi/JPL/606_01/images/0.0_214.233_52.44053.processed.fits", "/Users/bemi/JPL/606_01/ivar/0.0_214.233_52.44053.wht.fits")
#make_seg_map("/Users/bemi/JPL/606_01/images/0.0_214.233_52.44053.processed.mask.fits")

### Script to run ###

f = open("data_directories.txt")
directories = []
for line in f.readlines():
    directories.append(line.strip())
f.close()    

root = "/Users/bemi/JPL/"
for dir in ["606_01"]:
    os.chdir(root + dir)
    subprocess.call(["mkdir", "mask"])
    subprocess.call(["mkdir", "cats"])
    image_dir = root + dir + "/images/"
    images = os.listdir(image_dir)
    images = [str(image) for image in images]
    os.chdir(root)
    for image in images:
        if image[len(image)-15:len(image)] != ".processed.fits":
            continue
        if image == ".DS_Store":
            continue
        im = root + dir + "/images/" + image
        weight = root + dir + "/ivar/" + image[:len(image)-15] + ".wht.fits"
        print im, weight
        run_sextractor(im, weight)
        map = image[:len(image)-15] + ".processed.mask.fits"
        print root + dir + "/images/" + map
        make_seg_map(root + dir + "/images/" + map)
        subprocess.call(["mv", root + dir + "/images/" + map, "/Users/bemi/JPL/" + dir + "/mask"])
        subprocess.call(["mv", root + dir + "/images/" + image + ".cat", "/Users/bemi/JPL/" + dir + "/cats"])