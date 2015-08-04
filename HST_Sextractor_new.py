'''
Script Name: HST_sextractor.py
Author: Bradley Emi
Date: July 2015

v. 1.0 (7/20/2015) script does edge removal, star-diffraction spikes through run_pipeline. Manual masking and in-image overlap
has to be done separately.

v. 2.0 (8/4/2015) changed to object-oriented design. Manual masking, in-image overlap are included. A clean catalog is generated every time
a GalaxyCatalog instance is initiated. (No statistics)

########### Description ##########

This script takes a file, and uses the hot-cold method described in Leauthaud (2007) to select
galaxies and stars from the COSMOS sample using the SExtractor software. The catalogs are set to output the
parameters in output_params. It then cleans the catalogs for:

-Edge removal
-Star diffraction spikes
-Manual removal of image defects
-In-image overlap

The parameters for star-galaxy classification can be adjusted, and checked using make_classifier_plot.
The automatic star diffraction masking has been tested on AEGIS f606w and f814w filters, and assumes
that the size of each diffraction spike scales with flux. These parameters can also be adjusted, but the 
automatic parameters are included. 

Finally, statistics of how many objects are extracted and cleaned at each step are output into a text file.

########## Input ##########

image (.fits file)
background (.fits file)

########## Dependencies ##########

AstroAsciiData
PyFits/AstroPy
NumPy
MatPlotLib.PyPlot

########## Usage ##########

cat = GalaxyCatalog(file, background, filter, manual_mask_file, out_name)
cat.make_catalog()

'''

import asciidata
import subprocess
import pyfits
import numpy as np
import matplotlib.pyplot as plt
import time

'''
Steps:
1) Run sextractor for the bright objects
2) Run sextractor for the faint objects
3) Make a segmentation map out of the bright star catalog
4) Filter the faint catalog with the segmentation map
5) Merge the bright catalog and the filtered faint catalog
6) Fix the indexes of the merged catalog
7) Find the four parametric values (x,y,slope,intercept) on the mu/mag plot for S-G classification
8) Feed those values to alter_catalog_for_classification
9) Add the signal to noise ratio
10) Edge removal
11) Diffraction spike cleanup
12) Overlap
13) Manual cleanup
'''

### Some useful helper functions ###

#Renumbers the objects in a catalog
def renumber(catalog):
    cat = asciidata.open(catalog)
    for i in range(cat.nrows):
        cat['NUMBER'][i] = i
    cat.writeto(catalog)
    
#Takes two x-y pairs (tuples) and returns as a tuple the slope and intercept of the line connecting them
def points_to_line(point1, point2):
    x1 = point1[0]
    y1 = point1[1]
    x2 = point2[0]
    y2 = point2[1]
    slope = (y2-y1)/(x2-x1)
    intercept = -1*slope*x1+y1
    return (slope, intercept)
    
def is_below_boundary(input_x, input_y, flat_x_division, flat_y_division, slope, intercept):
    output = True
    if input_x < flat_x_division and input_y > flat_y_division:
        output = False
    if input_x > flat_x_division and input_y > slope*input_x + intercept:
        output = False
    return output
    
#rotate the point x,y through an angle theta in degrees about x0,y0
def rotate(x,y,x0,y0,theta):
    theta = np.radians(theta)
    x_rotated = []
    y_rotated = []
    for i in range(len(x)):
        x_rotated.append(np.cos(theta)*x[i]-np.sin(theta)*y[i]+(1-np.cos(theta))*x0+np.sin(theta)*y0)
        y_rotated.append(np.sin(theta)*x[i]+np.cos(theta)*y[i]-np.sin(theta)*x0+(1-np.cos(theta))*y0)
    return [x_rotated,y_rotated]

#Use the jordan curve theorem to determine whether px,py is inside a polygon with ordered vertex lists x,y
def inpoly(px,py,x,y):
    crossings = 0
    x = np.array(x)
    y = np.array(y)
    x = np.append(x,x[0])
    y = np.append(y,y[0])
    #Change to px,py coordinate system
    x = x - px
    y = y - py
    if sum(y) - sum(abs(y)) == 0 or sum(y) + sum(abs(y)) == 0:
        return 0
    for i in range(len(x)-1):
        Ax = x[i]
        Ay = y[i]
        Bx = x[i+1]
        By = y[i+1]
        if Ay*By < 0:
            if Ax>0 and Bx>0:
                crossings += 1
            else:
                c = Ax-(Ay*(Bx-Ax))/(By-Ay)
                if c > 0:
                    crossings += 1
    return crossings % 2
    
def delete_null(catalog):
    f = open(catalog)
    lines = f.readlines()
    f.close()
    g = open(catalog, "w")
    for line in lines:
        if (line.split())[0] != "Null":
            g.write(line)        

#tolerance in degrees
def delete_overlap(catalog, tolerance = 1./18000, clean=True):
    s = time.time()
    orig_name = catalog
    print "Deleting overlaps on file: ", orig_name
    catalog = asciidata.open(catalog)
    new_table = asciidata.create(catalog.ncols,catalog.nrows)
    delete_numbers = []
    for i in range(catalog.nrows-1):
        if (i+1) % 10 == 0:
            print "Working on object", i+1, "out of", catalog.nrows
        number = catalog['NUMBER'][i]
        ra, dec = catalog['ALPHA_SKY'][i], catalog['DELTA_SKY'][i]
        for j in range(i+1, catalog.nrows):
            number0 = catalog['NUMBER'][j]
            ra0, dec0 = catalog['ALPHA_SKY'][j], catalog['DELTA_SKY'][j]
            if abs(ra0-ra) < tolerance and abs(dec0-dec) < tolerance:
                delete_numbers.append(number0)
                delete_numbers.append(number)
    print "Delete numbers", delete_numbers
    new_table = asciidata.create(catalog.ncols,catalog.nrows)
    for i in range(catalog.nrows):
        if catalog['NUMBER'][i] not in delete_numbers:
             for k in range(catalog.ncols):
                 new_table[k][i] = catalog[k][i]
    #Get rid of empty rows
    row_number = 0
    while True:
        try:
            if new_table[0][row_number] is None:
                new_table.delete(row_number)
                if "606" in orig_name:
                    n_overlap606 += 1
                if "814" in orig_name:
                    n_overlap814 += 1
            else:
                row_number += 1
        except:
            break
    #Write out to another catalog
    new_table.writeto("overlap.cat")
    new_catalog = open("overlap.cat")
    old_catalog = open(orig_name)
    final_catalog = open("final_overlap.cat", "w")
    for line in old_catalog.readlines():
        if line[0] == "#":
            final_catalog.write(line)
        else:
            break
    for line in new_catalog:
        final_catalog.write(line)
    old_catalog.close()
    final_catalog.close()
    final = open("final_overlap.cat", "r")
    rewrite = open(orig_name, "w")
    for line in final.readlines():
        rewrite.write(line)
    #Optional clean
    if clean:
        subprocess.call(["rm", "overlap.cat"])
        subprocess.call(["rm", "final_overlap.cat"])
    e = time.time()
    print "Time:", e-s

def manual_mask(catalog, x_vertices, y_vertices, clean=True):
    orig_name = catalog
    catalog = asciidata.open(catalog)
    new_table = asciidata.create(catalog.ncols,catalog.nrows)
    delete_numbers = []
    for i in range(catalog.nrows):
        if (i+1) % 1000 == 0:
            print "Working on object", i+1, "out of", catalog.nrows
        number = catalog['NUMBER'][i]
        is_star = catalog['IS_STAR'][i]
        if is_star == 1:
            continue
        x_min,y_min = catalog['XMIN_IMAGE'][i], catalog['YMIN_IMAGE'][i]
        x_max,y_max = catalog['XMAX_IMAGE'][i], catalog['YMAX_IMAGE'][i]
        bottom_pixels = [(x,y_min) for x in range(x_min,x_max)]
        left_pixels = [(x_min,y) for y in range(y_min,y_max)]
        top_pixels = [(x, y_max) for x in range(x_min,x_max)]
        right_pixels = [(x_max,y) for y in range(y_min,y_max)]
        pixels = bottom_pixels + left_pixels + top_pixels + right_pixels
    	bools = [inpoly(pixel[0],pixel[1],x_vertices,y_vertices) for pixel in pixels]
        if max(bools) == 1:
            delete_numbers.append(number)
    print "Delete numbers", delete_numbers
    new_table = asciidata.create(catalog.ncols,catalog.nrows)
    for i in range(catalog.nrows):
        if catalog['NUMBER'][i] not in delete_numbers:
             for k in range(catalog.ncols):
                 new_table[k][i] = catalog[k][i]
    #Get rid of empty rows
    row_number = 0
    while True:
        try:
            if new_table[0][row_number] is None:
                new_table.delete(row_number)
            else:
                row_number += 1
        except:
            break
    #Write out to another catalog
    new_table.writeto("manual_filter.cat")
    new_catalog = open("manual_filter.cat")
    old_catalog = open(orig_name)
    final_catalog = open("final_catalog.cat", "w")
    for line in old_catalog.readlines():
        if line[0] == "#":
            final_catalog.write(line)
        else:
            break
    for line in new_catalog:
        final_catalog.write(line)
    old_catalog.close()
    final_catalog.close()
    final = open("final_catalog.cat", "r")
    rewrite = open(orig_name, "w")
    for line in final.readlines():
        rewrite.write(line)
    #Optional clean
    if clean:
        subprocess.call(["rm", "final_catalog.cat"])
        subprocess.call(["rm", "manual_filter.cat"])

### Class definition ###
class GalaxyCatalog:
    
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
    
    star_galaxy_weights = (19.0, -9.8, 0.9, -26.9)
    
    f606w_spike_params = (0.0350087,64.0863,40.0,2.614)
    f814w_spike_params = (0.0367020,77.7674,40.0,2.180)
    
    def __init__(self, file, weight_file, filter, out_name, manual_mask_file=None):  
        #Initial setup of attributes
        self.file = file
        self.weight_file = weight_file
        self.filter = filter
        self.spike_params = None
        if self.filter == 606:
            self.spike_params = self.f606w_spike_params
        if self.filter == 814:
            self.spike_params = self.f814w_spike_params
        self.out_name = out_name
        self.catalog_vertex_file = manual_mask_file
        
    def generate_catalog(self):
        #Runs sextractor for the bright catalog
        print self.output_params
        self.__run_sextractor(self.bright_config_dict, self.out_name + "_bright", self.output_params)
        self.bright_catalog = self.out_name + "_bright.cat"
        #Stores "bright.cat" as the bright catalog
        
        #Makes the segmentation map
        self.__make_segmentation_map(self.out_name)
        self.seg_map = self.out_name + "_seg_map.fits"
        
        #Runs sextractor for the faint catalog
        self.__run_sextractor(self.faint_config_dict, self.out_name + "_faint", self.output_params)
        self.faint_catalog = self.out_name + "_faint.cat"
        
        #Filters the faint catalog
        self.__filter_cat_with_segmentation_map(self.out_name + "_filteredfaint.cat")
        self.filtered_faint_catalog = self.out_name + "_filteredfaint.cat"
        
        #Merges the faint and bright catalogs
        self.__merge(self.out_name + "_merge.cat")
        self.merged_catalog = self.out_name + "_merge.cat"
        renumber(self.merged_catalog)
        
        #Add is_star
        self.__alter_catalog_for_classification(self.out_name + "_class.cat", self.star_galaxy_weights[0], self.star_galaxy_weights[1], self.star_galaxy_weights[2], self.star_galaxy_weights[3])
        self.class_catalog = self.out_name + "_class.cat"
        
        #Add S/N ratio
        self.__make_SNR(self.out_name + "_snr.cat")
        self.snr_catalog = self.out_name + "_snr.cat"   
        
        #Clean out edge objects
        self.__edge_overlap_clean(self.out_name + "_edge.cat")
        self.edge_catalog = self.out_name + "_edge.cat"
        
        #Clean star diffraction spikes/clean for overlap
        self.diffraction_mask_cleanup(self.out_name + "_diff.cat", self.spike_params)
        delete_overlap(self.out_name + "_diff.cat")
        delete_null(self.out_name + "_diff.cat")
        renumber(self.out_name + "_diff.cat")
        self.diff_catalog = self.out_name + ".cat"
        self.get_catalog = self.diff_catalog
        
        #Manual mask
        if self.catalog_vertex_file != None:
            self.manual_mask_catalogs()
        
        subprocess.call(["rm", self.out_name + "_*"])
                
    def run_sextractor(self, use_dict, out_name, output_params, clean=True):
        param_ascii = asciidata.create(1,len(output_params))
        row_counter = 0
        for param in output_params:
            param_ascii[0][row_counter] = param
            row_counter += 1
        param_fname = self.out_name + ".param"
        param_ascii.writeto(param_fname)
        #Create config newfiles[i] and write out to a file
        config_ascii = asciidata.create(2,4+len(use_dict))
        #File-Specific Configurations
        config_ascii[0][0] = 'CATALOG_NAME'
        config_ascii[1][0] = out_name + ".cat"
        config_ascii[0][1] = 'PARAMETERS_NAME'
        config_ascii[1][1] = self.out_name + ".param"
        config_ascii[0][2] = 'WEIGHT_TYPE'
        config_ascii[1][2] = 'MAP_WEIGHT'
        config_ascii[0][3] = 'WEIGHT_IMAGE'
        config_ascii[1][3] = self.weight_file
        row_counter = 4
        for key, value in use_dict.iteritems():
            config_ascii[0][row_counter] = key
            config_ascii[1][row_counter] = value
            row_counter += 1
        config_fname = out_name + ".config"
        config_ascii.writeto(config_fname)
                
        #Run sextractor and get the catalog
        subprocess.call(["sex", self.file , "-c", config_fname])
    
        #Optional Clean
        subprocess.call(["rm", config_fname])
        subprocess.call(["rm", param_fname])
    
    __run_sextractor = run_sextractor
    
    def make_segmentation_map(self, out_name, enlarge=20):
        hdulist = pyfits.open(self.file)
        x_dim = int(hdulist[0].header['NAXIS1'])
        y_dim = int(hdulist[0].header['NAXIS2'])
        hdulist.close()
        positions = []
        catalog = asciidata.open(self.bright_catalog)
        for i in range(catalog.nrows):
            x_min = catalog['XMIN_IMAGE'][i]
            y_min = catalog['YMIN_IMAGE'][i]
            x_max = catalog['XMAX_IMAGE'][i]
            y_max = catalog['YMAX_IMAGE'][i]
            pos_tuple = (x_min, x_max, y_min, y_max)
            positions.append(pos_tuple)
        #Make an empty numpy array of zeros in those dimensions
        segmentation_map_array = np.zeros((x_dim, y_dim))
        #Iterate through position tuples and switch flagged areas to 1's in the array
        for j in range(len(positions)):
            x_min = (positions[j])[0]
            x_max = (positions[j])[1]
            y_min = (positions[j])[2]
            y_max = (positions[j])[3]
            for x in range(x_min-enlarge,x_max+enlarge):
                for y in range(y_min-enlarge,y_max+enlarge):
                    try:
                        segmentation_map_array[y,x] = 1
                    except:
                        continue
        #Write out to a fits file
	    hdu_out = pyfits.PrimaryHDU(segmentation_map_array)
        hdu_out.writeto(out_name + "_seg_map.fits",clobber=True)
        
    __make_segmentation_map = make_segmentation_map
    
    def filter_cat_with_segmentation_map(self, out_name):
        catalog = asciidata.open(self.faint_catalog)
        new_table = asciidata.create(catalog.ncols,catalog.nrows) 
        segmentation_file = pyfits.open(self.seg_map)
        data = segmentation_file[0].data   
        for j in range(0,catalog.nrows):
            x_center = catalog['X_IMAGE'][j]
            y_center = catalog['Y_IMAGE'][j]
            if data[y_center,x_center] == 0:
    	        for k in range(0,catalog.ncols):
    	            new_table[k][j] = catalog[k][j]
    	#Get rid of empty rows
        row_number = 0
        while True:
            try:
                if new_table[0][row_number] is None:
                    new_table.delete(row_number)
                else:
                    row_number += 1
            except:
                break
        #Write out to another catalog
        new_table.writeto(out_name)  
        
    __filter_cat_with_segmentation_map = filter_cat_with_segmentation_map
    
    #Merges a bright catalog and a filtered faint catalog (with or without a header) into a single catalog with header from the bright catalog
    def merge(self, out_name):
        #Copy the header to a new file
        f = open(out_name, "w")
        faint_cat = open(self.filtered_faint_catalog)
        bright_cat = open(self.bright_catalog)
        for line in bright_cat.readlines():
            f.write(line)
        for line in faint_cat.readlines():
    	    if line[0] != "#":
                f.write(line)          
        faint_cat.close()
        bright_cat.close()
        
    __merge = merge
    
    def alter_catalog_for_classification(self, out_name, flat_x_division, flat_y_division, slope, intercept):
        catalog = asciidata.open(self.merged_catalog)
        for i in range(catalog.nrows):
            if is_below_boundary(catalog['MAG_AUTO'][i]+25, catalog['MU_MAX'][i], flat_x_division, flat_y_division, slope, intercept) and catalog['MAG_AUTO'][i]+25.0 < 25.0:
                catalog['IS_STAR'][i] = 1
            else:
                catalog['IS_STAR'][i] = 0
        catalog['IS_STAR'].set_colcomment("Revised Star-Galaxy Classifier")
        catalog.writeto(out_name)
    
    __alter_catalog_for_classification = alter_catalog_for_classification
    
    def make_SNR(self, out_name):
        catalog = asciidata.open(self.class_catalog)
        for i in range(catalog.nrows):
             catalog['SNR'][i] = catalog['FLUX_AUTO'][i]/catalog['FLUXERR_AUTO'][i]
        catalog['SNR'].set_colcomment("Signal to Noise Ratio")
        catalog.writeto(out_name)
        
    __make_SNR = make_SNR
    
    def edge_overlap_clean(self, out_name):
        A = (390.,321.)
        B = (498.,6725.)
        C = (6898.,7287.)
        D = (7002.,806.)
        (left_m, left_b) = points_to_line(A,B)
        (top_m, top_b) = points_to_line(B,C)
        (right_m, right_b) = points_to_line(C,D)
        (bottom_m, bottom_b) = points_to_line(D,A)
        #Open the original catalog
        table_catalog = asciidata.open(self.snr_catalog)
        #Make a new ascii table (the new catalog)
        new_table = asciidata.create(table_catalog.ncols,table_catalog.nrows)
        for i in range(table_catalog.nrows):
            x_min = table_catalog['XMIN_IMAGE'][i]
            x_max = table_catalog['XMAX_IMAGE'][i]
            y_min = table_catalog['YMIN_IMAGE'][i]
            y_max = table_catalog['YMAX_IMAGE'][i]
            if (x_min < left_m*x_min+left_b and x_max < right_m*x_max+right_b and \
            y_min > bottom_m*y_min+bottom_b and y_max < top_m*y_max+top_b):
                 for k in range(table_catalog.ncols):
                     new_table[k][i] = table_catalog[k][i]
        #Get rid of empty rows
        row_number = 0
        while True:
            try:
                if new_table[0][row_number] is None:
                    new_table.delete(row_number)
                else:
                    row_number += 1
            except:
                break
        #Write out to another catalog
        new_table.writeto("table_edge_clean.cat")
        new_catalog = open("table_edge_clean.cat")
        old_catalog = open(self.snr_catalog)
        final_catalog = open(out_name, "w")
        for line in old_catalog.readlines():
            if line[0] == "#":
                final_catalog.write(line)
            else:
                break
        for line in new_catalog:
            final_catalog.write(line)
    
    __edge_overlap_clean = edge_overlap_clean
    
	#Clean the merged catalog
    def diffraction_mask_cleanup(self, out_name, diff_spike_params, mag_cutoff = 19.0):
        #Get all stars for which we are masking
        #the length of the spike is given by l = m*Flux+b
        #the width of the spike is given by w
        m = diff_spike_params[0] 
        b = diff_spike_params[1]
        w = diff_spike_params[2]*0.5
        theta = diff_spike_params[3]
        star_numbers = []
        star_centroids = []
        star_radii = []
        star_spike_lengths = []
        orig_name = self.edge_catalog
        catalog = asciidata.open(self.edge_catalog)
        for i in range(catalog.nrows):
           if catalog['MAG_AUTO'][i] + 25 < mag_cutoff and catalog['IS_STAR'][i] == 1:
              star_numbers.append(catalog['NUMBER'][i])
              star_centroids.append((catalog['X_IMAGE'][i],catalog['Y_IMAGE'][i]))
              radius = np.mean([catalog['A_IMAGE'][i],catalog['B_IMAGE'][i]])
              star_radii.append(radius)
              flux = catalog['FLUX_AUTO'][i]
              length = m*flux+b
              star_spike_lengths.append(length)
        #Get the vertices of each star's diffraction mask
        x_vertex_sets = []
        y_vertex_sets = []
        for i in range(len(star_numbers)):
            centroid = (star_centroids[i])
            x0 = centroid[0]
            y0 = centroid[1]
            r = star_radii[i]
            l = star_spike_lengths[i]
            x_vertices = [x0-w,x0-w,x0+w,x0+w,x0+r,x0+l,x0+l,x0+r,x0+w,x0+w,x0-w,x0-w,x0-r,x0-l,x0-l,x0-r]
            y_vertices = [y0+r,y0+l,y0+l,y0+r,y0+w,y0+w,y0-w,y0-w,y0-r,y0-l,y0-l,y0-r,y0-w,y0-w,y0+w,y0+w]
            rotated = rotate(x_vertices,y_vertices,x0,y0,theta)
            x_rotated = rotated[0]
            y_rotated = rotated[1]
            x_vertex_sets.append(x_rotated)
            y_vertex_sets.append(y_rotated)
        delete_numbers = []
        print "Applying masks for", len(x_vertex_sets), "stars."
        for i in range(catalog.nrows):
            if (i+1) % 1000 == 0:
                print "Working on object", i+1, "out of", catalog.nrows
            number = catalog['NUMBER'][i]
            is_star = catalog['IS_STAR'][i]
            if is_star == 1:
                continue
            x_min,y_min = int(catalog['XMIN_IMAGE'][i]), int(catalog['YMIN_IMAGE'][i])
            x_max,y_max = int(catalog['XMAX_IMAGE'][i]), int(catalog['YMAX_IMAGE'][i])
            bottom_pixels = [(x,y_min) for x in range(x_min,x_max)]
            left_pixels = [(x_min,y) for y in range(y_min,y_max)]
            top_pixels = [(x, y_max) for x in range(x_min,x_max)]
            right_pixels = [(x_max,y) for y in range(y_min,y_max)]
            pixels = bottom_pixels + left_pixels + top_pixels + right_pixels
            bools = [inpoly(pixel[0],pixel[1],x_vertex_sets[j],y_vertex_sets[j]) for pixel in pixels for j in range(len(x_vertex_sets))]
            if max(bools) == 1:
                delete_numbers.append(number)
        print "Delete numbers", delete_numbers
        #Delete entries for which any pixel is within the mask
        #Make a new ascii table (the new catalog)
        new_table = asciidata.create(catalog.ncols,catalog.nrows)
        for i in range(catalog.nrows):
            if catalog['NUMBER'][i] not in delete_numbers:
                for k in range(catalog.ncols):
                    new_table[k][i] = catalog[k][i]
        #Get rid of empty rows
        row_number = 0
        while True:
            try:
                if new_table[0][row_number] is None:
                     new_table.delete(row_number)
                else:
                     row_number += 1
            except:
                break
        #Write out to another catalog
        new_table.writeto("diffraction_filter.cat")
        new_catalog = open("diffraction_filter.cat")
        old_catalog = open(orig_name)
        final_catalog = open(out_name, "w")
        for line in old_catalog.readlines():
            if line[0] == "#":
                final_catalog.write(line)
            else:
                break
        for line in new_catalog:
            final_catalog.write(line) 
    
	# Runs manual mask on a file where lines are:
	# # 'filename' (must be preceded by # and a space)
	# list of x-vertices for mask 1
	# list of y-vertices for mask 1
	# list of x-vertices for mask 2, etc...
	def manual_mask_catalogs(self, catalog_vertex_file=self.catalog_vertex_file):
		f = open(catalog_vertex_file)
		current_catalog = ""
		x_vertices = []
		y_vertices = []
		for line in f.readlines():
			split = line.split()
			if split[0] == '#':
				current_catalog = split[1]
				x_vertices = []
				y_vertices = []
			else:
				if len(x_vertices) == 0:
					for word in split:
						x_vertices.append(np.float32(word))
				else:
					for word in split:
						y_vertices.append(np.float32(word))
					print "Masking vertices on", current_catalog, "with vertices:"
					print x_vertices
					print y_vertices
					manual_mask(current_catalog, x_vertices, y_vertices)
					x_vertices = []
					y_vertices = []    
					 

'''
#Makes the mu_max vs. mag_auto plot
def make_classifier_plot(merged_catalog, flat_x_division=19.0, flat_y_division=-9.8, dividing_x1=19.0, dividing_y1 = -9.8, dividing_x2 = 21.0, dividing_y2 = -8.0):
    catalog = asciidata.open(merged_catalog)
    mu = []
    mag = []
    class_star = []
    for j in range(catalog.nrows):
	    mu.append( catalog['MU_MAX'][j] )
	    mag.append( catalog['MAG_AUTO'][j] + 25 )
	    class_star.append( catalog['CLASS_STAR'][j] )
    star_mu = []
    star_mag = []
    gal_mu = []
    gal_mag = []
    dividing_x1 = 19.0
    dividing_y1 = -9.8
    dividing_x2 = 21.0
    dividing_y2 = -8.0
    flat_x_division = 19.0
    flat_y_division = -9.8
    dividing_line = points_to_line((dividing_x1,dividing_y1),(dividing_x2,dividing_y2))
    star_detections = 0
    print dividing_line
    for i in range(len(class_star)):
        if mag[i] < flat_x_division:
            if mu[i] < flat_y_division:
                star_mu.append(mu[i])
                star_mag.append(mag[i])
                star_detections += 1
            else:
                gal_mu.append(mu[i])
                gal_mag.append(mag[i])
        if mag[i] > flat_x_division and mag[i] < 25:
            if mu[i] < dividing_line[0]*mag[i]+dividing_line[1]:
                star_mu.append(mu[i])
                star_mag.append(mag[i])
                star_detections += 1
            else:
                gal_mu.append(mu[i])
                gal_mag.append(mag[i])

    plt.scatter(gal_mag, gal_mu,c="BLUE")
    plt.scatter(star_mag, star_mu,c="RED")
    plt.plot([15.0,dividing_x1,dividing_x2,25.0],[dividing_y1,dividing_y1,dividing_y2,dividing_line[0]*25.0+dividing_line[1]])
    plt.xlim((15.0,30.0))
    plt.xlabel("MAG_AUTO")
    plt.ylabel("MU_MAX")
    plt.title("Surface brightness vs. apparent magnitude")
    plt.show()
    print dividing_line
'''

    
