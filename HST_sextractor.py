#Sextractor catalogs for HST Data

import asciidata
import subprocess
import pyfits
import numpy as np
import matplotlib.pyplot as plt
import time

'''
Steps so far:
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
12) Manual cleanup (after pipeline)
13) Overlap
'''


f606w_files = []
f814w_files = []
f606w_backgrounds = []
f814w_backgrounds = []

f606w = open("f606w_filenames.txt")
for line in f606w.readlines():
    f606w_files.append(line.strip())
f814w = open("f814w_filenames.txt")
for line in f814w.readlines():
    f814w_files.append(line.strip())
    
f606wb = open("f606w_backgrounds.txt")
for line in f606wb.readlines():
    f606w_backgrounds.append(line.strip())

f814wb = open("f814w_backgrounds.txt")
for line in f814wb.readlines():
    f814w_backgrounds.append(line.strip())
   
### Statistics

n_c1 = 0
n_bright = 0
n_faint = 0
n_hotcold = 0
n_border = 0
n_diff = 0
n_manual606 = 0
n_manual814 = 0
n_overlap606 = 0
n_overlap814 = 0
#Testing files

test_background = "EGS_10134_0a_acs_wfc_f814w_30mas_unrot_wht.fits"
test_image = "EGS_10134_0a_acs_wfc_f814w_30mas_unrot_drz.fits"
test_background2 = "EGS_10134_0a_acs_wfc_f606w_30mas_unrot_wht.fits"
test_image2 = "EGS_10134_0a_acs_wfc_f606w_30mas_unrot_drz.fits"
 
#Set the output parameters

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

#Runs sextractor on a list of filename strings, config_dict and 

def run_sextractor(file,weight_file,use_dict,output_params,out_name,clean=True): 
    #Make sure the file opens
    try:
        hdulist = pyfits.open(file)
        hdulist.close()
    except IOError:
        print "IO Error for filename: ", file
    try:
        hdulist2 = pyfits.open(weight_file)
        hdulist2.close()
    except IOError:
        print "IO Error for filename: ", weight_file

    #Create params_file and write out to a file
    param_ascii = asciidata.create(1,len(output_params))
    row_counter = 0
    for param in output_params:
        param_ascii[0][row_counter] = param
        row_counter += 1
    param_fname = out_name + ".param"
    param_ascii.writeto(param_fname)
    
    #Create config newfiles[i] and write out to a file
    config_ascii = asciidata.create(2,4+len(use_dict))
        
    #File-Specific Configurations
    config_ascii[0][0] = 'CATALOG_NAME'
    config_ascii[1][0] = out_name + ".cat"
    config_ascii[0][1] = 'PARAMETERS_NAME'
    config_ascii[1][1] = out_name + ".param"
    config_ascii[0][2] = 'WEIGHT_TYPE'
    config_ascii[1][2] = 'MAP_WEIGHT'
    config_ascii[0][3] = 'WEIGHT_IMAGE'
    config_ascii[1][3] = weight_file
    row_counter = 4
    for key, value in use_dict.iteritems():
        config_ascii[0][row_counter] = key
        config_ascii[1][row_counter] = value
        row_counter += 1
    config_fname = out_name + ".config"
    config_ascii.writeto(config_fname)
                
    #Run sextractor and get the catalog
    subprocess.call(["sex", file , "-c", config_fname])
    
    #Optional Clean
    subprocess.call(["rm", config_fname])
    subprocess.call(["rm", param_fname])

#Makes a segmentation map (.fits file) from a single sextractor catalog (string = to name of text file, image_file is a string with the name of the .fits file for the image)
def make_segmentation_map(sextractor_catalog,image_file,out_name,enlarge=20):
    #Get dimensions of the original image
    hdulist = pyfits.open(image_file)
    x_dim = int(hdulist[0].header['NAXIS1'])
    y_dim = int(hdulist[0].header['NAXIS2'])
    hdulist.close()
    
    #Get position tuples of the detected objects and round them to the nearest integer
    positions = []
    catalog = asciidata.open(sextractor_catalog)
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
    hdu_out.writeto(out_name + ".fits",clobber=True)

#Makes a second catalog with the entries whose center is in a flagged region of a segmentation map deleted
def filter_cat_with_segmentation_map(faint_catalog, segmentation_map, out_name, clean=False):
    #Open the original catalog
    catalog = asciidata.open(faint_catalog)
    #Make a new ascii table (the new catalog)
    new_table = asciidata.create(catalog.ncols,catalog.nrows)
    #Iterate through each line of the catalog, and copy it if its center is not flagged
    segmentation_file = pyfits.open(segmentation_map)
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

    #Optional clean
    if clean:
        subprocess.call(["rm", faint_catalog])

#Merges a bright catalog and a filtered faint catalog (with or without a header) into a single catalog with header from the bright catalog
def merge_bright_faint_catalogs(bright_catalog, faint_catalog, out_name):
    #Copy the header to a new file
    f = open(out_name, "w")
    faint_cat = open(faint_catalog)
    bright_cat = open(bright_catalog)
    for line in bright_cat.readlines():
        f.write(line)
    for line in faint_cat.readlines():
    	if line[0] != "#":
            f.write(line)          
    faint_cat.close()
    bright_cat.close()

def replace_indexes_of_merged_catalog(merged_catalog, out_name):
    catalog = asciidata.open(merged_catalog)
    for i in range(catalog.nrows):
        catalog['NUMBER'][i] = i
    catalog.writeto(out_name)

#Takes two x-y pairs (tuples) and returns as a tuple the slope and intercept of the line connecting them
def points_to_line(point1, point2):
    x1 = point1[0]
    y1 = point1[1]
    x2 = point2[0]
    y2 = point2[1]
    slope = (y2-y1)/(x2-x1)
    intercept = -1*slope*x1+y1
    return (slope, intercept)

#Makes the mu_max vs. mag_auto plot
def make_classifier_plot(merged_catalog):
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
    print star_detections



def is_below_boundary(input_x, input_y, flat_x_division, flat_y_division, slope, intercept):
    output = True
    if input_x < flat_x_division and input_y > flat_y_division:
        output = False
    if input_x > flat_x_division and input_y > slope*input_x + intercept:
        output = False
    return output

def alter_catalog_for_classification(merged_catalog, out_name, flat_x_division=19.0, flat_y_division=-9.8, slope=0.9, intercept=-26.9):
    catalog = asciidata.open(merged_catalog)
    for i in range(catalog.nrows):
        if is_below_boundary(catalog['MAG_AUTO'][i]+25, catalog['MU_MAX'][i], flat_x_division, flat_y_division, slope, intercept) and catalog['MAG_AUTO'][i]+25.0 < 25.0:
            catalog['IS_STAR'][i] = 1
        else:
            catalog['IS_STAR'][i] = 0
    catalog['IS_STAR'].set_colcomment("Revised Star-Galaxy Classifier")
    catalog.writeto(out_name)

def make_SNR(merged_catalog, out_name):
    catalog = asciidata.open(merged_catalog)
    for i in range(catalog.nrows):
        catalog['SNR'][i] = catalog['FLUX_AUTO'][i]/catalog['FLUXERR_AUTO'][i]
    catalog['SNR'].set_colcomment("Signal to Noise Ratio")
    catalog.writeto(out_name)

def edge_overlap_clean(catalog, out_name, clean=False):
    A = (390.,321.)
    B = (498.,6725.)
    C = (6898.,7287.)
    D = (7002.,806.)
    (left_m, left_b) = points_to_line(A,B)
    (top_m, top_b) = points_to_line(B,C)
    (right_m, right_b) = points_to_line(C,D)
    (bottom_m, bottom_b) = points_to_line(D,A)
    #Open the original catalog
    table_catalog = asciidata.open(catalog)
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
    old_catalog = open(catalog)
    final_catalog = open(out_name, "w")
    for line in old_catalog.readlines():
        if line[0] == "#":
            final_catalog.write(line)
        else:
            break
    for line in new_catalog:
        final_catalog.write(line)
    #Optional clean
    if clean:
        subprocess.call(["rm", catalog])
        subprocess.call(["rm", "table_edge_clean.cat"])
     
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
        
#Clean the merged catalog
#Spike parameters are: slope, intercept (of linear model mapping flux to spike length), width of spike, and theta (angle spikes are rotated) in degrees
f606w_spike_params = (0.0350087,64.0863,40.0,2.614)
f814w_spike_params = (0.0367020,77.7674,40.0,2.180)
def diffraction_mask_cleanup(catalog, out_name, mag_cutoff = 19.0, diff_spike_params = f606w_spike_params,clean=False):
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
    orig_name = catalog
    catalog = asciidata.open(catalog)
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
    #Optional clean
    if clean:
        subprocess.call(["rm", catalog])
        subprocess.call(["rm", "diffraction_filter.cat"])

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
    global n_manual606
    global n_manual814
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
                if "606" in orig_name:
                    n_manual606 += 1
                if "814" in orig_name:
                    n_manual814 += 1
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
        
# Runs manual mask on a file where lines are:
# # 'filename' (must be preceded by # and a space)
# list of x-vertices for mask 1
# list of y-vertices for mask 1
# list of x-vertices for mask 2, etc...
def manual_mask_catalogs(catalog_vertex_file):
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

#manual_mask_catalogs("manual_masks.txt")            

def run_pipeline(image,background,final_out_name,make_classifier=False,SG_params=None,diff_spike_params=f606w_spike_params,clean=False,debug=True):
    s = time.time()
    global n_bright
    global n_faint
    global n_c1
    global n_border
    global n_diff
    global n_manual
    global n_hotcold
    print "Running image pipeline on file:", image, "with background:", background
    run_sextractor(image, background, bright_config_dict, output_params, "bright")
    run_sextractor(image, background, faint_config_dict, output_params, "faint")
    bright = asciidata.open("bright.cat")
    n_bright += bright.nrows
    faint = asciidata.open("faint.cat")
    n_faint += faint.nrows
    if debug:
        print "SExtractor finished. Now making segmentation map..."
    make_segmentation_map("bright.cat", image, "seg_map")
    if debug:
        print "Segmentation map created. Now filtering the faint catalog..."
    filter_cat_with_segmentation_map("faint.cat", "seg_map.fits", "faint_filter.cat")
    if debug:
        print "Filtering complete. Now merging bright and faint catalogs..."
    merge_bright_faint_catalogs("bright.cat", "faint_filter.cat", "merged.cat")
    if debug:
        print "Catalogs merged. Now fixing catalog indices..."
    replace_indexes_of_merged_catalog("merged.cat", "merged_fix.cat")
    c1 = asciidata.open("merged_fix.cat")
    n_c1 += c1.nrows
    n_hotcold += (bright.nrows + faint.nrows - c1.nrows)
    if debug:
        print "Indices replaced."
    if make_classifier == True:
        if debug:
            print "Classifier plot generating..."
        make_classifier_plot("merged_fix.cat")
    if SG_params is None:
        if debug:
            print "No parameters for SG classification detected. Using the default settings for SG classification..."
        alter_catalog_for_classification("merged_fix.cat","merged_fix_class.cat")
    else:
        if debug:
            print "Using manual parameters for SG classification..."
        alter_catalog_for_classification("merged_fix.cat","merged_fix_class.cat",SG_params[0],SG_params[1],SG_params[2],SG_params[3])
    if debug:
        print "SG Classification complete. Now adding S/N ratios..."
    make_SNR("merged_fix_class.cat","merged_fix_class_SN.cat")
    if debug:
        print "S/N Ratios added. Now removing edge objects..."
    edge_overlap_clean("merged_fix_class_SN.cat","edge_clean.cat")
    replace_indexes_of_merged_catalog("edge_clean.cat","edge_clean_index.cat")
    border = asciidata.open("edge_clean_index.cat")
    n_border += (c1.nrows - border.nrows)
    if debug:
        print "Edge objects removed. Now cleaning for diffraction spikes..."
    diffraction_mask_cleanup("edge_clean_index.cat", "diffraction_cleanup.cat",diff_spike_params=diff_spike_params)
    replace_indexes_of_merged_catalog("diffraction_cleanup.cat",final_out_name)
    diff = asciidata.open("diffraction_cleanup.cat")
    n_diff += (c1.nrows - diff.nrows)
    if clean:
        print "Cleaning extraneous catalogs..."
        subprocess.call(["rm","bright.cat"])
        subprocess.call(["rm","faint.cat"])
        subprocess.call(["rm","faint_filter.cat"])
        subprocess.call(["rm","merged.cat"])
        subprocess.call(["rm","merged_fix.cat"])
        subprocess.call(["rm","merged_fix_class.cat"])
        subprocess.call(["rm","merged_fix_class_SN.cat"])
        subprocess.call(["rm","edge_clean.cat"])
        subprocess.call(["rm","edge_clean_index.cat"])
        subprocess.call(["rm","diffraction_cleanup.cat"])
    e = time.time()
    t = e-s
    print "Process complete. Time:", t, "seconds. File written out to :", final_out_name


make_classifier_plot("EGS_10134_0i_acs_wfc_f606w_30mas_unrot_drz.fits.cat")
'''
#####F606w pipeline
print "Running image pipeline on COSMOS files in the f606w filter."
s = time.time()    
for i in range(len(f606w_files)):
    run_pipeline(f606w_files[i],f606w_backgrounds[i],f606w_files[i]+".cat",make_classifier=False,clean=True)
    e = time.time()
    print "Total time elapsed is approximately", np.floor((e-s)/60), "minutes", ((e-s)/60-np.floor((e-s)/60))*60, "seconds."
'''
'''
#####F814w pipeline
print "Running image pipeline on COSMOS files."
s = time.time()    
for i in range(len(f814w_files)):
    run_pipeline(f814w_files[i],f814w_backgrounds[i],f814w_files[i].strip()+".cat",make_classifier=False,diff_spike_params=f814w_spike_params,clean=True)
    e = time.time()
    print "Total time elapsed is approximately", np.floor((e-s)/60), "minutes", ((e-s)/60-np.floor((e-s)/60))*60, "seconds."
'''
'''
def classification_map(catalog, out_name):
    xdim = 7500
    ydim = 7500
    #Get position tuples of the detected objects and round them to the nearest integer
    positions = []
    catalog = asciidata.open(catalog)
    for i in range(catalog.nrows):
        x_min = catalog['XMIN_IMAGE'][i]
        y_min = catalog['YMIN_IMAGE'][i]
        x_max = catalog['XMAX_IMAGE'][i]
        y_max = catalog['YMAX_IMAGE'][i]
        is_star = catalog['IS_STAR'][i]
        snr = catalog['SNR'][i]
        pos_tuple = (x_min, x_max, y_min, y_max, is_star,snr)
        positions.append(pos_tuple)
        
    classification_array = np.zeros((xdim,ydim))
    for j in range(len(positions)):
        x_min = (positions[j])[0]
        x_max = (positions[j])[1]
        y_min = (positions[j])[2]
        y_max = (positions[j])[3]
        is_star = (positions[j])[4]
        snr = (positions[j])[5]
        for x in range(x_min,x_max):
            for y in range(y_min,y_max):
                if is_star == 1:
                    classification_array[y,x] = snr
                else:
                    classification_array[y,x] = 0
    
    #Write out to a fits file
	hdu_out = pyfits.PrimaryHDU(classification_array)
    hdu_out.writeto(out_name + ".fits",clobber=True)
'''
'''
for file in f606w_files:
    delete_null(file + ".cat")
    
for file in f814w_files:
    delete_null(file + ".cat")
    delete_overlap(file + ".cat")
    delete_null(file + ".cat")
'''
'''
s = open("statistics_606.txt", "a")
#s.write("Items in c1: " + str(n_c1) + "\n")
#s.write("Items in cold: " + str(n_bright) + "\n")
#s.write("Items in hot: " + str(n_faint) + "\n")
#s.write("Items deleted in hot/cold step: " + str(n_hotcold) + "\n")
#s.write("Items in noisy border: " + str(n_border) + "\n")
#s.write("Items in a star diffraction mask: " + str(n_diff) + "\n")
s.write("Items in a manual mask: " + str(n_manual606))
s.close()

s2 = open("statistics.txt", "a")
s2.write("Items in a manual mask: " + str(n_manual814))
s2.close()      
'''    
    
    
