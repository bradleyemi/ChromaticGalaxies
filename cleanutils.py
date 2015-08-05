'''
This file needs to be imported for HST_Sextractor.py to run.
'''

import asciidata
import subprocess
import pyfits
import numpy as np
import matplotlib.pyplot as plt
import time

### Helper functions for HST_Sextractor.py ###

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