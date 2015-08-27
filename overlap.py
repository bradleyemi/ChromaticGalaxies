import asciidata
import subprocess
import pyfits
import numpy as np
import matplotlib.pyplot as plt
import time

image_map = {}
check_map = {}

#Takes an index of a catalog and determines which indexes to check (returns a list).
#This function will need to be modified for non AEGIS data. It just needs to take in an image number,
#and spit out a list of the indices of neighboring tiles.
def index_to_check(n):
   out = []
   if n % 3 != 2: #if a tile is not on the right side of the panel
       out.append(n+1) #Check the tile on the right
   if n % 3 != 0 and n < 60:
       out.append(n+2) #Check the tile on the top-left
   if n < 60: #If a tile is not on the top of the panel
       out.append(n+3) #Check the tile on the top
   if n % 3 != 2 and n < 60:
       out.append(n+4) #Check the tile on the top-right
   return out
   
def delete_items(catalog, delete_numbers):
    f = open(catalog, "r")
    lines = f.readlines()
    f.close()
    f = open(catalog, "w")
    for line in lines:
        split = line.split()
        if split[0] == "#":
            f.write(line)
        else:
            if int(split[0]) not in delete_numbers:
                f.write(line)
    f.close()
    
def check_adjacent(base_catalog, all_catalogs, check_indices, tolerance = 1./18000):
    nDeleted = 0
    base = asciidata.open(base_catalog)
    base_alpha = [base['ALPHA_SKY'][i] for i in range(base.nrows)]
    base_delta = [base['DELTA_SKY'][i] for i in range(base.nrows)]
    for index in check_indices:
        nDeletedInd = 0
        print "Checking for overlaps in the base catalog", base_catalog
        print "Checking against catalog", all_catalogs[index]
        check = asciidata.open(all_catalogs[index])
        check_alpha = [check['ALPHA_SKY'][i] for i in range(check.nrows)]
        check_delta = [check['DELTA_SKY'][i] for i in range(check.nrows)]
        check_number = [check['NUMBER'][i] for i in range(check.nrows)]
        delete_numbers = []
        for j in range(check.nrows):
            if j % 1000 == 0:
                print "Checking item", j
            item_alpha = check_alpha[j]
            item_delta = check_delta[j]
            item_bools = [(abs(item_alpha - base_alpha[k]) < tolerance and abs(item_delta - base_delta[k]) < tolerance) for k in range(base.nrows)]
            if True in item_bools:
                delete_numbers.append(check_number[j])
                nDeleted += 1
                nDeletedInd += 1
        print nDeletedInd, "objects were deleted in this check."
        delete_items(all_catalogs[index], delete_numbers)
    return nDeleted

def overlap_all(catalogs):
    nDeleted = 0
    s = time.time()
    for i in range(len(catalogs)):
        e = time.time()
        print "\n \n \nChecking catalog", i
        print "Time elapsed is", e-s, "seconds. \n \n \n"
        base_catalog = catalogs[i]
        n = check_adjacent(base_catalog, catalogs, index_to_check(i))
        nDeleted += n
        return nDeleted


