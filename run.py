import HST_Sextractor_new
import overlap
import assoc_catalogs
import postage_stamps

### Modify the filenames and paths below to point to your images, etc. ###
### Required inputs are described below. ##

#No functionality for more than 2 filters right now.
n_filters = 2

image_file = ["f606w_filenames.txt", "f814w_filenames.txt"] #Text file listing the names of the images, in the same order for both filters. 
background_file = ["f606w_backgrounds.txt", "f814w_backgrounds.txt"] #Text file listing the names of the inverse weight files, in the same order for both filters.
filter = [606, 814] #The corresponding filters. (only supports f814w and f606w, right now)
manual_mask_file = [None, None] #A manual masking file as described in HST_Sextractor.py. 
out_name = ["f606w_example", "f814w_example"] #The outname prefixes of the intermediary catalogs, in case you want to check if each step is working.
tt_root_dir = ["/Users/bemi/JPL/", "/Users/bemi/JPL"] #Where the TT files are stored.
tt_star_file = [tt_root_dir[0] + "F606W_TT/606_stars.txt", tt_root_dir[1] + "F814W_TT/814_stars.txt"] #The full path name of where the TT centroid lists are. (Text file, x-centroid and y-centroid are the columns)
catalog_list_file = ["f606w_catalogs.txt", "f814w_catalogs.txt"] #The name of the text file generated listing the names of all the catalogs.
postage_stamp_path = "/Users/bemi/" #Where the script will put the postage stamps.

################### No need to modify code below this line for most users. ############################

def get_focus_catalogs(image_file, background_file, filter, manual_mask_file, out_name, tt_root_dir, tt_star_file, catalog_list_file):
    #file import
    f = open(image_file)
    g = open(background_file)
    files = []
    backgrounds = []
    for line in f.readlines():
        files.append(line.strip())
    for line in g.readlines():
        backgrounds.append(line.strip())    
        
    #catalog generation and cleaning

    #you can adjust whether or not you only want to do the first n files here. To do them all, just leave it as len(files). 
    #n = len(files)
    n = 1

    catalogs = []
    f = open(catalog_list_file, "w")
    for i in range(n):
        #You should modify the names of the catalogs so each catalog has a unique identifier. This could just be the filename of the image. 
        #out_cat_name = files[i]
        out_cat_name = str(filter) + "_" + (files[i])[10:12]
        cat = HST_Sextractor_new.GalaxyCatalog(files[i], backgrounds[i], filter, out_cat_name, manual_mask_file)
        f.write(out_cat_name[:len(out_cat_name)-4] + ".focus.cat")
        cat.generate_catalog()
        catalogs.append(cat)
    
    #adding focus positions    
    cat_list = HST_Sextractor_new.GalaxyCatalogList(catalogs, filter)
    cat_list.add_focus(out_name + "_focus.txt", tt_root_dir)

    #deleting duplicates from overlap
    #nDeleted = overlap.overlap_all(cat_list.catalogs)
    #print nDeleted, "objects deleted"

for i in range(n_filters):
    get_focus_catalogs(image_file[i], background_file[i], filter[i], manual_mask_file[i], out_name[i], tt_root_dir[i], tt_star_file[i], catalog_list_file[i])

#associate catalogs
catalogs_1 = []
catalogs_2 = []
f = open(catalog_list_file[0])
g = open(catalog_list_file[1])
for line in f.readlines():
    catalogs_1.append(line.strip())
for line in g.readlines():
    catalogs_2.append(line.strip())
f.close()
g.close()

for i in range(len(catalogs_1)):
    s = time.time()
    print "Associating catalog", i
    assoc_catalogs.assoc_catalogs_opt(catalogs_1[i], catalogs_2[i])
    print "Time:", time.time()-s

#get postage stamps
for i in range(n_filters):
    postage_stamps.get_postage_stamps_all(catalog_list_file[i], image_file[i], filter[i], postage_stamp_path)
    



