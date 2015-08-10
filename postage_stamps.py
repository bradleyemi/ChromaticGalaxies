import galsim
import numpy as np
import pyfits
import asciidata
import matplotlib.pyplot as plt
import time
import subprocess

f606w_catalog_file = open("f606w_catalogs.txt")
f814w_catalog_file = open("f814w_catalogs.txt")

test_catalog = "EGS_10134_01_acs_wfc_f606w_30mas_unrot_drz.fits.focus.cat"
test_image = "EGS_10134_01_acs_wfc_f606w_30mas_unrot_drz.fits"

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
###### End .fits to GalSim import ####


class CatalogObject:
   x_axis_length = 7500
   y_axis_length = 7500
   
   def __init__(self, ident, fname, fweight, filter, x, y, ra, dec, radius, mag, focus, is_star, snr):
       self.ident = ident
       self.fname = fname
       self.fweight = fweight
       self.filter = filter
       self.x = x
       self.y = y
       self.ra = ra
       self.dec = dec
       self.stampL = 4.0*radius #11*np.sqrt((1.5*radius)**2 + (1.2/(0.03*2.35))**2)
       self.mag = mag
       self.focus = focus
       self.is_star = is_star
       self.snr = snr
       self.leftBound = int(self.x - self.stampL)
       self.rightBound = int(self.x + self.stampL)
       self.topBound = int(self.y + self.stampL)
       self.bottomBound = int(self.y - self.stampL)      
    
   def is_within_image(self):
       return self.x - self.stampL*0.5 > 0 and \
       self.x + self.stampL*0.5 < self.x_axis_length and \
       self.y - self.stampL*0.5 > 0 and \
       self.y + self.stampL*0.5 < self.y_axis_length
     
   def postage_stamp(self, image_data=None):
       if not self.is_within_image():
           raise ValueError("The postage stamp is outside the CCD boundary.")
       else:
           if image_data == None:
               f = pyfits.open(fname)
               image_data = f[0].data
               f.close()
           im = galsim.Image(image_data)
           boundDifference = (self.rightBound - self.leftBound) - (self.topBound - self.bottomBound)
           if boundDifference == 0:
              b = galsim.BoundsI(self.leftBound, self.rightBound, self.bottomBound, self.topBound)
           if boundDifference == 1:
              b = galsim.BoundsI(self.leftBound, self.rightBound, self.bottomBound, self.topBound+1)
           if boundDifference == -1:
              b = galsim.BoundsI(self.leftBound, self.rightBound+1, self.bottomBound, self.topBound)
           stamp = im.subImage(b)
           return stamp
           
   def get_TT_field(self):
        if self.filter == 606:
            tt_dict = tt_606
        elif self.filter == 814:
            tt_dict = tt_814
        else:
            raise Exception('Invalid filter')
        focus = int(np.round(self.focus))
        return tt_dict[focus]
        
   def find_nearest_centroid(self, tt_stars):
        x0 = self.x
        y0 = self.y
        if self.filter == 606:
            match_dist = 400.
        if self.filter == 814:
            match_dist = 400.
        f = open(tt_stars)
        lines = f.readlines()
        f.close()
        for line in lines:
            split = line.split()
            if split[0] == "#":
                continue
            x = np.float32(split[0])
            y = np.float32(split[1])
            if abs(x-x0) < match_dist and abs(y-y0) < match_dist:
                return (x,y)
        print "No star found"
    
   def PSF(self):
        tt_field = self.get_TT_field()
        if self.filter == 606:
            tt_stars = "606_stars.txt"
        if self.filter == 814:
            tt_stars = "/Users/bemi/JPL/F814W_TT/TinyTim_f-1.stars.dat"
        print "Matching centroid."
        s = time.time()
        (self.x_tt, self.y_tt) = self.find_nearest_centroid(tt_stars)
        print "Centroid found. Time =", time.time()-s
        self.left_tt = int(self.x_tt - self.stampL)
        self.right_tt = int(self.x_tt + self.stampL)
        self.bottom_tt = int(self.y_tt - self.stampL)
        self.top_tt = int(self.y_tt + self.stampL)
        boundDifference = (self.right_tt - self.left_tt) - (self.top_tt - self.bottom_tt)
        if boundDifference == 0:
              b = galsim.BoundsI(self.left_tt, self.right_tt, self.bottom_tt, self.top_tt)
        if boundDifference == 1:
              b = galsim.BoundsI(self.left_tt, self.right_tt, self.bottom_tt, self.top_tt+1)
        if boundDifference == -1:
              b = galsim.BoundsI(self.left_tt, self.right_tt+1, self.bottom_tt, self.top_tt)
        print "Getting tt-field"
        s = time.time()
        tt_file = pyfits.open(tt_field)
        tt_data = tt_file[0].data
        tt_image = galsim.Image(tt_data)
        tt_file.close()
        print "Image imported. Time=", time.time()-s
        sub = tt_image.subImage(b)  
        del tt_image
        return sub      

def snr_hist(catalog):
    cat = asciidata.open(catalog)
    snrs = []
    for i in range(cat.nrows):
        if cat['MAG_AUTO'][i] + 21.1 <= 22.5:
            snrs.append(cat['SNR'][i])
    for i in range(10):
        print i*10, "percentile is", np.percentile(snrs, i*10)
    plt.hist(snrs, bins=50, range=(0,50))
    plt.show()
    
#snr_hist(test_catalog)

def get_postage_stamps(catalog, file, filter, out_name):
    cat = asciidata.open(catalog)
    print "Catalog opened."
    f = pyfits.open(file)
    print "File opened."
    image_data = f[0].data
    f.close()
    image_hdulist = pyfits.HDUList()
    psf_hdulist = pyfits.HDUList()
    snr_list = []
    nSets = 0
    nTotal = 0
    for i in range(cat.nrows):
        print "Adding object", i
        ident = cat['NUMBER'][i]
        fname = cat['FILENAME'][i]
        fweight = fname[:len(fname)-8] + "wht.fits"
        x = cat['X_IMAGE'][i]
        y = cat['Y_IMAGE'][i]
        ra = cat['ALPHA_SKY'][i]
        dec = cat['DELTA_SKY'][i]
        radius = cat['FLUX_RADIUS'][i]
        mag = cat['MAG_AUTO'][i]
        focus = cat['FOCUS'][i]
        is_star = cat['IS_STAR'][i]
        snr = cat['SNR'][i]
        snr_list.append(snr)
        object = CatalogObject(ident, fname, fweight, filter, x, y, ra, dec, radius, mag, focus, is_star, snr)
        if object.stampL < 500 and \
           object.is_within_image() and \
           object.is_star == 0 and \
           object.snr > 20.0 and \
           object.mag + 21.1 < 22.5:
            try:
                img = object.postage_stamp(image_data=image_data)
                psf = object.PSF()
            except:
                continue
            img_data = img.array
            psf_data = psf.array
            if len(image_hdulist) == 0 and len(psf_hdulist) == 0:
                image_hdulist.append(pyfits.PrimaryHDU(data=img_data))
                psf_hdulist.append(pyfits.PrimaryHDU(data=psf_data))
            else:
                image_hdulist.append(pyfits.ImageHDU(data=img_data))
                psf_hdulist.append(pyfits.ImageHDU(data=psf_data))
            del img
            del psf
            del object
            nTotal += 1
        else:
            print "Object skipped."
        if len(image_hdulist) == 1:
            try: 
                 galsim.fits.writeFile(str(nSets) + ".fits", image_hdulist, dir="/Users/bemi/JPL/" + out_name)
                 galsim.fits.writeFile(str(nSets) + ".psf.fits", psf_hdulist, dir="/Users/bemi/JPL/" + out_name)
            except:
                 subprocess.call(["mkdir", "/Users/bemi/JPL/" + out_name])
                 galsim.fits.writeFile(str(nSets) + ".fits", image_hdulist, dir="/Users/bemi/JPL/" + out_name)
                 galsim.fits.writeFile(str(nSets) + ".psf.fits", psf_hdulist, dir="/Users/bemi/JPL/" + out_name)
            del image_hdulist
            del psf_hdulist
            image_hdulist = pyfits.HDUList()
            psf_hdulist = pyfits.HDUList()
            nSets += 1
    print "total objects counted", nTotal
        
def get_postage_stamps_all(catalog_list_file, image_list_file, filter):
    f = open(catalog_list_file)
    lines = f.readlines()
    g = open(image_list_file)
    image_lines = g.readlines()
    f.close()
    g.close()
    for i in range(len(lines)):
        get_postage_stamps(lines[i].strip(), image_lines[i].strip(), filter, str(filter) + "_" + (image_lines[i])[10:12])
    

get_postage_stamps_all("f814w_catalogs.txt", "f814w_filenames.txt", 814)
    
                