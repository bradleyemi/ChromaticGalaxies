import galsim
import numpy as np
import pyfits
import asciidata
import matplotlib.pyplot as plt
import time
import subprocess

class CatalogObject:
   x_axis_length = 7500
   y_axis_length = 7500
   
   def __init__(self, ident, fname, fweight, filter, x, y, ra, dec, radius, mag, focus, is_star, snr, assoc):
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
       self.assoc = assoc
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
               f = pyfits.open(self.fname)
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
        
   def get_TT_field(self, tt_dict):
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
    
   def PSF(self, tt_dict, tt_stars):
        tt_field = self.get_TT_field(tt_dict)
        (self.x_tt, self.y_tt) = self.find_nearest_centroid(tt_stars)
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

def get_postage_stamps(catalog, file, weight, filter, out_name, out_path):
    cat = asciidata.open(catalog)
    print "Catalog opened."
    f = pyfits.open(file)
    g = pyfits.open(weight)
    print "File opened."
    image_data = f[0].data
    weight_data = g[0].data
    f.close()
    g.close()
    image_hdulist = pyfits.HDUList()
    weight_hdulist = pyfits.HDUList()
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
        assoc = cat['ASSOC'][i]
        snr_list.append(snr)
        if assoc == -1:
            continue
        if assoc is None:
            continue
        object = CatalogObject(ident, fname, fweight, filter, x, y, ra, dec, radius, mag, focus, is_star, snr, assoc)
        if object.assoc != -1 and \
           object.is_within_image():
            try:
                img = object.postage_stamp(image_data=image_data)
                weight = object.postage_stamp(image_data=weight_data)
                psf = object.PSF()
            except:
                continue
            img_data = img.array
            wht_data = weight.array 
            psf_data = psf.array
            if len(image_hdulist) == 0 and len(psf_hdulist) == 0:
                image_hdulist.append(pyfits.PrimaryHDU(data=img_data))
                weight_hdulist.append(pyfits.PrimaryHDU(data=wht_data))
                psf_hdulist.append(pyfits.PrimaryHDU(data=psf_data))
            else:
                image_hdulist.append(pyfits.ImageHDU(data=img_data))
                weight_hdulist.append(pyfits.ImageHDU(data=wht_data))
                psf_hdulist.append(pyfits.ImageHDU(data=psf_data))
            del img
            del weight
            del psf
            del object
            nTotal += 1
        else:
            print "Object skipped."
        if len(image_hdulist) == 1:
            try: 
                 galsim.fits.writeFile(str(assoc) + ".0_" + str(ra) + "_" + str(dec) + ".processed.fits", image_hdulist, dir=out_path + out_name + "/images/")
                 galsim.fits.writeFile(str(assoc) + ".0_" + str(ra) + "_" + str(dec) + ".wht.fits", weight_hdulist, dir=out_path + out_name + "/ivar/")
                 galsim.fits.writeFile(str(assoc) + ".0_" + str(ra) + "_" + str(dec) + ".psf.fits", psf_hdulist, dir=out_path + out_name + "/psf/")
            except:
                 subprocess.call(["mkdir", out_path + out_name])
                 subprocess.call(["mkdir", out_path + out_name + "/images/"])
                 subprocess.call(["mkdir", out_path + out_name + "/ivar/"])
                 subprocess.call(["mkdir", out_path + out_name + "/psf/"])
                 galsim.fits.writeFile(str(assoc) + ".0_" + str(ra) + "_" + str(dec) + ".processed.fits", image_hdulist, dir=out_path + out_name + "/images/")
                 galsim.fits.writeFile(str(assoc) + ".0_" + str(ra) + "_" + str(dec) + ".wht.fits", weight_hdulist, dir=out_path+ out_name + "/ivar/")
                 galsim.fits.writeFile(str(assoc) + ".0_" + str(ra) + "_" + str(dec) + ".psf.fits", psf_hdulist, dir=out_path + out_name + "/psf/")
            del image_hdulist
            del weight_hdulist
            del psf_hdulist
            image_hdulist = pyfits.HDUList()
            psf_hdulist = pyfits.HDUList()
            weight_hdulist = pyfits.HDUList()
            nSets += 1
    print "total objects counted", nTotal
        
def get_postage_stamps_all(catalog_list_file, image_list_file, filter, out_path, start=0):
    f = open(catalog_list_file)
    lines = f.readlines()
    g = open(image_list_file)
    image_lines = g.readlines()
    f.close()
    g.close()
    for i in range(start,len(lines)):
        get_postage_stamps(lines[i].strip(), image_lines[i].strip(), (image_lines[i].strip())[:len(image_lines[i].strip())-8] + "wht.fits", filter, "stamps_" + image_lines[i].strip(), out_path)
    
                