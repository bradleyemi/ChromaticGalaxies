import asciidata
import time
import numpy as np
import matplotlib.pyplot as plt

def assoc_catalogs_opt(c1, c2, tolerance = 1/18000.):
    mag_tuples = []
    flux_tuples = []
    cat1 = asciidata.open(c1)
    cat2 = asciidata.open(c2)
    for i in range(cat2.nrows):
        cat2['ASSOC'][i] = -1
    for i in range(cat1.nrows):
        cat1['ASSOC'][i] = -1
    alpha1 = [cat1['ALPHA_SKY'][i] for i in range(cat1.nrows)]
    alpha2 = [cat2['ALPHA_SKY'][i] for i in range(cat2.nrows)]
    delta1 =   [cat1['DELTA_SKY'][i] for i in range(cat1.nrows)]
    delta2 =   [cat2['DELTA_SKY'][i] for i in range(cat2.nrows)]
    for i in range(cat1.nrows):
        match = [abs(alpha1[i]-alpha2[j]) < tolerance and abs(delta1[i]-delta2[j]) < tolerance for j in range(cat2.nrows)]
        try:
           idx = match.index(True)
           if cat2['ASSOC'][idx] == -1 \
           and cat1['SNR'][i] > 20 \
           and cat2['SNR'][idx] > 20 \
           and cat1['MAG_AUTO'][i] + 21.1 < 22.5 \
           and cat2['MAG_AUTO'][idx] + 21.1 < 22.5 \
           and cat1['FLUX_RADIUS'][i] > 0 \
           and cat2['FLUX_RADIUS'][idx] > 0 \
           and cat1['IS_STAR'][i] == 0 \
           and cat2['IS_STAR'][idx] == 0 \
           and cat1['FLUX_RADIUS'][i] < 500 \
           and cat2['FLUX_RADIUS'][idx] < 500:
               cat1['ASSOC'][i] = i
               cat2['ASSOC'][idx] = i
               mag_tuples.append((cat1['MAG_AUTO'][i], cat2['MAG_AUTO'][idx]))
               flux_tuples.append((cat1['FLUX_RADIUS'][i],cat2['FLUX_RADIUS'][idx]))               
        except ValueError:
           continue
    #cat1.writeto(c1)
    #cat2.writeto(c2)

