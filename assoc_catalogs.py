import asciidata
import time
import numpy as np
import matplotlib.pyplot as plt

def assoc_catalogs(c1, c2, tolerance = 1/18000.):
    cat1 = asciidata.open(c1)
    cat2 = asciidata.open(c2)
    for k in range(cat2.nrows):
        cat2['ASSOC'][k] = -1
    for i in range(cat1.nrows):
        if i % 500 == 0:
            print "Associating catalog item", i
        alpha1 = cat1['ALPHA_SKY'][i] 
        delta1 = cat1['DELTA_SKY'][i]
        cat1['ASSOC'][i] = -1
        for j in range(cat2.nrows):
            alpha2 = cat2['ALPHA_SKY'][j]
            delta2 = cat2['DELTA_SKY'][j]
            if abs(alpha1-alpha2) < tolerance and abs(delta1-delta2) < tolerance:
                cat1['ASSOC'][i] = i
                cat2['ASSOC'][j] = i
                break
    cat1.writeto(c1)
    cat2.writeto(c2)

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
    return mag_tuples, flux_tuples

test_catalog1 = "EGS_10134_01_acs_wfc_f606w_30mas_unrot_drz.fits.focus.cat"
test_catalog2 = "EGS_10134_01_acs_wfc_f814w_30mas_unrot_drz.fits.focus.cat"


catalogs_606 = []
catalogs_814 = []
f = open("f606w_catalogs.txt")
g = open("f814w_catalogs.txt")
for line in f.readlines():
    catalogs_606.append(line.strip())
for line in g.readlines():
    catalogs_814.append(line.strip())
f.close()
g.close()

mags1 = []
mags2 = []
fluxr1 = []
fluxr2 = []
for i in range(len(catalogs_606)):
    s = time.time()
    print "Associating catalog", i
    mag_tuples, flux_tuples = assoc_catalogs_opt(catalogs_606[i], catalogs_814[i])
    for j in range(len(mag_tuples)):
        mags1.append((mag_tuples[j])[0] + 26.508)
        mags2.append((mag_tuples[j])[1] + 25.955)
        fluxr1.append((flux_tuples[j])[0])
        fluxr2.append((flux_tuples[j])[1])
    print "Time:", time.time()-s

f = open("filter_statistics.txt", "w")
for i in range(len(mags1)):
    f.write(str(mags1[i]) + " " + str(mags2[i]) + " " + str(fluxr1[i]) + " " + str(fluxr2[i]) + "\n")
f.close()

'''
plt.scatter(mags1,mags2)
plt.xlabel("f606w magnitude")
plt.ylabel("f814w magnitude")
plt.title("Magnitudes for galaxies in f606w and f814w")
plt.show()
plt.scatter(fluxr1,fluxr2)
plt.xlabel("f606w half-light radius")
plt.ylabel("f814w half-light radius")
plt.title("Half-light radius for galaxies in f606w and f814w")
plt.show()
'''

'''
def check_retention(c):
     total = 0
     retained = 0
     f = open(c)
     for line in f.readlines():
         if line[0] == "#":
             continue
         else:
             split = line.split()
             total += 1
             if int(split[22]) != -1:
                 retained += 1
     return retained, total
    
total606 = 0
ret606 = 0
total814 = 0
ret814 = 0

for i in range(len(catalogs_606)):
    s = time.time()
    r1, t1 = check_retention(catalogs_606[i])
    r2, t2 = check_retention(catalogs_814[i])
    total606 += t1
    ret606 += r1
    total814 += t2
    ret814 += r2
    e = time.time()
    print "Time:", np.floor(e-s)

print "606 total", total606
print "814 total", total814
print "606 retained", ret606
print "814 retained", ret814

print "606 percent retained", float(ret606)/float(total606)
print "814 percent retained", float(ret814)/float(total814)
'''   
