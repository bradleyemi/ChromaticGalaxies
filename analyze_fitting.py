import pyfits
import asciidata
import numpy as np
import matplotlib.pyplot as plt

hdulist = pyfits.open("/Users/bemi/JPL/assoc_606_01/output/RAWFIT00000.00471.fits")
hdulist2 = pyfits.open("/Users/bemi/JPL/assoc_814_01/output/RAWFIT00000.00471.fits")
stats = open("filter_statistics.txt")
data1 = hdulist[1].data
cols1 = hdulist[1].columns
data2 = hdulist2[1].data
cols2 = hdulist2[1].columns
#cols.info()

#Magnitudes and half-light-radii

m1 = []
m2 = []
r1 = []
r2 = []
for line in stats.readlines():
    s = line.split()
    m1.append(float(s[0]))
    m2.append(float(s[1]))
    r1.append(float(s[2]))
    r2.append(float(s[3]))

m1 = np.asarray(m1)
m2 = np.asarray(m2)
r1 = np.asarray(r1)*0.03
r2 = np.asarray(r2)*0.03
mresid = m1-m2
rresid = r1-r2
print np.mean(rresid)

plt.hist([m1,m2],label=["f606w","f814w"], bins=30)
plt.title("Magnitudes of galaxies in final sample")
plt.xlabel("Magnitude")
plt.ylabel("Frequency")
plt.legend()
plt.show()

plt.hist(mresid, bins=30)
plt.title("Difference in magnitude for individual galaxies across f606w and f814w filters")
plt.xlabel("Residuals, f606w magnitude - f814w magnitude")
plt.ylabel("Frequency")
plt.show()

plt.hist([r1,r2],label=["f606w","f814w"], bins=30)
plt.title("Half-light radius of galaxies in final sample")
plt.xlabel("Half-light radius, arcseconds")
plt.ylabel("Frequency")
plt.legend()
plt.show()

plt.hist(rresid, bins=30, range=(-0.3,0.3))
plt.title("Difference in half light radius for individual galaxies across f606w and f814w filters")
plt.xlabel("Residuals, f606w radius - f814w radius, arcseconds")
plt.ylabel("Frequency")
plt.show()


#Fits and chi-square
fit_dvc = (data1.field('FIT_DVC'), data2.field('FIT_DVC'))
chisq_dvc = (data1.field('CHISQ_DVC'), data2.field('CHISQ_DVC'))

fit_ser = (data1.field('FIT_SER'), data2.field('FIT_SER'))
chisq_ser = (data1.field('CHISQ_SER'), data2.field('CHISQ_SER'))

fit_expdvc = (data1.field('FIT_EXPDVC'), data2.field('FIT_EXPDVC'))
chisq_expdvc = (data1.field('CHISQ_EXPDVC'), data2.field('CHISQ_EXPDVC'))

plt.hist([chisq_dvc[0], chisq_ser[0], chisq_expdvc[0]], bins=10, range=(0.0,3.0), label=["De Vacoleurs", "Sersic", "Bulge-disk"])
plt.xlabel("Chi-square value")
plt.ylabel("Frequency")
plt.legend()
plt.title("Goodness of fit for galaxy profiles in f606w filter")
plt.show()
plt.hist([chisq_dvc[1], chisq_ser[1], chisq_expdvc[1]], bins=10, range=(0.0,3.0), label=["De Vacoleurs", "Sersic", "Bulge-disk"])
plt.xlabel("Chi-square value")
plt.ylabel("Frequency")
plt.legend()
plt.title("Goodness of fit for galaxy profiles in f814w filter")
plt.show() 

#Flux ratio for b/d

fluxratio_expdvc = (data1.field('FLUXRATIO_EXPDVC'), data2.field('FLUXRATIO_EXPDVC'))
plt.hist(fluxratio_expdvc, bins=20, range=(0.0,1.0), label=("f606w", "f814w"))
plt.title("Bulge-to-disk ratios from best bulge-disk fit")
plt.xlabel("Bulge-to-disk ratio")
plt.ylabel("Frequency")
plt.legend()
plt.show()


#Sersic index

fit_ser606 = fit_ser[0]
fit_ser814 = fit_ser[1]
ser_index606 = [(fit_ser606[i])[2] for i in range(len(fit_ser606))]
ser_index814 = [(fit_ser814[i])[2] for i in range(len(fit_ser814))]
plt.hist([ser_index606, ser_index814], bins=20, range=(0.0, 10.0), label=["f606w", "f814w"])
plt.title("Best-fit Sersic index for fitting galaxy profiles")
plt.xlabel("Sersic index")
plt.ylabel("Frequency")
plt.legend()
plt.show()


#Fitting status check

name = (data1.field('NAME'), data2.field('NAME'))
stat_dvc = (data1.field('STAT_DVC'), data2.field('STAT_DVC'))
stat_ser = (data1.field('STAT_SER'), data2.field('STAT_SER'))
stat_expdvc = (data1.field('STAT_EXPDVC'), data2.field('STAT_EXPDVC'))
for i in range(np.size(name[0])):
  for j in range(2):
    if (stat_dvc[j])[i] == 0 or (stat_dvc[j])[i] == 5:
        print j, (name[j])[i], "dvc"
    if (stat_ser[j])[i] == 0 or (stat_ser[j])[i] == 5:
        print j, (name[j])[i], "ser"
    if (stat_expdvc[j])[i] == 0 or (stat_expdvc[j])[i] == 5:
        print j, (name[j])[i], "expdvc"
