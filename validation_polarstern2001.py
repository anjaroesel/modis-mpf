from scipy import *
import sys,pipes,struct,os,glob,fnmatch,pickle
from pylab import *
from scipy.stats import nanmean, nanstd
from mpl_toolkits.basemap import Basemap
from polar_projection import *
import scipy.io as io
from mpl_toolkits.basemap import Basemap
import Nio
from sattools import *
from gmttools import *


#validation_polarstern2001.py

#local
#folder='/scratch/local3/Ship_Observ/'
#snow
folder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/SHIP_OBS/'
daylyfile=folder+'POLsternARK172_2001_DAYLY.csv'

###read ship data 2001
a=loadtxt(daylyfile,skiprows=1)
a[a==999]=nan
day=a[:,0]
lat=a[:,1]
lon=a[:,2]
mpPS2001=a[:,3]
mpstd=a[:,4]
icek=a[:,5]
icekstd=a[:,6]

x,y=mapll(lat,lon,1)




#read mp data
####---snow1:
folder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/PRODUCTS/'
plotfolder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/PRODUCTS/plots/'
dataset0=folder+'2001_209/2001_209_rmp_125masked.nc'
dataset1=folder+'2001_217/2001_217_rmp_125masked.nc' #no data available
dataset2=folder+'2001_225/2001_225_rmp_125masked.nc'
dataset3=folder+'2001_233/2001_233_rmp_125masked.nc'
dataset4=folder+'2001_241/2001_241_rmp_125masked.nc'
dataset5=folder+'2001_249/2001_249_rmp_125masked.nc'
dataset6=folder+'2001_257/2001_257_rmp_125masked.nc'

date0=dataset0[-17:-9]
date1=dataset1[-17:-9]#e.g.'2008_169'
print date1
date2=dataset2[-17:-9]
date3=dataset3[-17:-9]
date4=dataset4[-17:-9]
date5=dataset5[-17:-9]
date6=dataset6[-17:-9]



jahr=date0[0:4]
julday0=date0[5:8]
julday1=date1[5:8]
julday2=date0[5:8]
julday3=date0[5:8]
julday4=date0[5:8]
julday5=date0[5:8]
julday6=date0[5:8]

MP_obj0=NetCDFFile(dataset0)
MP_obj1=NetCDFFile(dataset1)
MP_obj2=NetCDFFile(dataset2)
MP_obj3=NetCDFFile(dataset3)
MP_obj4=NetCDFFile(dataset4)
MP_obj5=NetCDFFile(dataset5)
MP_obj6=NetCDFFile(dataset6)
#MP_SD=NetCDFFile(dataset_sd)
#MP_FL=NetCDFFile(dataset_flag)

MPX=array(MP_obj0.variables['x'][:])
MPY=array(MP_obj0.variables['y'][:])
MP0=array(MP_obj0.variables['z'][:])
MP1=array(MP_obj1.variables['z'][:])
MP2=array(MP_obj2.variables['z'][:])
MP3=array(MP_obj3.variables['z'][:])
MP4=array(MP_obj4.variables['z'][:])
MP5=array(MP_obj5.variables['z'][:])
MP6=array(MP_obj6.variables['z'][:])

#STDEV=array(MP_SD.variables['z'][:])
#FLAG=array(MP_FL.variables['z'][:])



###getting indices of MPset
iX=[]
iY=[]
for n in range(len(y)):
    print n
    iyy=nonzero(MPY>=y[n])[0]
    iyy=iyy[0]
    iY.append(iyy)

for n in range(len(x)):
    print n
    ixx=nonzero(MPX>=x[n])[0]
    ixx=ixx[0]
    iX.append(ixx)

#I=array((iY,iX))


MPval=[]
for n in range(len(x)):
    print n
    mpval=MP6[iY[n],iX[n]]
    MPval.append(mpval)


MPval[41:]


MPval=array([0.27984005, 0.33492154, 0.32458195, 0.30340204, 0.28097454, 0.2706652, 0.21701634, 0.023085214, 0.21095981, 0.12607606], dtype=float32)


mpPS2001=array([ 24.41,  22.92,  24.  ,   9.38,   3.88,   6.71,  10.42,  20.56,   1.17,  16.25,])




MPval=MPval*100

#from scipy.stats import linregress
cc=corrcoef([MPval,mpPS2001])
print cc
RMSE=sqrt((1./n)*sum((MPval-mpPS2001)**2))
print RMSE

gx=linspace(0.,35,100)

figure()
plot(MPval,mpPS2001,'k+')
title('mp MODIS 0.5km grid vs Ship Observation PS 2001')
ylabel('MP from POLARSTERN [%]')
xlabel('MP from MODIS [%]')
#show()
m=1.
b=0.
gy=m*gx+b
plot(gx,gy,'gray')



annotate('RMSE='+str(round(RMSE,3)),xy=(1,30))
#annotate('corr='+str(round(cc[1][0],3)),xy=(1,30))
savefig('/pf/u/u241127/plots/validation_ps2001.png')








