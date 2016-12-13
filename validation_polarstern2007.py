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


#run validation_polarstern2007.py


#snow
folder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/SHIP_OBS/'
daylyfile=folder+'POLsternARK222_2007_DAYLY.csv'

###read ship data 2007
a=loadtxt(daylyfile,skiprows=1)
a[a==999]=nan
day=a[:,0]
lat=a[:,1]
lon=a[:,2]
mpPS2007=a[:,3]
mpstd=a[:,4]
icek=a[:,5]
icekstd=a[:,6]

x,y=mapll(lat,lon,1)




#read mp data
####---snow1:
folder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/PRODUCTS/'
plotfolder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/PRODUCTS/plots/'

dataset0=folder+'2007_201/2007_201_rmp_125masked.nc'
dataset1=folder+'2007_209/2007_209_rmp_125masked.nc'
dataset2=folder+'2007_217/2007_217_rmp_125masked.nc'
dataset3=folder+'2007_225/2007_225_rmp_125masked.nc'
dataset4=folder+'2007_233/2007_233_rmp_125masked.nc'
dataset5=folder+'2007_241/2007_241_rmp_125masked.nc'
dataset6=folder+'2007_249/2007_249_rmp_125masked.nc'
dataset7=folder+'2007_257/2007_257_rmp_125masked.nc'

date0=dataset0[-17:-9]
date1=dataset1[-17:-9]#e.g.'2008_169'
print date1
date2=dataset2[-17:-9]
date3=dataset3[-17:-9]
date4=dataset4[-17:-9]
date5=dataset5[-17:-9]
date6=dataset6[-17:-9]



jahr=date1[0:4]
julday1=date1[5:8]
julday2=date1[5:8]
julday3=date1[5:8]
julday4=date1[5:8]
julday5=date1[5:8]
julday6=date1[5:8]

MP_obj0=NetCDFFile(dataset0)
MP_obj1=NetCDFFile(dataset1)
MP_obj2=NetCDFFile(dataset2)
MP_obj3=NetCDFFile(dataset3)
MP_obj4=NetCDFFile(dataset4)
MP_obj5=NetCDFFile(dataset5)
MP_obj6=NetCDFFile(dataset6)
MP_obj7=NetCDFFile(dataset7)
#MP_SD=NetCDFFile(dataset_sd)
#MP_FL=NetCDFFile(dataset_flag)

MPX=array(MP_obj1.variables['x'][:])
MPY=array(MP_obj1.variables['y'][:])
MP0=array(MP_obj0.variables['z'][:])
MP1=array(MP_obj1.variables['z'][:])
MP2=array(MP_obj2.variables['z'][:])
MP3=array(MP_obj3.variables['z'][:])
MP4=array(MP_obj4.variables['z'][:])
MP5=array(MP_obj5.variables['z'][:])
MP6=array(MP_obj6.variables['z'][:])
MP7=array(MP_obj7.variables['z'][:])
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
    mpval=MP7[iY[n],iX[n]]
    MPval.append(mpval)

MPval=array(MPval)*100
MPval[44:]

MPval=array([ 5.34935856,  34.6866951 ,  37.59857559,  35.63520432,
                29.64470673,  30.76460075,  26.07561111,  27.84197235, 29.52030373,  29.20651817,  30.03935242,  27.90864754,
                24.28595924,  23.60974312,  30.97301483,  18.07652473, 18.93989182,  27.41495705,  23.37692833,
                         23.54593658,  27.13404846,  21.43231201, 19.53557777,  28.02973747,  11.1829567 ,  12.17084122,
                        19.09560966,  20.04397774,
                        13.71156406,  16.42700958,  31.85304451,           37.20304108,  38.25751495,  
                28.8565731 ,  39.40811539], dtype=float32)

mpPS2007=array([ 0.3 ,  0.38,  0.27,  0.35,  
        0.25,  0.24,  0.3 ,  0.31,  0.28,  0.35,  0.33,  0.39,  
        0.31,  0.29,  0.28,  0.29,  0.33,  0.24,  0.25,  
          0.24,  0.26,  0.23,  0.26,  0.23,  0.28,  0.32,  
         0.33,  0.34,
         0.47,  0.4 ,  0.17,    0.29,  0.3 ,   
          0.17,  0.06])




MPval=MPval*0.01
cc=corrcoef([MPval,mpPS2007])
print cc
RMSE=sqrt((1./n)*sum((MPval-mpPS2007)**2))
print RMSE

figure()
plot(MPval,mpPS2007,'k+')
title('mp MODIS 12.5km grid vs Ship Observation PS 2007')
ylabel('MP from POLARSTERN [%]')
xlabel('MP from MODIS [%]')




#from scipy.stats import linregress
cc=corrcoef([MPval,mpPS2007])
print cc

gx=linspace(0.,.5,100)
m=1.
b=0.
gy=m*gx+b
plot(gx,gy,'gray')



annotate('RMSE='+str(round(RMSE,3)),xy=(0.1,0.030))
#annotate('corr='+str(round(cc[1][0],3)),xy=(0.1,0.030))
savefig('/pf/u/u241127/plots/validation_ps2007.png')










