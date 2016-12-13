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
from scipy.stats import linregress

#validation_polarstern2004.py


#snow
folder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/SHIP_OBS/'
daylyfile=folder+'POLsternARK202_2004_DAYLY.csv'

###read ship data 2004
a=loadtxt(daylyfile,skiprows=1)
a[a==999]=nan
day=a[:,0]
lat=a[:,1]
lon=a[:,2]
mpPS2004=a[:,3]
mpstd=a[:,4]
icek=a[:,5]
icekstd=a[:,6]

x,y=mapll(lat,lon,1)




#read mp data
####---snow1:
folder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/PRODUCTS/'
plotfolder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/PRODUCTS/plots/'

dataset0=folder+'2004_201/2004_201_rmp_125masked.nc'
dataset1=folder+'2004_209/2004_209_rmp_125masked.nc'
dataset2=folder+'2004_217/2004_217_rmp_125masked.nc'
dataset3=folder+'2004_225/2004_225_rmp_125masked.nc'

date0=dataset0[-17:-9]
date1=dataset1[-17:-9]#e.g.'2008_169'
print date1
date2=dataset2[-17:-9]
date3=dataset3[-17:-9]



jahr=date1[0:4]
julday1=date1[5:8]
julday2=date1[5:8]
julday3=date1[5:8]

MP_obj0=NetCDFFile(dataset0)
MP_obj1=NetCDFFile(dataset1)
MP_obj2=NetCDFFile(dataset2)
MP_obj3=NetCDFFile(dataset3)


MPX=array(MP_obj1.variables['x'][:])
MPY=array(MP_obj1.variables['y'][:])
MP0=array(MP_obj0.variables['z'][:])
MP1=array(MP_obj1.variables['z'][:])
MP2=array(MP_obj2.variables['z'][:])
MP3=array(MP_obj3.variables['z'][:])


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
    mpval=MP1[iY[n],iX[n]]
    MPval.append(mpval)

MPval=array(MPval)*100
MPval[3:11]


MPval=array([29.22357559,            35.45396805, 
                  33.10559464,  37.23913574,
                  33.10559464,
                17.77221107,  19.98047256,  18.72385788,  21.08803177,
                19.70838928,  16.03284645,  23.66284752,  18.77792931,
                23.22715569,  20.8874321 ,  21.70614243,  22.68392563,
                23.9473381 ,  23.9473381 ,  23.39551163], dtype=float32)


mpPS2004=array([ 26.67,    40.  ,  
          23.53,  25.63,  
          23.85,  
        32.07,  19.12,  23.82,  18.08,  
        22.  ,  20.31,  23.05,  16.76,  
        13.13,  13.  ,  17.86,  14.  ,  
        10.  ,  15.  ,  15.  ])


#MPval=MPval*0.01
cc=corrcoef([MPval,mpPS2004])
print cc
RMSE=sqrt((1./n)*sum((MPval-mpPS2004)**2))
print RMSE

figure()
plot(MPval,mpPS2004,'k+')
title('mp MODIS 12.5km grid vs Ship Observation PS 2004')
ylabel('MP from POLARSTERN [%]')
xlabel('MP from MODIS [%]')

lg=linregress(mpval,mp)
m=lg[0]
b=lg[1]
#gy=m*gx+b
gx=linspace(0.,40,100)

m=1.
b=0.
gy=m*gx+b
plot(gx,gy,'gray')



annotate('RMSE='+str(round(RMSE,3)),xy=(5,0.5))
#annotate('corr='+str(round(cc[1][0],3)),xy=(5,0.5))
savefig('/pf/u/u241127/plots/validation_ps2004.png')













"""

###playground

MODIStile='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/GRIDDATAMOD09/MOD09A1/2001_177_h13v01/2001_177_h13v01_1.nc'

MP_o=NetCDFFile(MODIStile)

MPX=array(MP_o.variables['x'][:])
MPY=array(MP_o.variables['y'][:])
MP=array(MP_o.variables['z'][:])

"""


