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


#run validation_hotrax2005.py


#snow
folder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/SHIP_OBS/'
daylyfile=folder+'HotraxPonds_DAILY.txt'

###read ship data 2005
a=loadtxt(daylyfile,skiprows=1,delimiter=',')
a[a==999]=nan
day=a[:,0]
lat=a[:,1]
lon=a[:,2]
mpHotrax=a[:,4]
mpstd=a[:,6]
icek=a[:,3]
icekstd=a[:,5]

x,y=mapll(lat,lon,1)


figure()
errorbar(day,mpHotrax,yerr=mpstd,fmt='ko')
grid()
axis([220,270,-0.05,0.5])
ylabel('melt pond fraction')
xlabel('day of the year')
savefig('/pf/u/u241127/plots/mpfraction_HOTRAX2005.png')

#read mp data
####---snow1:
folder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/PRODUCTS/'
plotfolder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/PRODUCTS/plots/'


dataset0=folder+'2005_217/2005_217_rmp_125masked.nc'
dataset1=folder+'2005_225/2005_225_rmp_125masked.nc'
dataset2=folder+'2005_233/2005_233_rmp_125masked.nc'
dataset3=folder+'2005_241/2005_241_rmp_125masked.nc'
dataset4=folder+'2005_249/2005_249_rmp_125masked.nc'
dataset5=folder+'2005_257/2005_257_rmp_125masked.nc'

dataset0s=folder+'2005_217/2005_217_mp_125SD.nc'
dataset1s=folder+'2005_225/2005_225_mp_125SD.nc'
dataset2s=folder+'2005_233/2005_233_mp_125SD.nc'
dataset3s=folder+'2005_241/2005_241_mp_125SD.nc'
dataset4s=folder+'2005_249/2005_249_mp_125SD.nc'
dataset5s=folder+'2005_257/2005_257_mp_125SD.nc'


date0=dataset0[-17:-9]
date1=dataset1[-17:-9]#e.g.'2008_169'
print date1
date2=dataset2[-17:-9]
date3=dataset3[-17:-9]
date4=dataset4[-17:-9]
date5=dataset5[-17:-9]




jahr=date1[0:4]
julday1=date1[5:8]
julday2=date1[5:8]
julday3=date1[5:8]
julday4=date1[5:8]
julday5=date1[5:8]


MP_obj0=NetCDFFile(dataset0)
MP_obj1=NetCDFFile(dataset1)
MP_obj2=NetCDFFile(dataset2)
MP_obj3=NetCDFFile(dataset3)
MP_obj4=NetCDFFile(dataset4)
MP_obj5=NetCDFFile(dataset5)


MP_obj0s=NetCDFFile(dataset0s)
MP_obj1s=NetCDFFile(dataset1s)
MP_obj2s=NetCDFFile(dataset2s)
MP_obj3s=NetCDFFile(dataset3s)
MP_obj4s=NetCDFFile(dataset4s)
MP_obj5s=NetCDFFile(dataset5s)
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

SD0=array(MP_obj0s.variables['z'][:])
SD1=array(MP_obj1s.variables['z'][:])
SD2=array(MP_obj2s.variables['z'][:])
SD3=array(MP_obj3s.variables['z'][:])
SD4=array(MP_obj4s.variables['z'][:])
SD5=array(MP_obj5s.variables['z'][:])

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
    mpval=SD5[iY[n],iX[n]]
    MPval.append(mpval)

MPval=array(MPval)*100
MPval[37:]

MPval=array([23.68135834,  24.09506989,  24.3747406 ,
         27.69227982,  
         25.10074997,  21.6918335 ,  20.0218811 ,  21.55732727,
        10.7405901 ,   7.52149725,  12.45158958,  12.88492203,   8.32215977], dtype=float32)*0.01

SDval=array([2.61308599,  4.07714462,  5.19341803,
         1.3627516 ,  
        11.3501482 ,   9.18190575,   6.56585789,   5.10458469,
        10.43531036,   10.95269299,   9.77672291,     7.37194252,   6.70928431], dtype=float32)*0.01

mpHotrax=array([   0.19,  0.24,  0.38,  
                    0.33,  
                  0.22,  0.25,  0.2 ,  0.08,  
                0.18,    0.15,  0.01,    0.11,  0.1 ], dtype=float32)

mpstd=array([  0.12,  0.1 ,  0.03,  0.15,  
          0.03,  0.08,  0.09,  0.06,  0.08,  0.  ,  
0.  ,  0.04,    0.03 ])



#MPval.shape


cc=corrcoef([MPval,mpHotrax])
print cc

cc2=cc**2
print cc2


RMSE=sqrt((1./n)*sum((MPval-mpHotrax)**2))
print RMSE



figure(figsize=(6,6))
plot(mpHotrax,MPval,'k+')
axis([0,0.5,0,0.5])
#title('mp MODIS 12.5km grid vs Ship Observation HOTRAX 2005')
xlabel('Melt pond fraction from HOTRAX')
ylabel('Melt pond fraction from MODIS')


gx=linspace(0.,.5,100)
m=1.
b=0.
gy=m*gx+b
plot(gx,gy,'gray')

errorbar(mpHotrax,MPval,xerr=mpstd,yerr=SDval,fmt='ko')
annotate('RMSE='+str(round(RMSE,3)),xy=(0.01,0.47))
#annotate('RMSE='+str(round(RMSE,3)),xy=(0.4,0.030))

#annotate('corr='+str(round(cc[1][0],3)),xy=(0.1,0.030))
savefig('/pf/u/u241127/plots/validation_HOTRAX2005.png')










