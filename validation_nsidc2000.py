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

#validation_nsidc2000.py

#local
#folder='/scratch/local3/Ship_Observ/'
#snow
folder='/scratch/clisap/seaice/AERIAL/NSIDC_meltpond/G02159_Cell_coverages_txt_files/'
beau_file=folder+'Beau2000.txt'
cana_file=folder+'Cana2000.txt'
fram_file=folder+'Fram2000.txt'
sibi_file=folder+'Esiber2000.txt'

###read spydata
a=loadtxt(beau_file)
day=a[:,0]
lat=73.
lon=-150.
mp_beauf=a[:,1]
mpstd=a[:,2]
icek=a[:,3]
icekstd=a[:,4]
x,y=mapll(lat,lon,1)

b=loadtxt(cana_file)
day=b[:,0]
lat=85.
lon=-120.
mp_cana=b[:,1]
mpstd=b[:,2]
icek=b[:,3]
icekstd=b[:,4]
x,y=mapll(lat,lon,1)

c=loadtxt(fram_file)
day=c[:,0]
lat=85.
lon=0.
mp_fram=c[:,1]
mpstd=c[:,2]
icek=c[:,3]
icekstd=c[:,4]
x,y=mapll(lat,lon,1)

d=loadtxt(sibi_file)
day=d[:,0]
lat=82.
lon=150.
mp_sibi=d[:,1]
mpstd=d[:,2]
icek=d[:,3]
icekstd=d[:,4]
x,y=mapll(lat,lon,1)

#read mp data
####---snow1:
folder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/PRODUCTS/'
plotfolder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/PRODUCTS/plots/'
dataset0=folder+'2000_169/2000_169_mp_125.nc'
dataset1=folder+'2000_177/2000_177_mp_125.nc'
dataset2=folder+'2000_185/2000_185_mp_125.nc'
dataset3=folder+'2000_193/2000_193_mp_125.nc'
dataset4=folder+'2000_201/2000_201_mp_125.nc'
dataset5=folder+'2000_209/2000_209_mp_125.nc'
dataset6=folder+'2000_217/2000_217_mp_125.nc'
dataset7=folder+'2000_225/2000_225_mp_125.nc'
dataset8=folder+'2000_233/2000_233_mp_125.nc'
dataset9=folder+'2000_241/2000_241_mp_125.nc'
dataset10=folder+'2000_249/2000_249_mp_125.nc'

date0=dataset0[-17:-9]
date1=dataset1[-17:-9]#e.g.'2008_169'
print date1
date2=dataset2[-17:-9]
date3=dataset3[-17:-9]



jahr=date1[0:4]
julday0=date0[5:8]
julday1=date1[5:8]
julday2=date1[5:8]
julday3=date1[5:8]

MP_obj0=NetCDFFile(dataset0)
MP_obj1=NetCDFFile(dataset1)
MP_obj2=NetCDFFile(dataset2)
MP_obj3=NetCDFFile(dataset3)
MP_obj4=NetCDFFile(dataset4)
MP_obj5=NetCDFFile(dataset5)
MP_obj6=NetCDFFile(dataset6)
MP_obj7=NetCDFFile(dataset7)
MP_obj8=NetCDFFile(dataset8)
MP_obj9=NetCDFFile(dataset9)
MP_obj10=NetCDFFile(dataset10)
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
MP8=array(MP_obj8.variables['z'][:])
MP9=array(MP_obj9.variables['z'][:])
MP10=array(MP_obj10.variables['z'][:])

#STDEV=array(MP_SD.variables['z'][:])
#FLAG=array(MP_FL.variables['z'][:])

###getting indices of MPset
iX=[]
iY=[]

iyy=nonzero(MPY>=y)[0]
iyy=iyy[0]
iY.append(iyy)


ixx=nonzero(MPX>=x)[0]
ixx=ixx[0]
iX.append(ixx)

#I=array((iY,iX))


MPval=[]

mpval=MP8[iY,iX]
MPval.append(mpval)

#MPval=array(MPval)*100

MPval_Beauf=array([nan,  0.05956559, nan,0.05956559,nan,nan,   0.05956559,  nan, 0.05956559, 0.32390824,nan,0.41820455, nan ], dtype=float32)

MPval_cana=array([nan,nan,nan,nan,nan,0.15311033,nan,nan,0.15311033], dtype=float32)

MPval_fram=array([nan,nan,nan,nan,nan,nan, 0.17900483,nan], dtype=float32)

MPval_sibi=array([nan,nan,nan, 0.17900483,nan,nan], dtype=float32)


mpval=array([0.05956559, 0.05956559,   0.05956559,   0.05956559, 0.32390824,0.41820455,0.15311033,0.15311033,0.17900483,0.17900483], dtype=float32)
mp=array([0.05862319,   0.06986486,     0.08242754,    0.11344595,  0.2339,   0.34048327,0.17047431,0.13488024, 0.07093664, 0.1594702 ], dtype=float32)
cc=corrcoef([mpval,mp])
print cc

figure()

lg=linregress(mpval,mp)
m=lg[0]
b=lg[1]
#gy=m*gx+b
gx=linspace(0.,0.5,100)

m=1.
b=0.
gy=m*gx+b
plot(gx,gy,'gray')

plot(MPval_Beauf,mp_beauf,'k+',label='Beaufort')
plot(MPval_cana,mp_cana,'r+',label='Canada')
plot(MPval_fram,mp_fram,'b+',label='Fram')
plot(MPval_sibi,mp_sibi,'g+',label='Siberia')

legend(loc=2)
title('mp MODIS 12.5km grid vs NSIDC 2000')
xlabel('MP from NSIDC [%]')
ylabel('MP from MODIS [%]')
annotate('corr='+str(round(cc[1][0],3)),xy=(0.01,0.35))
savefig('/pf/u/u241127/plots/validation_nsidc2000')










show()







###playgrpund

MODIStile='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/GRIDDATAMOD09/MOD09A1/2001_177_h13v01/2001_177_h13v01_1.nc'

MP_o=NetCDFFile(MODIStile)

MPX=array(MP_o.variables['x'][:])
MPY=array(MP_o.variables['y'][:])
MP=array(MP_o.variables['z'][:])




