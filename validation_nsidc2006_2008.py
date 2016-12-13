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

# run validation_nsidc2006_2008.py

#local
#folder='/scratch/local3/Ship_Observ/'
#snow
folder='/scratch/clisap/seaice/AERIAL/NSIDC_meltpond/NSIDC_newMP_text_files/'
sibi_file6=folder+'Esiber2006.txt'
sibi_file7=folder+'Esiber2007.txt'
sibi_file8=folder+'Esiber2008.txt'

###read spydata
a=loadtxt(sibi_file6, dtype='S')
day6=array(a[:,1], dtype='float32')
lat=82.
lon=150.
mp_sib6=array(a[:,2], dtype='float32')
mpstd6=array(a[:,3], dtype='float32')
icek6=array(a[:,4], dtype='float32')
icekstd6=array(a[:,5], dtype='float32')
x,y=mapll(lat,lon,1)


b=loadtxt(sibi_file7, dtype='S')
day7=array(b[:,1], dtype='float32')
lat=82.
lon=150.
mp_sib7=array(b[:,2], dtype='float32')
mpstd7=array(b[:,3], dtype='float32')
icek7=array(b[:,4], dtype='float32')
icekstd7=array(b[:,5]
#x,y=mapll(lat,lon,1)

c=loadtxt(sibi_file8, dtype='S')
day8=array(c[:,1], dtype='float32')
lat=82.
lon=150.
mp_sib8=array(c[:,2], dtype='float32')
mpstd8=array(c[:,3], dtype='float32')
icek8=array(c[:,4], dtype='float32')
icekstd8=array(c[:,5], dtype='float32')
#x,y=mapll(lat,lon,1)



#read mp data
####---snow1:
folder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/PRODUCTS/'
plotfolder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/PRODUCTS/plots/'


dataset0=folder+'2006_249/2006_249_rmp_125masked.nc'
#dataset1=folder+'2006_257/2006_257_rmp_125.nc'
#dataset2=folder+'2007_249/2007_249_rmp_125.nc'
#dataset3=folder+'2007_257/2007_257_rmp_125.nc'
dataset4=folder+'2008_193/2008_193_rmp_125masked.nc'
dataset5=folder+'2008_209/2008_209_rmp_125masked.nc'
dataset6=folder+'2008_249/2008_249_rmp_125masked.nc'

dataset0s=folder+'2006_249/2006_249_mp_125SD.nc'
dataset4s=folder+'2008_193/2008_193_mp_125SD.nc'
dataset5s=folder+'2008_209/2008_209_mp_125SD.nc'
dataset6s=folder+'2008_249/2008_249_mp_125SD.nc'



date0=dataset0[-18:-10]
date1=dataset1[-18:-10]
date2=dataset2[-18:-10]
date3=dataset3[-18:-10]
date4=dataset4[-18:-10]
date5=dataset5[-18:-10]
date6=dataset6[-18:-10]

jahr6=date0[0:4]
jahr7=date2[0:4]
jahr8=date5[0:4]

julday0=date0[5:8]
julday1=date1[5:8]
julday2=date2[5:8]
julday3=date3[5:8]
julday4=date4[5:8]
julday5=date5[5:8]
julday6=date6[5:8]


MP_obj0=NetCDFFile(dataset0)
#MP_obj1=NetCDFFile(dataset1)
#MP_obj2=NetCDFFile(dataset2)
#MP_obj3=NetCDFFile(dataset3)
MP_obj4=NetCDFFile(dataset4)
MP_obj5=NetCDFFile(dataset5)
MP_obj6=NetCDFFile(dataset6)

MP_obj0s=NetCDFFile(dataset0s)
MP_obj4s=NetCDFFile(dataset4s)
MP_obj5s=NetCDFFile(dataset5s)
MP_obj6s=NetCDFFile(dataset6s)



MPX=array(MP_obj2.variables['x'][:])
MPY=array(MP_obj2.variables['y'][:])
MP0=array(MP_obj0.variables['z'][:])
#MP1=array(MP_obj1.variables['z'][:])
#MP2=array(MP_obj2.variables['z'][:])
#MP3=array(MP_obj3.variables['z'][:])
MP4=array(MP_obj4.variables['z'][:])
MP5=array(MP_obj5.variables['z'][:])
MP6=array(MP_obj6.variables['z'][:])


SD0=array(MP_obj0s.variables['z'][:])
SD4=array(MP_obj4s.variables['z'][:])
SD5=array(MP_obj5s.variables['z'][:])
SD6=array(MP_obj6s.variables['z'][:])


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

mpval=SD6[iY,iX]
MPval.append(mpval)


MPval0=array([ 0.14584674], dtype=float32)
MPval1=array([nan], dtype=float32) #####non masked
MPval2=array([nan], dtype=float32)
MPval3=array([ nan], dtype=float32)
MPval4=array([ 0.33032274], dtype=float32)
MPval5=array([ 0.21965837], dtype=float32)
MPval6=array([ 0.25250077], dtype=float32)



#MPval=array(MPval)*100

sdval=array([ 0.05675108, 0.05024811, 0.05024811,0.03449684,0.03449684,0.09089205], dtype=float32)
mpval=array([0.14584674,0.33032274,0.33032274, 0.21965837,0.21965837,0.25250077], dtype=float32)
mp=array([0.23933518,  0.09333333,0.03557692,  0.04257143,  0.0528, 0.04177257], dtype=float32)
sd=array([0.03255937, 0.06027714,  0.01964523,  0.01956393,  0.01860107,  0.05224364], dtype=float32)

cc=corrcoef([mpval,mp])
print cc

n=mpval.shape[0]
RMSE=sqrt((1./n)*sum((mpval-mp)**2))

print 'RMSE: '+str(RMSE)



figure(figsize=(6,6))
lg=linregress(mpval,mp)
m=lg[0]
b=lg[1]
#gy=m*gx+b
gx=linspace(0.,0.5,100)
m=1.
b=0.
gy=m*gx+b
plot(gx,gy,'gray')
#plot(mp_sib6[0],mpval[0],'k+',label='Siberia 2006')
#plot(mp_sib8,mpval[1:],'kv',label='Siberia 2008')

errorbar(mp_sib6[0],mpval[0],xerr=sd[0],yerr=sdval[0],fmt='go',label='Siberia 2006')
errorbar(mp_sib8,mpval[1:],xerr=sd[1:],yerr=sdval[1:],fmt='bo',label='Siberia 2008')

axis([0,0.5,0,0.5])
legend(loc=1, numpoints=1)
title('mp MODIS 12.5km grid vs NSIDC 2006-2008')
xlabel('MP fraction from NSIDC')
ylabel('MP fraction from MODIS')
#annotate('corr='+str(round(cc[1][0],3)),xy=(0.01,0.31))
annotate('RMSE='+str(round(RMSE,3)),xy=(0.01,0.47))
savefig('/pf/u/u241127/plots/validation_nsidc2006_2008.png')






