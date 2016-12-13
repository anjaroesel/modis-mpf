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

# run validation_nsidc2000_2001.py

#local
#folder='/scratch/local3/Ship_Observ/'
#snow
folder='/scratch/clisap/seaice/AERIAL/NSIDC_meltpond/G02159_Cell_coverages_txt_files/'
beau_file=folder+'Beau2000.txt'
cana_file=folder+'Cana2000.txt'
fram_file=folder+'Fram2000.txt'
sibi_file=folder+'Esiber2000.txt'
#2001
beau_file2001=folder+'Beau2001.txt'
cana_file2001=folder+'Cana2001.txt'
fram_file2001=folder+'Fram2001.txt'
sibi_file2001=folder+'Esiber2001.txt'

###read spydata
a=loadtxt(beau_file)
day=a[:,0]
mp_beauf=a[:,1]
mpstd_b=a[:,2]
icekb=a[:,3]
icekstdb=a[:,4]

aa=loadtxt(beau_file2001)  
day=aa[:,0]
latb=73.
lonb=-150.
mp_beauf2001=aa[:,1]
mpstd2001_b=aa[:,2]
#icek2001=aa[:,3]
#icekstd2001=aa[:,4]
xb,yb=mapll(latb,lonb,1)

b=loadtxt(cana_file)
day=b[:,0]
latc=85.
lonc=-120.
mp_cana=b[:,1]
mpstd_c=b[:,2]
icekc=b[:,3]
icekstdc=b[:,4]
xc,yc=mapll(latc,lonc,1)

bb=loadtxt(cana_file2001)
day2001cana=bb[:,0]
mp_cana2001=bb[:,1]
mpstd_c2001=bb[:,2]
icekc2001=bb[:,3]
icekstdc2001=bb[:,4]


c=loadtxt(fram_file)
day=c[:,0]
latf=85.
lonf=0.
mp_fram=c[:,1]
mpstd_f=c[:,2]
icek=c[:,3]
icekstd=c[:,4]
xf,yf=mapll(latf,lonf,1)

cc=loadtxt(fram_file2001)
dayf2001=cc[:,0]
mp_fram2001=cc[:,1]
mpstd_f2001=cc[:,2]
icekf2001=cc[:,3]
icekstdf2001=cc[:,4]
xf,yf=mapll(latf,lonf,1)


d=loadtxt(sibi_file)
days=d[:,0]
lats=82.
lons=150.
mp_sibi=d[:,1]
mpstd_s=d[:,2]
iceks=d[:,3]
icekstds=d[:,4]
xs,ys=mapll(lats,lons,1)

dd=loadtxt(sibi_file2001)
days2001=dd[:,0]
mp_sibi2001=d[:,1]
mpstd_s2001=d[:,2]
iceks2001=d[:,3]
icekstds2001=d[:,4]

#read mp data
####---snow1:
folder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/PRODUCTS/'
plotfolder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/PRODUCTS/plots/'

#2000
dataset0=folder+'2000_169/2000_169_rmp_125masked.nc'   #2989809893047
dataset1=folder+'2000_177/2000_177_rmp_125masked.nc'
dataset2=folder+'2000_185/2000_185_rmp_125.nc'
dataset3=folder+'2000_193/2000_193_rmp_125.nc'
dataset4=folder+'2000_201/2000_201_rmp_125masked.nc'
dataset5=folder+'2000_209/2000_209_rmp_125masked.nc'
dataset6=folder+'2000_217/2000_217_rmp_125masked.nc'
dataset7=folder+'2000_225/2000_225_rmp_125masked.nc'
dataset8=folder+'2000_233/2000_233_rmp_125masked.nc'
dataset9=folder+'2000_241/2000_241_rmp_125masked.nc'
dataset10=folder+'2000_249/2000_249_rmp_125masked.nc'

dataset0s=folder+'2000_169/2000_169_mp_125SD.nc'
dataset1s=folder+'2000_177/2000_177_mp_125SD.nc'
dataset2s=folder+'2000_185/2000_185_mp_125SD.nc'
dataset3s=folder+'2000_193/2000_193_mp_125SD.nc'
dataset4s=folder+'2000_201/2000_201_mp_125SD.nc'
dataset5s=folder+'2000_209/2000_209_mp_125SD.nc'
dataset6s=folder+'2000_217/2000_217_mp_125SD.nc'
dataset7s=folder+'2000_225/2000_225_mp_125SD.nc'
dataset8s=folder+'2000_233/2000_233_mp_125SD.nc'
dataset9s=folder+'2000_241/2000_241_mp_125SD.nc'
dataset10s=folder+'2000_249/2000_249_mp_125SD.nc'

#2001
dataset0=folder+'2001_153/2001_153_rmp_125masked.nc'  #nan,8,4,1,2,nan,5
dataset1=folder+'2001_161/2001_161_rmp_125masked.nc'
dataset2=folder+'2001_185/2001_185_rmp_125masked.nc'
dataset3=folder+'2001_193/2001_193_rmp_125masked.nc'
dataset4=folder+'2001_201/2001_201_rmp_125masked.nc'
dataset5=folder+'2001_209/2001_209_rmp_125masked.nc'
dataset7=folder+'2001_225/2001_225_rmp_125masked.nc'
dataset8=folder+'2001_233/2001_233_rmp_125masked.nc'
dataset9=folder+'2001_241/2001_241_rmp_125masked.nc'
dataset10=folder+'2001_249/2001_249_rmp_125masked.nc'
dataset11=folder+'2001_217/2001_217_rmp_125masked.nc'

dataset0s=folder+'2001_153/2001_153_mp_125SD.nc'
dataset1s=folder+'2001_161/2001_161_mp_125SD.nc'
dataset2s=folder+'2001_185/2001_185_mp_125SD.nc'
dataset3s=folder+'2001_193/2001_193_mp_125SD.nc'
dataset4s=folder+'2001_201/2001_201_mp_125SD.nc'
dataset5s=folder+'2001_209/2001_209_mp_125SD.nc'
dataset7s=folder+'2001_225/2001_225_mp_125SD.nc'
dataset8s=folder+'2001_233/2001_233_mp_125SD.nc'
dataset9s=folder+'2001_241/2001_241_mp_125SD.nc'
dataset10s=folder+'2001_249/2001_249_mp_125SD.nc'
dataset11s=folder+'2001_217/2001_217_mp_125SD.nc'

date2=dataset2[-18:-10]
date3=dataset3[-18:-10]


jahr=date2[0:4]
julday2=date2[5:8]
julday3=date2[5:8]

MP_obj0=NetCDFFile(dataset0)
MP_obj1=NetCDFFile(dataset1)
MP_obj2=NetCDFFile(dataset2)
MP_obj3=NetCDFFile(dataset3)
MP_obj4=NetCDFFile(dataset4)
MP_obj5=NetCDFFile(dataset5)
#MP_obj6=NetCDFFile(dataset6)
MP_obj7=NetCDFFile(dataset7)
MP_obj8=NetCDFFile(dataset8)
MP_obj9=NetCDFFile(dataset9)
MP_obj10=NetCDFFile(dataset10)
MP_obj11=NetCDFFile(dataset11)


MP_obj0s=NetCDFFile(dataset0s)
MP_obj1s=NetCDFFile(dataset1s)
MP_obj2s=NetCDFFile(dataset2s)
MP_obj3s=NetCDFFile(dataset3s)
MP_obj4s=NetCDFFile(dataset4s)
MP_obj5s=NetCDFFile(dataset5s)
MP_obj6s=NetCDFFile(dataset6s)
MP_obj7s=NetCDFFile(dataset7s)
MP_obj8s=NetCDFFile(dataset8s)
MP_obj9s=NetCDFFile(dataset9s)
MP_obj10s=NetCDFFile(dataset10s)
MP_obj11s=NetCDFFile(dataset11s)


MPX=array(MP_obj2.variables['x'][:])
MPY=array(MP_obj2.variables['y'][:])
MP0=array(MP_obj0.variables['z'][:])
MP1=array(MP_obj1.variables['z'][:])
MP2=array(MP_obj2.variables['z'][:])
MP3=array(MP_obj3.variables['z'][:])
MP4=array(MP_obj4.variables['z'][:])
MP5=array(MP_obj5.variables['z'][:])
#MP6=array(MP_obj6.variables['z'][:])
MP7=array(MP_obj7.variables['z'][:])
MP8=array(MP_obj8.variables['z'][:])
MP9=array(MP_obj9.variables['z'][:])
MP10=array(MP_obj10.variables['z'][:])
MP11=array(MP_obj11.variables['z'][:])


SD0=array(MP_obj0s.variables['z'][:])
SD1=array(MP_obj1s.variables['z'][:])
SD2=array(MP_obj2s.variables['z'][:])
SD3=array(MP_obj3s.variables['z'][:])
SD4=array(MP_obj4s.variables['z'][:])
SD5=array(MP_obj5s.variables['z'][:])
SD6=array(MP_obj6s.variables['z'][:])
SD7=array(MP_obj7s.variables['z'][:])
SD8=array(MP_obj8s.variables['z'][:])
SD9=array(MP_obj9s.variables['z'][:])
SD10=array(MP_obj10s.variables['z'][:])
SD11=array(MP_obj11s.variables['z'][:])

"""
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

mpval=SD2[iY,iX] 
MPval.append(mpval)
"""
#Beauf data sets MP# :   #2989809893047
MPval_Beauf=array([0.20419148,0.1738494,0.15157218,0.1738494,0.15157218,nan,0.1738494,0.15157218,0.1738494,0.2188144,nan,0.33029205,nan], dtype=float32) 
mp_beauf=array([        nan,  0.05862319, nan,  0.06986486, nan,nan,0.08242754,0.09571429,0.11344595,0.2339, nan, 0.34048327,nan])
SDval_Beauf=array([0.0342136,0.12014791,0.15299113,0.12014791,0.15299113,nan,0.12014791 ,0.15299113,  0.12014791,0.03691563, nan,0.04642674,nan], dtype=float32)
mpstd_b=array([ nan,  0.0764788 , nan,  0.0947574 , nan, nan,  0.1060407 ,  0.09193438,  0.12925189,  0.04657773,nan,  0.09570705, nan])

#225711110
MPval_cana=array([nan,nan,nan,nan,nan,nan,nan,nan,0.13720128], dtype=float32)
mp_cana=array([ 0.16989796,  0.1841896 ,  0.20445283,  0.16690852,  0.17666667,0.17047431,  0.2957868 ,  0.36021739,  0.13488024])
SDval_cana=array([ nan,nan,nan,nan,nan,nan,nan,nan,0.035614], dtype=float32)
mpstd_c=array([ 0.05840659,  0.04227286,  0.03406465,  0.04504465,  0.0484768 , 0.05759544,  0.10254777,  0.13554107,  0.04656183])

#00680610
MPval_fram=array([0.13720128,0.13720128, 0.1430507,nan,0.13720128,0.1430507, nan,0.13720128], dtype=float32)
mp_fram=array([        nan,      nan,      nan,    nan,       nan,0.07093664, 0.18772222,nan])
SDval_fram=array([0.02921063,0.02921063, 0.,nan,0.02921063,0., nan,0.02921063], dtype=float32)
mpstd_f=array([        nan,         nan,         nan,         nan,         nan, 0.03924993,  0.05079129,         nan])

#359273
MPval_sibi=array([nan,nan,nan, nan,nan,nan], dtype=float32)
mp_sibi=array([ 0.20555556,  0.11229323,  0.13053846,  0.1594702 ,  0.10003876, 0.1162963 ])
SDval_sibi=array([nan,nan,nan, nan,nan,nan], dtype=float32)
mpstd_s=array([ 0.06652366,  0.04526239,  0.02912658,  0.09728849,  0.02673545, 0.08181349])

#11,4,1,3,1,2,9,1,5,8,9
MPval_Beauf2001=array([0.25575307,0.33029205,0.23484522 ,0.2188144,0.23484522, 0.20419148,0.1738494,0.23484522,nan,0.15157218,0.1738494 ], dtype=float32)
mp_beauf2001=array([nan,  0.17775216, nan,  nan,  0.25647541, 0.19464674,  0.13601824,nan,nan,  nan,   nan])
SDval_Beauf2001=array([0.01750151,0.04642674,0.09380564,0.03691563,0.09380564,0.0342136,0.12014791,0.09380564,0.03960273,0.15299113, 0.12014791], dtype=float32)
mpstd2001_b=array([        nan,  0.03927989,         nan,         nan,  0.05543575, 0.04919389,  0.09867543,  nan,    nan,  nan,  nan])


#2,8,2,11,2,11,10,5
MPval_cana2001=array([0.20419148,0.15157218,0.20419148,0.25575307,0.20419148,0.25575307,nan,nan], dtype=float32)
mp_cana2001=array([ 0.04538462,         nan,  0.07493827,  0.12819767,         nan,    0.02      ,         nan,  0.10492063])
SDval_cana2001=array([0.0342136, 0.15299113,0.0342136 ,0.01750151 ,0.0342136 ,0.01750151 ,nan,nan], dtype=float32)
mpstd_c2001=array([ 0.01748846,         nan,  0.03479701,  0.05884908,         nan,   0.01      ,         nan,  0.03647053])




# 9,4,nan,0,3,3,4,5,8,1,11,11,11,nan,1
MPval_fram2001=array([nan,nan,nan,nan,nan,nan,nan,nan,nan,nan, 0.20540446, 0.20540446, 0.20540446,nan,nan], dtype=float32)
mp_fram2001=array([nan,0.23059259,0.06582031,nan, 0.23875,nan,0.16307292,0.21179592,nan,nan,0.05640288,0.05727575,nan,0.04053381,nan])
SDval_fram2001=array([nan,nan,nan,nan,nan,nan,nan,nan,nan,nan, 0.04347075,0.04347075,0.04347075,nan,nan], dtype=float32)
mpstd_f2001=array([nan,0.14436253,0.03229932,nan,0.076844,nan,0.06192746,0.06456696,nan,nan,0.05782291,0.03972263,nan, 0.0365777,nan])


#nan,8,4,1,2,nan,5
MPval_sibi2001=array([nan,nan,nan,nan,nan,nan,nan], dtype=float32)
SDval_sibi2001=array([nan,nan,nan,nan,nan,nan,nan], dtype=float32)






mpval_MOD=array([0.1738494,0.1738494,0.1738494,0.15157218,0.1738494,0.2188144,0.33029205,0.13720128,0.1430507,0.33029205,0.23484522, 0.20419148,0.1738494,0.20419148,0.20419148,0.25575307,0.25575307,0.20540446, 0.20540446], dtype=float32)
mp_NSIDC=array([0.05862319,   0.06986486, 0.08242754,0.09571429,0.11344595,0.2339,  0.34048327,0.13488024,0.07093664, 0.17775216, 0.25647541, 0.19464674,  0.13601824,0.04538462, 0.07493827,  0.12819767,    0.02 , 0.05640288,0.05727575    ], dtype=float32)
sdval_MOD=array([0.12014791,0.12014791,0.12014791 ,0.15299113,  0.12014791,0.03691563,0.04642674,0.035614,0., 0.04642674,0.09380564,0.0342136,0.12014791,0.0342136, 0.0342136 ,0.01750151 ,0.01750151 , 0.04347075,0.04347075], dtype=float32)
sd_NSIDC=array([0.0764788 ,   0.0947574 ,  0.1060407 ,  0.09193438,  0.12925189,  0.04657773, 0.09570705,0.04656183,0.03924993, 0.03927989, 0.05543575, 0.04919389,  0.09867543, 0.01748846,  0.03479701,0.05884908,0.01 ,0.05782291,0.03972263], dtype=float32)



cc=corrcoef([mpval_MOD,mp_NSIDC])
print cc

cc2=cc**2
print cc2

n=mpval_MOD.shape[0]
RMSE=sqrt((1./n)*sum((mpval_MOD-mp_NSIDC)**2))

print 'RMSE: '+str(RMSE)



figure(figsize=(6,6))
lg=linregress(mpval_MOD,mp_NSIDC)
m=lg[0]
b=lg[1]
#gy=m*gx+b
gx=linspace(0.,0.5,100)
m=1.
b=0.
gy=m*gx+b
plot(gx,gy,'gray')
axis([0,0.5,0,0.5])
#plot(mp_NSIDC,mpval_MOD,'k+',label='All')
#errorbar(mp_NSIDC,mpval_MOD,xerr=sd_NSIDC,yerr=sdval_MOD,fmt='go',label='All')

errorbar(mp_beauf,MPval_Beauf,xerr=mpstd_b,yerr=SDval_Beauf,fmt='ko',label='Beaufort 2000')
errorbar(mp_cana,MPval_cana,xerr=mpstd_c,yerr=SDval_cana,fmt='o',color='aqua',label='Canada 2000')
errorbar(mp_fram,MPval_fram,xerr=mpstd_f,yerr=SDval_fram,fmt='bo',label='Fram 2000')
errorbar(mp_beauf2001,MPval_Beauf2001,xerr=mpstd2001_b,yerr=SDval_Beauf2001,fmt='kv',label='Beaufort 2001')
errorbar(mp_cana2001,MPval_cana2001,xerr=mpstd_c2001,yerr=SDval_cana2001,fmt='v',color='aqua',label='Canada 2001')
errorbar(mp_fram2001,MPval_fram2001,xerr=mpstd_f2001,yerr=SDval_fram2001,fmt='bv',label='Fram 2001')


legend(loc=4, numpoints=1)
#title('mp MODIS 12.5km grid vs NSIDC 2000-2001')
xlabel('Melt pond fraction from NSIDC')
ylabel('Melt pond fraction from MODIS')
#annotate('corr='+str(round(cc[1][0],3)),xy=(0.01,0.31))
annotate('RMSE='+str(round(RMSE,3)),xy=(0.01,0.47))
savefig('/pf/u/u241127/plots/validation_nsidc2000_2001_years.png')
#savefig('/pf/u/u241127/plots/validation_nsidc2000_2001_all.png')


