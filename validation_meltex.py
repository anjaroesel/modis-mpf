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


#run validation_meltex.py


#snow
folder='/scratch/clisap/seaice/AERIAL/MELTEX/meltponds/'
daylyfile0=folder+'meltex_20080526.dat'
daylyfile1=folder+'meltex_20080603.dat'
daylyfile2=folder+'meltex_20080604.dat'
daylyfile3=folder+'meltex_20080606.dat'
daylyfile4=folder+'meltex_20080607.dat'

###read meltex data 2008
a0=loadtxt(daylyfile0,skiprows=1)
lat0=a0[:,1]
lon0=a0[:,2]
mp0=a0[:,4]
ow0=a0[:,5]
ice0=a0[:,6]
x0,y0=mapll(lat0,lon0,1)


a1=loadtxt(daylyfile1,skiprows=1)
lat1=a1[:,1]
lon1=a1[:,2]
mp1=a1[:,4]
ow1=a1[:,5]
ice1=a1[:,6]
x1,y1=mapll(lat1,lon1,1)

a2=loadtxt(daylyfile2,skiprows=1)
lat2=a2[:,1]
lon2=a2[:,2]
mp2=a2[:,4]
ow2=a2[:,5]
ice2=a2[:,6]
x2,y2=mapll(lat2,lon2,1)

a3=loadtxt(daylyfile3,skiprows=1)
lat3=a3[:,1]
lon3=a3[:,2]
mp3=a3[:,4]
ow3=a3[:,5]
ice3=a3[:,6]
x3,y3=mapll(lat3,lon3,1)

a4=loadtxt(daylyfile4,skiprows=1)
lat4=a4[:,1]
lon4=a4[:,2]
mp4=a4[:,4]
ow4=a4[:,5]
ice4=a4[:,6]
x4,y4=mapll(lat4,lon4,1)


#read mp data
####---snow1:
folder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/PRODUCTS/daily/'
plotfolder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/PRODUCTS/daily/plots/'


dataset0=folder+'2008_147/2008_147_rmp_125masked.nc'
dataset1=folder+'2008_155/2008_155_rmp_125masked.nc'
dataset2=folder+'2008_156/2008_156_rmp_125masked.nc'
dataset3=folder+'2008_158/2008_158_rmp_125masked.nc'
dataset4=folder+'2008_159/2008_159_rmp_125masked.nc'


date0=dataset0[-17:-9]
date1=dataset1[-17:-9]#e.g.'2008_169'
print date1
date2=dataset2[-17:-9]
date3=dataset3[-17:-9]
date4=dataset4[-17:-9]





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


#MP_SD=NetCDFFile(dataset_sd)
#MP_FL=NetCDFFile(dataset_flag)

MPX=array(MP_obj1.variables['x'][:])
MPY=array(MP_obj1.variables['y'][:])
MP0=array(MP_obj0.variables['z'][:])
MP1=array(MP_obj1.variables['z'][:])
MP2=array(MP_obj2.variables['z'][:])
MP3=array(MP_obj3.variables['z'][:])
MP4=array(MP_obj4.variables['z'][:])


#STDEV=array(MP_SD.variables['z'][:])
#FLAG=array(MP_FL.variables['z'][:])

###getting indices of MPset
###dataset0
iX0=[]
iY0=[]
for n in range(len(y0)):
    print n
    iyy=nonzero(MPY>=y0[n])[0]
    iyy=iyy[0]
    iY0.append(iyy)

for n in range(len(x0)):
    print n
    ixx=nonzero(MPX>=x0[n])[0]
    ixx=ixx[0]
    iX0.append(ixx)

MPval0=[]
for n in range(len(x0)):
    print n
    mpval=MP0[iY0[n],iX0[n]]
    MPval0.append(mpval)

MPval0=array(MPval0)*100

mp00=[]
SD0=[]
for n in range(len(MPval0)):
    print n
    i=nonzero(MPval0==MPval0[n])
    mittel=mean(mp0[i])
    mp00.append(mittel)
    STDEV=std(mp0[i])
    SD0.append(STDEV)

mp00=array(mp00)
SD0=array(SD0)

cc0=corrcoef([MPval0,mp00])
print cc0
RMSE0=sqrt((1./n)*sum((MPval0-mp00)**2))
print RMSE0



###dataset1
iX1=[]
iY1=[]
for n in range(len(y1)):
    print n
    iyy=nonzero(MPY>=y1[n])[0]
    iyy=iyy[0]
    iY1.append(iyy)

for n in range(len(x1)):
    print n
    ixx=nonzero(MPX>=x1[n])[0]
    ixx=ixx[0]
    iX1.append(ixx)

MPval1=[]
for n in range(len(x1)):
    print n
    mpval=MP1[iY1[n],iX1[n]]
    MPval1.append(mpval)

MPval1=array(MPval1)*100

mp11=[]
SD1=[]
for n in range(len(MPval1)):
    print n
    i=nonzero(MPval1==MPval1[n])
    mittel=mean(mp1[i])
    mp11.append(mittel)
    STDEV=std(mp1[i])
    SD1.append(STDEV)

mp11=array(mp11)
SD1=array(SD1)

cc1=corrcoef([MPval1,mp11])
print cc1
RMSE1=sqrt((1./n)*sum((MPval1-mp11)**2))
print RMSE1


###dataset2
iX2=[]
iY2=[]
for n in range(len(y2)):
    print n
    iyy=nonzero(MPY>=y2[n])[0]
    iyy=iyy[0]
    iY2.append(iyy)

for n in range(len(x2)):
    print n
    ixx=nonzero(MPX>=x2[n])[0]
    ixx=ixx[0]
    iX2.append(ixx)

MPval2=[]
for n in range(len(x2)):
    print n
    mpval=MP2[iY2[n],iX2[n]]
    MPval2.append(mpval)

MPval2=array(MPval2)*100

mp22=[]
SD2=[]
for n in range(len(MPval2)):
    print n
    i=nonzero(MPval2==MPval2[n])
    mittel=mean(mp2[i])
    mp22.append(mittel)
    STDEV=std(mp2[i])
    SD2.append(STDEV)

mp22=array(mp22)
SD2=array(SD2)

cc2=corrcoef([MPval2,mp22])
print cc2
RMSE2=sqrt((1./n)*sum((MPval2-mp22)**2))
print RMSE2


###dataset3
iX3=[]
iY3=[]
for n in range(len(y3)):
    print n
    iyy=nonzero(MPY>=y3[n])[0]
    iyy=iyy[0]
    iY3.append(iyy)

for n in range(len(x3)):
    print n
    ixx=nonzero(MPX>=x3[n])[0]
    ixx=ixx[0]
    iX3.append(ixx)

MPval3=[]
for n in range(len(x3)):
    print n
    mpval=MP3[iY3[n],iX3[n]]
    MPval3.append(mpval)

MPval3=array(MPval3)*100

mp33=[]
SD3=[]
for n in range(len(MPval3)):
    print n
    i=nonzero(MPval3==MPval3[n])
    mittel=mean(mp3[i])
    mp33.append(mittel)
    STDEV=std(mp3[i])
    SD3.append(STDEV)

mp33=array(mp33)
SD33=array(SD3)

cc3=corrcoef([MPval3,mp33])
print cc3
RMSE3=sqrt((1./n)*sum((MPval3-mp33)**2))
print RMSE3

###dataset4
iX4=[]
iY4=[]
for n in range(len(y4)):
    print n
    iyy=nonzero(MPY>=y4[n])[0]
    iyy=iyy[0]
    iY4.append(iyy)

for n in range(len(x4)):
    print n
    ixx=nonzero(MPX>=x4[n])[0]
    ixx=ixx[0]
    iX4.append(ixx)

MPval4=[]
for n in range(len(x4)):
    print n
    mpval=MP4[iY4[n],iX4[n]]
    MPval4.append(mpval)

MPval4=array(MPval4)*100

mp44=[]
SD4=[]
for n in range(len(MPval4)):
    print n
    i=nonzero(MPval4==MPval4[n])
    mittel=mean(mp4[i])
    mp44.append(mittel)
    STDEV=std(mp4[i])
    SD4.append(STDEV)

mp44=array(mp44)
SD4=array(SD4)

cc4=corrcoef([MPval4,mp44])
print cc4
RMSE4=sqrt((1./n)*sum((MPval4-mp44)**2))
print RMSE4








figure(figsize=(6,6))
#plot(MPval0,mp00,'k+')
gx=linspace(0.,60,100)
m=1.
b=0.
gy=m*gx+b
plot(gx,gy,'gray')
errorbar(MPval1,mp11,yerr=SD1,xerr=None,fmt='mo',label='3. June (Drift Ice)')
errorbar(MPval2,mp22,yerr=SD2,xerr=None,fmt='go',label='4. June (Drift Ice)')
errorbar(MPval3,mp33,yerr=SD3,xerr=None,fmt='yo',label='6. June (Fast Ice)')
errorbar(MPval4,mp44,yerr=SD4,xerr=None,fmt='bo',label='7. June (Drift Ice)')
errorbar(MPval0,mp00,yerr=SD0,xerr=None,fmt='ko',label='26. May (Drift Ice)')
legend(numpoints=1)
title('mp MODIS 12.5km grid vs MELTEX2008')
ylim([0,80])
xlim([0,80])
ylabel('MP from MELTEX [%]')
xlabel('MP from MODIS [%]')

#annotate('RMSE='+str(round(RMSE,3)),xy=(0.1,0.030))

#annotate('corr='+str(round(cc[1][0],3)),xy=(0.1,0.030))
savefig('/pf/u/u241127/plots/validation_MELTEX2008.png')


###############spatial plotting
      
plotsize= 5000000 #in meter
lon_0=-45. #center coordinates of the plot
lat_0=90.0
#lon_0=-110.0 
#lat_0=75.0
llcrnrlon=-90.
llcrnrlat=55.
urcrnrlon=90
urcrnrlat=55.

####cut region MELTEX
latmin,lonmin=68,-122.
latmax,lonmax=74.2,-145.

Xmax,Ymax=mapll(latmax,lonmax,1)
Xmin,Ymin=mapll(latmin,lonmin,1)

iXmax=nonzero(MPX>=Xmax)[0][0]
iXmin=nonzero(MPX>=Xmin)[0][0]
iYmax=nonzero(MPY>=Ymax)[0][0]
iYmin=nonzero(MPY>=Ymin)[0][0]

MPX=MPX[iXmin:iXmax]
MPY=MPY[iYmin:iYmax]
MP0=MP0[iYmin:iYmax,iXmin:iXmax]
MP1=MP1[iYmin:iYmax,iXmin:iXmax]
MP2=MP2[iYmin:iYmax,iXmin:iXmax]
MP3=MP3[iYmin:iYmax,iXmin:iXmax]
MP4=MP4[iYmin:iYmax,iXmin:iXmax]


#Tmap=Basemap(resolution='h',projection='stere', boundinglat=70., lon_0=lon_0, lat_0=lat_0, width=plotsize, height=plotsize, rsphere=(6378273.,6356889.))
Tmap=Basemap(resolution='h',projection='stere', boundinglat=60.,llcrnrlon=lonmin, llcrnrlat=latmin, urcrnrlon=lonmax, urcrnrlat=latmax, lon_0=lon_0, lat_0=lat_0)

figure(10)
Tmap.drawcoastlines()
Tmap.drawmeridians(np.arange(0,360,10),labels=[1,0,0,0])
Tmap.drawparallels(np.arange(-90,90,5),labels=[0,0,0,1])
x0,y0=Tmap(lon0,lat0)
#Tmap.imshow(MP0)
Tmap.imshow(MP0, interpolation='nearest')
Tmap.scatter(x0,y0,s=15,c=mp00,cmap='jet', marker='o',vmin=0.,vmax=40.,edgecolors= 'None')
colorbar()
text(830000.0,630000.0,'melt pond fraction [%]', rotation='vertical',fontsize='x-large')
title('26 May 2008',fontsize='x-large')
savefig('/pf/u/u241127/plots/meltex_26_5_2008.png')


figure(11)
Tmap.drawcoastlines()
Tmap.drawmeridians(np.arange(0,360,10),labels=[1,0,0,0])
Tmap.drawparallels(np.arange(-90,90,5),labels=[0,0,0,1])
x1,y1=Tmap(lon1,lat1)
#Tmap.imshow(MP1)
Tmap.imshow(MP1, interpolation='nearest')
Tmap.scatter(x1,y1,s=15,c=mp11,cmap='jet', marker='o',vmin=0.,vmax=40.,edgecolors= 'None')
colorbar()
text(830000.0,630000.0,'melt pond fraction [%]', rotation='vertical',fontsize='x-large')
title('3 June 2008',fontsize='x-large')
savefig('/pf/u/u241127/plots/meltex_3_6_2008.png')

figure(12)
Tmap.drawcoastlines()
Tmap.drawmeridians(np.arange(0,360,10),labels=[1,0,0,0])
Tmap.drawparallels(np.arange(-90,90,5),labels=[0,0,0,1])
x2,y2=Tmap(lon2,lat2)
#Tmap.imshow(MP2)
Tmap.imshow(MP2, interpolation='nearest')
Tmap.scatter(x2,y2,s=15,c=mp22,cmap='jet', marker='o',vmin=0.,vmax=40.,edgecolors= 'None')
colorbar()
text(830000.0,630000.0,'melt pond fraction [%]', rotation='vertical',fontsize='x-large')
title('4 June 2008',fontsize='x-large')
savefig('/pf/u/u241127/plots/meltex_4_6_2008.png')

figure(13)
Tmap.drawcoastlines()
Tmap.drawmeridians(np.arange(0,360,10),labels=[1,0,0,0])
Tmap.drawparallels(np.arange(-90,90,5),labels=[0,0,0,1])
x3,y3=Tmap(lon3,lat3)
#Tmap.imshow(MP3)
Tmap.imshow(MP3, interpolation='nearest')
Tmap.scatter(x3,y3,s=15,c=mp33,cmap='jet', marker='o',vmin=0.,vmax=40.,edgecolors= 'None')
colorbar()
text(830000.0,630000.0,'melt pond fraction [%]', rotation='vertical',fontsize='x-large')
title('6 June 2008')
savefig('/pf/u/u241127/plots/meltex_6_6_2008.png')

figure(14)
Tmap.drawcoastlines()
Tmap.drawmeridians(np.arange(0,360,10),labels=[1,0,0,0])
Tmap.drawparallels(np.arange(-90,90,5),labels=[0,0,0,1])
x4,y4=Tmap(lon4,lat4)
#Tmap.imshow(MP4)
Tmap.imshow(MP4, interpolation='nearest')
Tmap.scatter(x4,y4,s=15,c=mp44,cmap='jet', marker='o',vmin=0.,vmax=40.,edgecolors= 'None')
colorbar()
text(830000.0,630000.0,'melt pond fraction [%]', rotation='vertical',fontsize='x-large')
title('7 June 2008',fontsize='x-large')
savefig('/pf/u/u241127/plots/meltex_7_6_2008.png')
show()







####gridden airplane koord auf 12.5 km gitter und plot with MODIS
G={}
cs1=12.5
x0,y0,x2,y2=X4_airplane.max(),Y4_airplane.max(),X4_airplane.min(),Y4_airplane.min()
x0,y0,x2,y2=-1893.75, 293.75, -2168.75,81.25
grid_par_neu(1,x0,y0,x2,y2,cs1,G)


table_airplane=folder+'2008_159/grd_airplane.xyz'
table_airplaneSD=folder+'2008_159/grd_airplaneSD.xyz'
table_modis=folder+'2008_159/grd_modis.xyz'

ctabfilemp=folder+'mp.ctab'
ctabfilestdev=folder+'stdev.ctab'

write_table(table_airplane,X4_airplane,Y4_airplane,mp4)



outfile_table_125=folder+'2008_159/2008_159_mp_125.xyz'
outfile_mp_125=folder+'2008_159/2008_159_mp_125.nc'
mapfile125=plotfolder+'2008_159_mp_125meltex.ps'
SDfile125=folder+'2008_159/2008_159_mp_SD_125.xyz'
SDoutfile_mp125=folder+'2008_159/2008_159_mp_SD_125.nc'
SDmapfile125=plotfolder+'2008_159_mp_125meltexSD.ps'
outfile_mp_125modis=folder+'2008_159/2008_159_mp_125modis.nc'
outfile_table_125modis=folder+'2008_159/2008_159_mp_125modis.xyz'



os.system('blockmean '+table_airplane+opts(G,['Rx','I'])+' -Eslh -V > '+outfile_table_125)
os.system('xyz2grd '+outfile_table_125+opts(G,['Rx','I'])+' -V -G'+outfile_mp_125) 
cmd('grdimage '+outfile_mp_125+' '+opts(G,['Rx','Jx','Bx'])+' -C'+ctabfilemp+' -P -K -V> '+mapfile125)
cmd('psscale -D8c/-0.5c/12c/0.4ch -C'+ctabfilemp+' -Ef -O -B0.1 -V -K >> '+mapfile125)
cmd('pscoast '+opts(G,['Rll','Bll','coast'])+' -G128/128/128 -O -V  >> '+mapfile125)   
cmd('gv '+mapfile125+'&')

cmd('gmtconvert '+outfile_table_125+' -F0,1,3 -V > '+SDfile125)
cmd('xyz2grd '+SDfile125+opts(G,['Rx','I'])+' -V -G'+SDoutfile_mp125)
cmd('grdimage '+SDoutfile_mp125+' '+opts(G,['Rx','Jx','Bx'])+' -C'+ctabfilestdev+' -P -K -V> '+SDmapfile125)
cmd('psscale -D8c/-0.5c/12c/0.4ch -C'+ctabfilestdev+' -B0.1 -V -K -O>> '+SDmapfile125)
cmd('pscoast '+opts(G,['Rll','Bll','coast'])+' -G128/128/128 -O -V  >> '+SDmapfile125)
cmd('gv '+SDmapfile125+'&')




mp_airpl=NetCDFFile(outfile_mp_125)
sd_airpl=NetCDFFile(SDoutfile_mp125)

MPX_airpl=array(mp_airpl.variables['x'][:])
MPY_airpl=array(mp_airpl.variables['y'][:])
MP_airpl=array(mp_airpl.variables['z'][:,:])
SD_airpl=array(sd_airpl.variables['z'][:,:])

MP=MP4[49:,14:37]

#figure()
#imshow(MP_airpl,origin='lower',interpolation='nearest')
#colorbar()

#figure()
#imshow(MP,origin='lower',interpolation='nearest')
#colorbar()


i=isfinite(MP_airpl)
MP_airpl_i=MP_airpl[i]
MP_i=MP[i]*100
SD_airpl_i=SD_airpl[i]

x=linspace(0.,MP_airpl_i.shape[0],MP_airpl_i.shape[0])

figure()
#plot(x,MP_airpl_i,'.')
errorbar(x, MP_airpl_i, yerr=SD_airpl_i, fmt='o',label='MELTEX w. STD')
plot(x,MP_i,'ro',label='MODIS')
legend(numpoints=1)
title('Profile MODIS vs. MELTEX 12.5 km grid (7. June 2008)' )
xlabel='n'
ylabel='mp fraction [%]'





#####scale to 50km pixel
#folder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/PRODUCTS/daily/'
#plotfolder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/PRODUCTS/daily/plots/'
#dataset0=folder+'2008_147/2008_147_rmp_125masked.nc'
#dataset1=folder+'2008_155/2008_155_rmp_125masked.nc'
#dataset2=folder+'2008_156/2008_156_rmp_125masked.nc'
#dataset3=folder+'2008_158/2008_158_rmp_125masked.nc'
#dataset4=folder+'2008_159/2008_159_rmp_125masked.nc'
#7 June 2008
i=isfinite(MP4)
mpi=MP4[i]

def XYgrid(x0,y0,x2,y2,cs):
    y=linspace(y2+cs/2,y0-cs/2,(y0-y2)/cs)
    x=linspace(x2+cs/2,x0-cs/2,(x0-x2)/cs)
    return meshgrid(x,y) 

x44,y44=XYgrid(Xmax,Ymax,Xmin,Ymin,12.5)
x4i=x44[i]
y4i=y44[i]

FILENAME=folder+'2008_159/2008_159_mp.xyz'
ctabfilemp=folder+'mp.ctab'
ctabfilestdev=folder+'stdev.ctab'
write_table(FILENAME,x4i,y4i,mpi)
outfile_table_50=folder+'2008_159/2008_159_mp_50.xyz'
outfile_mp_50=folder+'2008_159/2008_159_mp_50.nc'
mapfile50=plotfolder+'2008_159_mp_50meltex.ps'
SDfile50=folder+'2008_159/2008_159_mp_SD_50.xyz'
SDoutfile_mp50=folder+'2008_159/2008_159_mp_SD_50.nc'
SDmapfile=plotfolder+'2008_159_mp_50meltexSD.ps'
G={}
cs=50.
grid_par_neu(1,Xmax,Ymax,Xmin,Ymin,cs,G)

cmd('gmtset OBLIQUE_ANOTATION 1 PAPER_MEDIA A4+ PAGE_ORIENTATION portrait PLOT_DEGREE_FORMAT ddd:mm:ss -V')
cmd('gmtset COLOR_NAN 200')
os.system('blockmean '+FILENAME+opts(G,['Rx','I'])+' -Eslh -V > '+outfile_table_50)
os.system('xyz2grd '+outfile_table_50+opts(G,['Rx','I'])+' -V -G'+outfile_mp_50) 
cmd('grdimage '+outfile_mp_50+' '+opts(G,['Rx','Jx','Bx'])+' -C'+ctabfilemp+' -P -K -V> '+mapfile50)
cmd('psscale -D8c/-0.5c/12c/0.4ch -C'+ctabfilemp+' -Ef -O -B0.1 -V -K >> '+mapfile50)
cmd('pscoast '+opts(G,['Rll','Bll','coast'])+' -G128/128/128 -O -V  >> '+mapfile50)   
cmd('gv '+mapfile50+'&')

#stdev


cmd('gmtconvert '+outfile_table_50+' -F0,1,3 -V > '+SDfile50)
cmd('xyz2grd '+SDfile50+opts(G,['Rx','I'])+' -V -G'+SDoutfile_mp50)
cmd('grdimage '+SDoutfile_mp50+' '+opts(G,['Rx','Jx','Bx'])+' -C'+ctabfilestdev+' -P -K -V> '+SDmapfile)
cmd('psscale -D8c/-0.5c/12c/0.4ch -C'+ctabfilestdev+' -B0.1 -V -K -O>> '+SDmapfile)
cmd('pscoast '+opts(G,['Rll','Bll','coast'])+' -G128/128/128 -O -V  >> '+SDmapfile)
cmd('gv '+SDmapfile+'&')

x44,y44=XYgrid(Xmax,Ymax,Xmin,Ymin,cs)
x4i=x44[i]
y4i=y44[i]
X4_airplane,Y4_airplane=mapll(lat4,lon4,1)





