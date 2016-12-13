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


#run validation_meltex_4and7_6_v2.py


#snow
afolder='/scratch/clisap/seaice/AERIAL/MELTEX/meltponds/'
daylyfile4=afolder+'meltex_20080607.dat'
daylyfile2=afolder+'meltex_20080604.dat'

###read meltex data 2008
a2=loadtxt(daylyfile2,skiprows=1)
lat2=a2[:,1]
lon2=a2[:,2]
mp2=a2[:,4]
ow2=a2[:,5]
ice2=a2[:,6]
x2,y2=mapll(lat2,lon2,1)

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
wfolder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/PRODUCTS/2008_153/'

##weekly modis
dataset_week=wfolder+'2008_153_rmp_125masked.nc'
sd_mod_w=wfolder+'2008_153_mp_125SD.nc'
MP_objw=NetCDFFile(dataset_week)
SD_objw=NetCDFFile(sd_mod_w)
MPw=array(MP_objw.variables['z'][:])
SD_MODw=array(SD_objw.variables['z'][:])

#modis 4 june 2008
dataset2=folder+'2008_156/2008_156_rmp_125masked.nc'
sd_mod_2=folder+'2008_156/2008_156_mp_125SD.nc'
MP_obj2=NetCDFFile(dataset2)
SD_obj2=NetCDFFile(sd_mod_2)
MPX=array(MP_obj2.variables['x'][:])
MPY=array(MP_obj2.variables['y'][:])

MP2=array(MP_obj2.variables['z'][:])
SD_MOD2=array(SD_obj2.variables['z'][:])

#modis 7 june 2008
dataset4=folder+'2008_159/2008_159_rmp_125masked.nc'
sd_mod_4=folder+'2008_159/2008_159_mp_125SD.nc'
date4=dataset4[-17:-9]
jahr=date4[0:4]
julday4=date4[5:8]
MP_obj4=NetCDFFile(dataset4)
SD_obj4=NetCDFFile(sd_mod_4)

MP4=array(MP_obj4.variables['z'][:])
SD_MOD4=array(SD_obj4.variables['z'][:])









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
latmax,lonmax=74.2,-150.

Xmax,Ymax=mapll(latmax,lonmax,1)
Xmin,Ymin=mapll(latmin,lonmin,1)

iXmax=nonzero(MPX>=Xmax)[0][0]
iXmin=nonzero(MPX>=Xmin)[0][0]
iYmax=nonzero(MPY>=Ymax)[0][0]
iYmin=nonzero(MPY>=Ymin)[0][0]

MPX=MPX[iXmin:iXmax]
MPY=MPY[iYmin:iYmax]
MP2=MP2[iYmin:iYmax,iXmin:iXmax]
SD_MOD2=SD_MOD2[iYmin:iYmax,iXmin:iXmax]
MP4=MP4[iYmin:iYmax,iXmin:iXmax]
SD_MOD4=SD_MOD4[iYmin:iYmax,iXmin:iXmax]
MPw=MPw[iYmin:iYmax,iXmin:iXmax]
SD_MODw=SD_MODw[iYmin:iYmax,iXmin:iXmax]


#Tmap=Basemap(resolution='h',projection='stere', boundinglat=70., lon_0=lon_0, lat_0=lat_0, width=plotsize, height=plotsize, rsphere=(6378273.,6356889.))
Tmap=Basemap(resolution='h',projection='stere', boundinglat=60.,llcrnrlon=lonmin, llcrnrlat=latmin, urcrnrlon=lonmax, urcrnrlat=latmax, lon_0=lon_0, lat_0=lat_0)


figure(14)
Tmap.drawcoastlines()
Tmap.drawmeridians(np.arange(0,360,10),labels=[1,0,0,0])
Tmap.drawparallels(np.arange(-90,90,5),labels=[0,0,0,1])
x4,y4=Tmap(lon4,lat4)
#Tmap.imshow(MP4)
Tmap.imshow(MP4, interpolation='nearest')
Tmap.scatter(x4,y4,s=15,c=mp4,cmap='jet', marker='o',vmin=0.,vmax=40.,edgecolors= 'None')
colorbar()
text(890000.0,760000.0,'melt pond fraction [%]', rotation='vertical',fontsize='x-large')
title('7 June 2008 on daily modis',fontsize='x-large')
#savefig('/pf/u/u241127/plots/meltex_7_6_2008.png')

figure(15)
Tmap.drawcoastlines()
Tmap.drawmeridians(np.arange(0,360,10),labels=[1,0,0,0])
Tmap.drawparallels(np.arange(-90,90,5),labels=[0,0,0,1])
x2,y2=Tmap(lon2,lat2)
#Tmap.imshow(MP4)
Tmap.imshow(MP2, interpolation='nearest')
Tmap.scatter(x2,y2,s=15,c=mp2,cmap='jet', marker='o',vmin=0.,vmax=40.,edgecolors= 'None')
colorbar()
text(890000.0,760000.0,'melt pond fraction [%]', rotation='vertical',fontsize='x-large')
title('4 June 2008 on daily modis',fontsize='x-large')
savefig('/pf/u/u241127/plots/meltex_4_6_2008.png')


figure(16) ####on weekly modis
Tmap.drawcoastlines()
Tmap.drawmeridians(np.arange(0,360,10),labels=[1,0,0,0])
Tmap.drawparallels(np.arange(-90,90,5),labels=[0,0,0,1])
x2,y2=Tmap(lon2,lat2)
Tmap.imshow(MPw, interpolation='nearest')
Tmap.scatter(x2,y2,s=15,c=mp2,cmap='jet', marker='o',vmin=0.,vmax=40.,edgecolors= 'None')
colorbar()
text(890000.0,760000.0,'melt pond fraction [%]', rotation='vertical',fontsize='x-large')
title('4 June 2008 on weekly modis',fontsize='x-large')
#savefig('/pf/u/u241127/plots/meltex_4_6_2008_weekly.png')

figure(17) ####on weekly modis
Tmap.drawcoastlines()
Tmap.drawmeridians(np.arange(0,360,10),labels=[1,0,0,0])
Tmap.drawparallels(np.arange(-90,90,5),labels=[0,0,0,1])
x2,y2=Tmap(lon2,lat2)
Tmap.imshow(MPw, interpolation='nearest')
Tmap.scatter(x4,y4,s=15,c=mp4,cmap='jet', marker='o',vmin=0.,vmax=40.,edgecolors= 'None')
colorbar()
text(890000.0,760000.0,'melt pond fraction [%]', rotation='vertical',fontsize='x-large')
title('7 June 2008 on weekly modis',fontsize='x-large')
savefig('/pf/u/u241127/plots/meltex_7_6_2008_weekly.png')


####gridden airplane koord auf 12.5 km gitter und plot with MODIS

####7.6.
X4_airplane,Y4_airplane=mapll(lat4,lon4,1)
G={}
cs1=12.5
x0,y0,x2,y2=X4_airplane.max(),Y4_airplane.max(),X4_airplane.min(),Y4_airplane.min()
x0,y0,x2,y2=-1893.75, 293.75, -2168.75,81.25 #per Hand deie naechstliegenden Koodinaten aus MPX /MPY rausgesucht
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
#cmd('gv '+mapfile125+'&')

cmd('gmtconvert '+outfile_table_125+' -F0,1,3 -V > '+SDfile125)
cmd('xyz2grd '+SDfile125+opts(G,['Rx','I'])+' -V -G'+SDoutfile_mp125)
cmd('grdimage '+SDoutfile_mp125+' '+opts(G,['Rx','Jx','Bx'])+' -C'+ctabfilestdev+' -P -K -V> '+SDmapfile125)
cmd('psscale -D8c/-0.5c/12c/0.4ch -C'+ctabfilestdev+' -B0.1 -V -K -O>> '+SDmapfile125)
cmd('pscoast '+opts(G,['Rll','Bll','coast'])+' -G128/128/128 -O -V  >> '+SDmapfile125)
#cmd('gv '+SDmapfile125+'&')

mp_airpl4=NetCDFFile(outfile_mp_125)
sd_airpl4=NetCDFFile(SDoutfile_mp125)

MPX_airpl4=array(mp_airpl4.variables['x'][:])
MPY_airpl4=array(mp_airpl4.variables['y'][:])
MP_airpl4=array(mp_airpl4.variables['z'][:,:])
SD_airpl4=array(sd_airpl4.variables['z'][:,:])

ix14=nonzero(MPX==MPX_airpl4[0])[0][0]-1
ix24=nonzero(MPX==MPX_airpl4[-1])[0][0]
iy14=nonzero(MPY==MPY_airpl4[0])[0][0]-1
iy24=nonzero(MPY==MPY_airpl4[-1])[0][0]


MP_MOD4=MP4[iy14:iy24,ix14:ix24]*100
SD_MOD4=SD_MOD4[iy14:iy24,ix14:ix24]*100


MPw=MPw[iy14:iy24,ix14:ix24]*100
SD_MODw=SD_MODw[iy14:iy24,ix14:ix24]*100

####plot
#####spatial plot with gridded airplane-data
#######new subset...........
###7.6.
latmax4,lonmax4=mapxy(MPX_airpl4.max(),MPY_airpl4.max(),1)
latmin4,lonmin4=mapxy(MPX_airpl4.min(),MPY_airpl4.min(),1)


Tmap=Basemap(resolution='h',projection='stere', boundinglat=60.,llcrnrlon=lonmin4, llcrnrlat=latmin4, urcrnrlon=lonmax4, urcrnrlat=latmax4, lon_0=lon_0, lat_0=lat_0)

MPX_a4,MPY_a4=meshgrid(MPX_airpl4,MPY_airpl4)
lat_air4,lon_air4=mapxy(MPX_a4,MPY_a4, 1)


figure(1)
Tmap.drawcoastlines()
Tmap.drawmeridians(np.arange(0,360,5),labels=[1,0,0,0])
Tmap.drawparallels(np.arange(-90,90,2),labels=[0,0,0,1])
x4a,y4a=Tmap(lon_air4,lat_air4)
#Tmap.imshow(MP4)
Tmap.imshow(MP_MOD4, interpolation='nearest')
Tmap.scatter(x4a,y4a,s=35,c=MP_airpl4,cmap='jet', marker='o',vmin=0.,vmax=40.,edgecolors= 'None')
colorbar()
text(830000.0,630000.0,'melt pond fraction [%]', rotation='vertical',fontsize='x-large')
title('7 June 2008',fontsize='x-large')
savefig('/pf/u/u241127/plots/meltex_7_6_2008_zoom.png')

figure(2)
Tmap.drawcoastlines()
Tmap.drawmeridians(np.arange(0,360,5),labels=[1,0,0,0])
Tmap.drawparallels(np.arange(-90,90,2),labels=[0,0,0,1])
x4a,y4a=Tmap(lon_air4,lat_air4)
#Tmap.imshow(MP4)
Tmap.imshow(MPw, interpolation='nearest')
Tmap.scatter(x4a,y4a,s=35,c=MP_airpl4,cmap='jet', marker='o',vmin=0.,vmax=40.,edgecolors= 'None')
colorbar()
text(830000.0,630000.0,'melt pond fraction [%]', rotation='vertical',fontsize='x-large')
title('7 June 2008 on weekly modis',fontsize='x-large')
savefig('/pf/u/u241127/plots/meltex_7_6_2008_zoom_weekly.png')








####4.6.
X2_airplane,Y2_airplane=mapll(lat2,lon2,1)
G={}
cs1=12.5
x0,y0,x2,y2=X2_airplane.max(),Y2_airplane.max(),X2_airplane.min(),Y2_airplane.min()
x0,y0,x2,y2=-2068.75,343.75, -2181.25,18.75 #per Hand die naechstliegenden Koodinaten aus MPX /MPY rausgesucht
grid_par_neu(1,x0,y0,x2,y2,cs1,G)

table_airplane=folder+'2008_156/grd_airplane.xyz'
table_airplaneSD=folder+'2008_156/grd_airplaneSD.xyz'
table_modis=folder+'2008_156/grd_modis.xyz'
ctabfilemp=folder+'mp.ctab'
ctabfilestdev=folder+'stdev.ctab'
write_table(table_airplane,X2_airplane,Y2_airplane,mp2)


outfile_table_125=folder+'2008_156/2008_156_mp_125.xyz'
outfile_mp_125=folder+'2008_156/2008_156_mp_125.nc'
mapfile125=plotfolder+'2008_156_mp_125meltex.ps'
SDfile125=folder+'2008_156/2008_156_mp_SD_125.xyz'
SDoutfile_mp125=folder+'2008_156/2008_156_mp_SD_125.nc'
SDmapfile125=plotfolder+'2008_156_mp_125meltexSD.ps'
outfile_mp_125modis=folder+'2008_156/2008_156_mp_125modis.nc'
outfile_table_125modis=folder+'2008_156/2008_156_mp_125modis.xyz'


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

mp_airpl2=NetCDFFile(outfile_mp_125)
sd_airpl2=NetCDFFile(SDoutfile_mp125)

MPX_airpl2=array(mp_airpl2.variables['x'][:])
MPY_airpl2=array(mp_airpl2.variables['y'][:])
MP_airpl2=array(mp_airpl2.variables['z'][:,:])
SD_airpl2=array(sd_airpl2.variables['z'][:,:])

ix12=nonzero(MPX==MPX_airpl2[0])[0][0]-1
ix22=nonzero(MPX==MPX_airpl2[-1])[0][0]
iy12=nonzero(MPY==MPY_airpl2[0])[0][0]-1
iy22=nonzero(MPY==MPY_airpl2[-1])[0][0]

MP_MOD2=MP2[iy12:iy22,ix12:ix22]*100
SD_MOD2=SD_MOD2[iy12:iy22,ix12:ix22]*100



MPw=MPw[iy12:iy22,ix12:ix22]*100
SD_MODw=SD_MODw[iy12:iy22,ix12:ix22]*100
#####spatial plot with gridded airplane-data
#######new subset...........
###7.6.
latmax4,lonmax4=mapxy(MPX_airpl4.max(),MPY_airpl4.max(),1)
latmin4,lonmin4=mapxy(MPX_airpl4.min(),MPY_airpl4.min(),1)


Tmap=Basemap(resolution='h',projection='stere', boundinglat=60.,llcrnrlon=lonmin4, llcrnrlat=latmin4, urcrnrlon=lonmax4, urcrnrlat=latmax4, lon_0=lon_0, lat_0=lat_0)

MPX_a4,MPY_a4=meshgrid(MPX_airpl4,MPY_airpl4)
lat_air4,lon_air4=mapxy(MPX_a4,MPY_a4, 1)


figure(1)
Tmap.drawcoastlines()
Tmap.drawmeridians(np.arange(0,360,5),labels=[1,0,0,0])
Tmap.drawparallels(np.arange(-90,90,2),labels=[0,0,0,1])
x4a,y4a=Tmap(lon_air4,lat_air4)
#Tmap.imshow(MP4)
Tmap.imshow(MP_MOD4, interpolation='nearest')
Tmap.scatter(x4a,y4a,s=35,c=MP_airpl4,cmap='jet', marker='o',vmin=0.,vmax=40.,edgecolors= 'None')
colorbar()
text(830000.0,630000.0,'melt pond fraction [%]', rotation='vertical',fontsize='x-large')
title('7 June 2008',fontsize='x-large')
savefig('/pf/u/u241127/plots/meltex_7_6_2008_zoom.png')

figure(2)
Tmap.drawcoastlines()
Tmap.drawmeridians(np.arange(0,360,5),labels=[1,0,0,0])
Tmap.drawparallels(np.arange(-90,90,2),labels=[0,0,0,1])
x4a,y4a=Tmap(lon_air,lat_air)
#Tmap.imshow(MP4)
Tmap.imshow(MPw, interpolation='nearest')
Tmap.scatter(x4a,y4a,s=35,c=MP_airpl,cmap='jet', marker='o',vmin=0.,vmax=40.,edgecolors= 'None')
colorbar()
text(830000.0,630000.0,'melt pond fraction [%]', rotation='vertical',fontsize='x-large')
title('7 June 2008 on weekly modis',fontsize='x-large')
savefig('/pf/u/u241127/plots/meltex_7_6_2008_zoom_weekly.png')

###########4.6.
latmax2,lonmax2=mapxy(MPX_airpl2.max(),MPY_airpl2.max(),1)
latmin2,lonmin2=mapxy(MPX_airpl2.min(),MPY_airpl2.min(),1)


Tmap=Basemap(resolution='h',projection='stere', boundinglat=60.,llcrnrlon=lonmin2, llcrnrlat=latmin2, urcrnrlon=lonmax2, urcrnrlat=latmax2, lon_0=lon_0, lat_0=lat_0)

MPX_a2,MPY_a2=meshgrid(MPX_airpl2,MPY_airpl2)
lat_air2,lon_air2=mapxy(MPX_a2,MPY_a2, 1)


figure(3)
Tmap.drawcoastlines()
Tmap.drawmeridians(np.arange(0,360,5),labels=[1,0,0,0])
Tmap.drawparallels(np.arange(-90,90,2),labels=[0,0,0,1])
x2,y2=Tmap(lon2,lat2)

Tmap.imshow(MP2, interpolation='nearest')
Tmap.scatter(x2,y2,s=15,c=mp2,cmap='jet', marker='o',vmin=0.,vmax=40.,edgecolors= 'None')
colorbar()
text(830000.0,630000.0,'melt pond fraction [%]', rotation='vertical',fontsize='x-large')
title('4 June 2008',fontsize='x-large')
savefig('/pf/u/u241127/plots/meltex_4_6_2008_zoom.png')






#cut modis pixels 
ix1=nonzero(MPX==MPX_airpl[0])[0][0]-1
ix2=nonzero(MPX==MPX_airpl[-1])[0][0]
iy1=nonzero(MPY==MPY_airpl[0])[0][0]-1
iy2=nonzero(MPY==MPY_airpl[-1])[0][0]

MP_MOD=MP4[iy1:iy2,ix1:ix2]*100
SD_MOD=SD_MOD4[iy1:iy2,ix1:ix2]*100
#figure()
#imshow(MP_airpl,origin='lower',interpolation='nearest')
#colorbar()

#figure()
#imshow(MP_MOD,origin='lower',interpolation='nearest')
#colorbar()


i=isfinite(MP_airpl)
MP_airpl_i=MP_airpl[i]
MP_MOD_i=MP_MOD[i]
SD_airpl_i=SD_airpl[i]
SD_MOD_i=SD_MOD[i]

x=linspace(0.,MP_airpl_i.shape[0],MP_airpl_i.shape[0])

figure()
#plot(x,MP_airpl_i,'.')
errorbar(x, MP_airpl_i, yerr=SD_airpl_i, fmt='o',label='MELTEX w. STD')
#plot(x,MP_i,'ro',label='MODIS')
errorbar(x, MP_MOD_i, yerr=SD_MOD_i, fmt='ro',label='MODIS w. STD')

legend(numpoints=1)
title('Profile MODIS vs. MELTEX 12.5 km grid (7. June 2008)' )
xlabel('n')
ylabel('mp fraction [%]')


STOP


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





