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
import scipy.ndimage as ndi
from pyhdf import SD

#run validation_meltex_4_7_6_v7_ice_c.py

#snow
afolder='/scratch/clisap/seaice/AERIAL/MELTEX/meltponds/'
daylyfile7=afolder+'meltex_20080607.dat'
daylyfile4=afolder+'meltex_20080604.dat'


###read meltex data 2008
a7=loadtxt(daylyfile7,skiprows=1)
lat7=a7[:,1]
lon7=a7[:,2]
mp7=a7[:,4]
ow7=a7[:,5]
ice7=a7[:,6]
x7,y7=mapll(lat7,lon7,1)

###read meltex data 2008
a4=loadtxt(daylyfile4,skiprows=1)
lat4=a4[:,1]
lon4=a4[:,2]
mp4=a4[:,4]
ow4=a4[:,5]
ice4=a4[:,6]
x4,y4=mapll(lat4,lon4,1)

###gauss filter over the aircraft data
mp7_filtered=ndi.gaussian_filter1d(mp7,10)
mp7_filtered10=mp7_filtered[range(0,mp7_filtered.shape[0],10)] #nur jeden 20. wert
lat7_10=lat7[range(0,mp7_filtered.shape[0],10)]
lon7_10=lon7[range(0,mp7_filtered.shape[0],10)]
Xair710,Yair710=mapll(lat7_10,lon7_10,1)

mp7_10=mp7[range(0,mp7.shape[0],10)]
diff7=mp7_filtered10-mp7_10

ice7_filtered=ndi.gaussian_filter1d(ice7,10)
ice7_filtered10=ice7_filtered[range(0,ice7_filtered.shape[0],10)] #nur jeden 20. wert
ice7_10=ice7[range(0,ice7.shape[0],10)]


#x7=linspace(0.,mp7_filtered10.shape[0],mp7_filtered10.shape[0])
#figure()
#errorbar(x7,mp7_filtered10,yerr=diff7,fmt='o',label='MELTEX w. STD')

###gauss filter over the aircraft data
mp4_filtered=ndi.gaussian_filter1d(mp4,10)
mp4_10=mp4[range(0,mp4.shape[0],10)]
mp4_filtered10=mp4_filtered[range(0,mp4_filtered.shape[0],10)] #nur jeden 10. wert
lat4_10=lat4[range(0,mp4_filtered.shape[0],10)]
lon4_10=lon4[range(0,mp4_filtered.shape[0],10)]
Xair410,Yair410=mapll(lat4_10,lon4_10,1)

diff4=mp4_filtered10-mp4_10
#x4=linspace(0.,mp4_filtered10.shape[0],mp4_filtered10.shape[0])
#figure()
#errorbar(x4,mp4_filtered10,yerr=diff,fmt='o',label='MELTEX w. STD')
ice4_filtered=ndi.gaussian_filter1d(ice4,10)
ice4_filtered10=ice4_filtered[range(0,ice4_filtered.shape[0],10)] #nur jeden 20. wert
ice4_10=ice4[range(0,ice7.shape[0],10)]

###################read aircraft height
gpsfile7=afolder+'20080607_gps.asc'
gpsfile4=afolder+'20080604_gps.asc'

###read meltex data 2008
g7=loadtxt(gpsfile7,skiprows=1)
glat7=g7[:,3]
glat7=glat7[7000:]
glon7=g7[:,4]
h7=g7[:,5]

###read meltex data 2008
g4=loadtxt(gpsfile4,skiprows=1)
glat4=g4[:,3]
glon4=g4[:,4]
h4=g4[:,5]
gx4,gy4=mapll(glat4,glon4,1)



glat7_track=[]
for n in range(glat7.shape[0]):
    if glat7[n]<=lat7[0] and glat7[n]>=lat7[-1]:
        glat7_track.append(glat7[n])

glat7_track=array(glat7_track)
istart=nonzero(glat7==glat7_track[0])[0][0]
iend=nonzero(glat7==glat7_track[-1])[0][0]
h7=h7[istart:iend]



nonzero(lat4[0]<=glat4)
istart=nonzero(lat4[0]<=glat4)[0][0]
iend=nonzero(lat4[-1]>=glat4)[0][-1]
h4=h4[istart:iend]

#set up map projection
lon_0=-45. #center coordinates of the plot --> polarstereographic with greenland pointing downwards
lat_0=90.0
####cut region MELTEX
latmin,lonmin=69,-135.
latmax,lonmax=73,-146.

Xmax,Ymax=mapll(latmax,lonmax,1)
Xmin,Ymin=mapll(latmin,lonmin,1)
###map set-up
Tmap=Basemap(resolution='h',projection='stere', boundinglat=60.,llcrnrlon=lonmin, llcrnrlat=latmin, urcrnrlon=lonmax, urcrnrlat=latmax, lon_0=lon_0, lat_0=lat_0)



############################################
#read AMSRE data
####---snow1:
AMSREfolder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/GRID_AMSRE/'
AMSREfile4n=AMSREfolder+'AMSR_E_L3_SeaIce12km_V12_20080604.hdf'
AMSREfile7n=AMSREfolder+'AMSR_E_L3_SeaIce12km_V12_20080607.hdf'

def XYgrid2(region,cs):
    x0,y0,x2,y2=x0y0x2y2(region)
    y=linspace(y2+cs/2,y0-cs/2,(y0-y2)/cs)
    x=linspace(x2+cs/2,x0-cs/2,(x0-x2)/cs)
    return x,y


###NSIDC NT2 (Bootstrap key [29])
sd_amsre_data4n=SD.SD(AMSREfile4n)
k4=sd_amsre_data4n.datasets().keys()[26]
icecon4=sd_amsre_data4n.select(k4)
Aice4n=icecon4.get()[::-1]

sd_amsre_data7n=SD.SD(AMSREfile7n)
k7=sd_amsre_data7n.datasets().keys()[26]
icecon7=sd_amsre_data7n.select(k7)
Aice7n=icecon7.get()[::-1]

figure()
imshow(Aice4n)



region='Arc'
cs=12.5
AX,AY=XYgrid2(region,cs)
x0,x2=AX.max(),AX.min()
y0,y2=AY.max(),AY.min()


iXmax=nonzero(AX>=Xmax)[0][0]
iXmin=nonzero(AX>=Xmin)[0][0]
iYmax=nonzero(AY>=Ymax)[0][0]
iYmin=nonzero(AY>=Ymin)[0][0]

#iXmax=304
#iXmin=0
#iYmax=448
#iYmin=0

AX=AX[iXmin:iXmax]
AY=AY[iYmin:iYmax]
Aice4n=Aice4n[iYmin:iYmax,iXmin:iXmax]*100
Aice7n=Aice7n[iYmin:iYmax,iXmin:iXmax]*100

figure()
imshow(Aice4n)



#plot the image
figure()
Tmap.imshow(Aice4n, interpolation='nearest',aspect='auto',cmap='bone')
x4,y4=Tmap(lon4_10,lat4_10)
Tmap.scatter(x4,y4,s=20,c=ice4_filtered10,cmap='bone', marker='o',vmin=0.,vmax=100.,edgecolors='black')
colorbar()
title('AMSRE - NT2 Ice Concentration 4.6.2008')
Tmap.drawcoastlines(color='black')
Tmap.drawmeridians(np.arange(0,360,10),labels=[1,0,0,0])
Tmap.drawparallels(np.arange(-90,90,5),labels=[0,0,0,1])
text(580000.0,300000.0,'ice concentration [%]', rotation='vertical',fontsize='x-large')
savefig('/pf/u/u241127/plots/meltex_amsre_iceconc_ice4_nt2.png')

figure()
Tmap.imshow(Aice7n, interpolation='nearest',aspect='auto',cmap='bone')
x7,y7=Tmap(lon7_10,lat7_10)
Tmap.scatter(x7,y7,s=20,c=ice7_filtered10,cmap='bone', marker='o',vmin=0.,vmax=100.,edgecolors='black')
colorbar()
title('AMSRE - NT2 Ice Concentration 7.6.2008')
Tmap.drawcoastlines(color='black')
Tmap.drawmeridians(np.arange(0,360,10),labels=[1,0,0,0])
Tmap.drawparallels(np.arange(-90,90,5),labels=[0,0,0,1])
text(580000.0,300000.0,'ice concentration [%]', rotation='vertical',fontsize='x-large')
savefig('/pf/u/u241127/plots/meltex_amsre_iceconc_ice7_nt2.png')


###ASI Algo
AMSREfile4=AMSREfolder+'asi-n6250-20080604-v5i.hdf'
AMSREfile7=AMSREfolder+'asi-n6250-20080607-v5i.hdf'

sd_amsre_data4=SD.SD(AMSREfile4)
k4=sd_amsre_data4.datasets().keys()[0]
icecon4=sd_amsre_data4.select(k4)
Aice4=icecon4.get()

sd_amsre_data7=SD.SD(AMSREfile7)
k7=sd_amsre_data7.datasets().keys()[0]
icecon7=sd_amsre_data7.select(k7)
Aice7=icecon7.get()

figure()
imshow(Aice7)

region='Arc'
cs=6.25
AX,AY=XYgrid2(region,cs)
x0,x2=AX.max(),AX.min()
y0,y2=AY.max(),AY.min()


iXmax=nonzero(AX>=Xmax)[0][0]
iXmin=nonzero(AX>=Xmin)[0][0]
iYmax=nonzero(AY>=Ymax)[0][0]
iYmin=nonzero(AY>=Ymin)[0][0]

#iXmax=608
#iXmin=0
#iYmax=896
#iYmin=0
AX=AX[iXmin:iXmax]
AY=AY[iYmin:iYmax]
Aice4=Aice4[iYmin:iYmax,iXmin:iXmax]*100
Aice7=Aice7[iYmin:iYmax,iXmin:iXmax]*100

figure()
imshow(Aice4)

#plot the image
figure()
Tmap.imshow(Aice4, interpolation='nearest',aspect='auto',cmap='bone')
x4,y4=Tmap(lon4_10,lat4_10)
Tmap.scatter(x4,y4,s=20,c=ice4_filtered10,cmap='bone', marker='o',vmin=0.,vmax=100.,edgecolors='black')
colorbar()
title('AMSRE - ASI Ice Concentration 4.6.2008')
Tmap.drawcoastlines(color='black')
Tmap.drawmeridians(np.arange(0,360,10),labels=[1,0,0,0])
Tmap.drawparallels(np.arange(-90,90,5),labels=[0,0,0,1])
text(580000.0,300000.0,'ice concentration [%]', rotation='vertical',fontsize='x-large')
savefig('/pf/u/u241127/plots/meltex_amsre_iceconc_ice4_asi.png')

figure()
Tmap.imshow(Aice7, interpolation='nearest',aspect='auto',cmap='bone')
x7,y7=Tmap(lon7_10,lat7_10)
Tmap.scatter(x7,y7,s=20,c=ice7_filtered10,cmap='bone', marker='o',vmin=0.,vmax=100.,edgecolors='black')
colorbar()
title('AMSRE - ASI Ice Concentration 7.6.2008')
Tmap.drawcoastlines(color='black')
Tmap.drawmeridians(np.arange(0,360,10),labels=[1,0,0,0])
Tmap.drawparallels(np.arange(-90,90,5),labels=[0,0,0,1])
text(580000.0,300000.0,'ice concentration [%]', rotation='vertical',fontsize='x-large')
savefig('/pf/u/u241127/plots/meltex_amsre_iceconc_ice7_asi.png')

############################################
#read MODIS data
####---snow1:
folder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/PRODUCTS/daily/'
plotfolder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/PRODUCTS/daily/plots/'



#modis 5 june 2008
sddataset5=folder+'2008_157/2008_157_mp_125SD.nc'
SD_obj5=NetCDFFile(sddataset5)
SD5=array(SD_obj5.variables['z'][:])
dataset5=folder+'2008_157/2008_157_ice_125.nc'
ic_obj5=NetCDFFile(dataset5)
ICE5=array(ic_obj5.variables['z'][:])
#koordinaten
ICEX=array(ic_obj5.variables['x'][:])
ICEY=array(ic_obj5.variables['y'][:])


#modis 7 june 2008
sd_mod_7=folder+'2008_159/2008_159_mp_125SD.nc'
SD_obj7=NetCDFFile(sd_mod_7)
SD7=array(SD_obj7.variables['z'][:])
icedataset7=folder+'2008_159/2008_159_ice_125.nc'
ICE_obj7=NetCDFFile(icedataset7)
ICE7=array(ICE_obj7.variables['z'][:])


####cut region MELTEX
x0,x2=max(ICEX),min(ICEX)
y0,y2=max(ICEY),min(ICEY)
#set up map projection
lon_0=-45. #center coordinates of the plot --> polarstereographic with greenland pointing downwards
lat_0=90.0

latmin,lonmin=69.,-135.
latmax,lonmax=73.,-146.

Xmax,Ymax=mapll(latmax,lonmax,1)
Xmin,Ymin=mapll(latmin,lonmin,1)

iXmax=nonzero(ICEX>=Xmax)[0][0]
iXmin=nonzero(ICEX>=Xmin)[0][0]
iYmax=nonzero(ICEY>=Ymax)[0][0]
iYmin=nonzero(ICEY>=Ymin)[0][0]

ICEX=ICEX[iXmin:iXmax]
ICEY=ICEY[iYmin:iYmax]
ICE5=ICE5[iYmin:iYmax,iXmin:iXmax]*100.
SD5=SD5[iYmin:iYmax,iXmin:iXmax]*100.
ICE7=ICE7[iYmin:iYmax,iXmin:iXmax]*100.
SD7=SD7[iYmin:iYmax,iXmin:iXmax]*100.


#plot the image
figure()
Tmap.imshow(ICE7, interpolation='nearest',cmap='bone')
#colorbar(orientation='horizontal',shrink=0.5)
x7,y7=Tmap(lon7_10,lat7_10)
Tmap.scatter(x7,y7,s=20,c=ice7_filtered10,cmap='bone', marker='o',vmin=0.,vmax=100.,edgecolors= 'black')
colorbar()
title('MODIS & MELTEX ICE CONCENTRATION 7.6.2008')
Tmap.drawcoastlines(color='black')
Tmap.drawmeridians(np.arange(0,360,10),labels=[1,0,0,0])
Tmap.drawparallels(np.arange(-90,90,5),labels=[0,0,0,1])
text(580000.0,300000.0,'ice concentration [%]', rotation='vertical',fontsize='large')
#text(210000.0,-150000.0,'ice concentration [%]', fontsize='large')
savefig('/pf/u/u241127/plots/meltex_modis_onlyice_7_6_125.png')


figure()
Tmap.imshow(ICE5, interpolation='nearest',cmap='bone')
#colorbar(orientation='horizontal',shrink=0.5)
x4,y4=Tmap(lon4_10,lat4_10)
Tmap.scatter(x4,y4,s=20,c=ice4_filtered10,cmap='bone', marker='o',vmin=0.,vmax=100.,edgecolors= 'black')
colorbar()
title('MODIS & MELTEX ICE CONCENTRATION 4.6.2008')
Tmap.drawcoastlines(color='black')
Tmap.drawmeridians(np.arange(0,360,10),labels=[1,0,0,0])
Tmap.drawparallels(np.arange(-90,90,5),labels=[0,0,0,1])
text(580000.0,300000.0,'ice concentration [%]', rotation='vertical',fontsize='large')
#text(210000.0,-150000.0,'ice concentration [%]', fontsize='large')
savefig('/pf/u/u241127/plots/meltex_modis_onlyice_4_6_125.png')


##########################
#######################################
#PROFILE MODIS
xmin=ICEX.min()
ymin=ICEY.min()
cs=12.5
ix7=[]
iy7=[]
for n in range(Xair710.shape[0]):
    x=int(abs(xmin+cs/2-Xair710[n])/cs)
    ix7.append(x+1)


for n in range(Yair710.shape[0]):
    y=int(abs(ymin+cs/2-Yair710[n])/cs)
    iy7.append(y+1)


ix7=array(ix7)
iy7=array(iy7)


#Test
ICEXi7=[]
for n in range(Yair710.shape[0]):
    x=ICEX[ix7[n]]
    ICEXi7.append(x)

ICEXi7=array(ICEXi7)




ICEi7=[]
for n in range(Yair710.shape[0]):
    mp=ICE7[iy7[n]][ix7[n]]
    ICEi7.append(mp)

ICEi7=array(ICEi7)

SDi7=[]
for n in range(Yair710.shape[0]):
    sd=SD7[iy7[n]][ix7[n]]
    SDi7.append(sd)

SDi7=array(SDi7)

#########height calc
n_glat7_track=len(glat7_track)
n_ice7_filtered10=len(ice7_filtered10)

f=n_glat7_track/n_ice7_filtered10
h7=h7[::f]



x7=linspace(0.,ice7_filtered10.shape[0],ice7_filtered10.shape[0])





########################twinx plot
figure()
ax1 = subplot(111)
errorbar(x7,ice7_filtered10,yerr=diff7,fmt='go',label='MELTEX w. STD')
errorbar(x7,ICEi7,yerr=SDi7,fmt='ro',label='MODIS w. STD')
axis([0,100,0,100])
xlabel('# of pixel')
title('7.6.2008')
ax1.set_ylabel('sea ice concentration  [%]')
legend(numpoints=1)
# 2nd y-axes 
ax2 = twinx()
plot(x7, h7[:-1], 'b-')
#ax2.yaxis.tick_right()
ax2.set_ylabel('height [m]',color='b')
for tl in ax2.get_yticklabels():
    tl.set_color('b')

##########RMSE calc
n=ICEi7.shape[0]
RMSE7=round(sqrt((1./n)*nansum((ice7_filtered10-ICEi7)**2)),1)
print 'RMSE7: '+str(RMSE7)
annotate('RMSE='+str(RMSE7)+'%',xy=(3,3300))

savefig('/pf/u/u241127/plots/meltex_7_6_2008_profile_height_ic.png')





###4.6.
ix4=[]
iy4=[]
for n in range(Xair410.shape[0]):
    x=int(abs(xmin+cs/2-Xair410[n])/cs)
    ix4.append(x+1)


for n in range(Yair410.shape[0]):
    y=int(abs(ymin+cs/2-Yair410[n])/cs)
    iy4.append(y+1)


ix4=array(ix4)
iy4=array(iy4)

ICEi4=[]
for n in range(Yair410.shape[0]):
    mp=ICE5[iy4[n]][ix4[n]]
    ICEi4.append(mp)

ICEi4=array(ICEi4)

SDi4=[]
for n in range(Yair410.shape[0]):
    sd=SD5[iy4[n]][ix4[n]]
    SDi4.append(sd)

SDi4=array(SDi4)


#########height calc
n_h4=len(h4)
n_ice4_filtered10=len(ice4_filtered10)

f=n_h4/n_ice4_filtered10
h4=h4[::f][1:]

#figure()
#plot(h4)
#title('aircraft height 4.6.2008')
#savefig('/pf/u/u241127/plots/meltex_4_6_2008_height.png')


x4=linspace(0.,mp4_filtered10.shape[0],mp4_filtered10.shape[0])


figure()
ax1 = subplot(111)
errorbar(x4,ice4_filtered10,yerr=diff4,fmt='go',label='MELTEX w. STD')
errorbar(x4,ICEi4,yerr=SDi4,fmt='ro',label='MODIS w. STD')
axis([0,100,0,100])
xlabel('# of pixel')
title('4.6.2008')
ax1.set_ylabel('sea ice concentration [%]')
legend(numpoints=1)
# 2nd y-axes 
ax2 = twinx()
plot(x4, h4, 'b-')
#ax2.yaxis.tick_right()
ax2.set_ylabel('height [m]',color='b')
for tl in ax2.get_yticklabels():
    tl.set_color('b')

##########RMSE calc
n=ICEi4.shape[0]
RMSE4=round(sqrt((1./n)*nansum((ice4_filtered10-ICEi4)**2)),1)

print 'RMSE4: '+str(RMSE4)
annotate('RMSE='+str(RMSE4)+'%',xy=(3,2350))

savefig('/pf/u/u241127/plots/meltex_4_6_2008_profile_height_ic.png')

##########PROFILE AMSR-E##################profile amsre
##########################################
#4.6.
xmin=AX.min()
ymin=AY.min()
cs=6.25
iAx4=[]
iAy4=[]
for n in range(Xair410.shape[0]):
    x=int(abs(xmin+cs/2-Xair410[n])/cs)
    iAx4.append(x+1)


for n in range(Yair410.shape[0]):
    y=int(abs(ymin+cs/2-Yair410[n])/cs)
    iAy4.append(y+1)


iAx4=array(iAx4)
iAy4=array(iAy4)

ICEAi4=[]
for n in range(Yair410.shape[0]):
    mp=Aice4[iy4[n]][ix4[n]]
    ICEAi4.append(mp)

ICEAi4=array(ICEAi4)



########################twinx plot
figure()
ax1 = subplot(111)
errorbar(x4,ice4_filtered10,yerr=diff4,fmt='go',label='MELTEX w. STD')
plot(x4,ICEAi4,'ro',label='AMSR-E')
axis([0,100,0,100])
xlabel('# of pixel')
title('4.6.2008')
ax1.set_ylabel('ice concentration [%]')
legend(numpoints=1, loc=4)
# 2nd y-axes 
ax2 = twinx()
plot(x4, h4, 'b-')
#ax2.yaxis.tick_right()
ax2.set_ylabel('height [m]',color='b')
for tl in ax2.get_yticklabels():
    tl.set_color('b')
##########RMSE calc
n=ICEAi4.shape[0]
RMSE4=round(sqrt((1./n)*nansum((ice4_filtered10-ICEAi4)**2)),1)
print 'RMSE7: '+str(RMSE4)
annotate('RMSE='+str(RMSE4)+'%',xy=(3,3300))
savefig('/pf/u/u241127/plots/meltex_7_6_2008_profile_height_ic_amsre.png')


#7.6.
xmin=AX.min()
ymin=AY.min()
cs=6.25
iAx7=[]
iAy7=[]
for n in range(Xair710.shape[0]):
    x=int(abs(xmin+cs/2-Xair710[n])/cs)
    iAx7.append(x+1)


for n in range(Yair710.shape[0]):
    y=int(abs(ymin+cs/2-Yair710[n])/cs)
    iAy7.append(y+1)


iAx7=array(iAx7)
iAy7=array(iAy7)

"""
#Test
ICEAXi7=[]
for n in range(Yair710.shape[0]):
    x=ICEAX[iAx7[n]]
    ICEAXi7.append(x)

ICEAXi7=array(ICEAXi7)

ICEAYi7=[]
for n in range(Yair710.shape[0]):
    y=ICEAY[iAy7[n]]
    ICEAYi7.append(y)

ICEAYi7=array(ICEAYi7)
"""

ICEAi7=[]
for n in range(Yair710.shape[0]):
    mp=Aice7[iy7[n]][ix7[n]]
    ICEAi7.append(mp)

ICEAi7=array(ICEAi7)


ice7_filtered=ndi.gaussian_filter1d(ice7,10)
ice7_filtered10=ice7_filtered[range(0,ice7_filtered.shape[0],10)] #nur jeden 20. wert


########################twinx plot
figure()
ax1 = subplot(111)
errorbar(x7,ice7_filtered10,yerr=diff7,fmt='go',label='MELTEX w. STD')
plot(x7,ICEAi7,'ro',label='AMSR-E')
axis([0,100,0,100])
xlabel('# of pixel')
title('7.6.2008')
ax1.set_ylabel('ice concentration [%]')
legend(numpoints=1)
# 2nd y-axes 
ax2 = twinx()
plot(x7, h7[:-1], 'b-')
#ax2.yaxis.tick_right()
ax2.set_ylabel('height [m]',color='b')
for tl in ax2.get_yticklabels():
    tl.set_color('b')

##########RMSE calc
n=ICEAi7.shape[0]
RMSE7=round(sqrt((1./n)*nansum((ice7_filtered10-ICEAi7)**2)),1)
print 'RMSE7: '+str(RMSE7)
annotate('RMSE='+str(RMSE7)+'%',xy=(3,3300))
savefig('/pf/u/u241127/plots/meltex_7_6_2008_profile_height_ic_amsre.png')






##means
meanMOD4=str(mean(ICEi4))
meanMELTEX4=str(mean(ice4_filtered10))
meanMOD7=str(mean(ICEi7))
meanMELTEX7=str(mean(ice7_filtered10))

print('4.6. mean MODIS: '+meanMOD4)
print('4.6. mean MELTEX: '+meanMELTEX4)
print('7.6. mean MODIS: '+meanMOD7)
print('7.6. mean MELTEX: '+meanMELTEX7)

