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

#run plot_AMSRE_MODIS_IC.py



def XYgrid2(region,cs):
    x0,y0,x2,y2=x0y0x2y2(region)
    y=linspace(y2+cs/2,y0-cs/2,(y0-y2)/cs)
    x=linspace(x2+cs/2,x0-cs/2,(x0-x2)/cs)
    return x,y


####map setup
lon_0=-45. #center coordinates of the plot
lat_0=90.0

####cut plot region #canadian achipelago
latmin,lonmin=65,-100.
latmax,lonmax=80,-110.

##Meltex
#latmin,lonmin=69.,-135.
#latmax,lonmax=73.,-146.


Xmax,Ymax=mapll(latmax,lonmax,1)
Xmin,Ymin=mapll(latmin,lonmin,1)


############################################
#read MODIS data
####---snow1:
folder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/PRODUCTS/'

icedataset=folder+'2011_169/2011_169_ice_125.nc'
ICE_obj=NetCDFFile(icedataset)
ICE=array(ICE_obj.variables['z'][:])
ICEX=array(ICE_obj.variables['x'][:])
ICEY=array(ICE_obj.variables['y'][:])


owdataset=folder+'2011_169/2011_169_ow_125.nc'
OW_obj=NetCDFFile(owdataset)
OW=array(OW_obj.variables['z'][:])


iXmax=nonzero(ICEX>=Xmax)[0][0]
iXmin=nonzero(ICEX>=Xmin)[0][0]
iYmax=nonzero(ICEY>=Ymax)[0][0]
iYmin=nonzero(ICEY>=Ymin)[0][0]

ICEX=ICEX[iXmin:iXmax]
ICEY=ICEY[iYmin:iYmax]
ICE=ICE[iYmin:iYmax,iXmin:iXmax]*100.
OW=OW[iYmin:iYmax,iXmin:iXmax]*100.
ICE2=100-OW



mpdataset=folder+'2011_169/2011_169_mp_125.nc'
m_obj=NetCDFFile(mpdataset)
MP=array(m_obj.variables['z'][:])
MX=array(m_obj.variables['x'][:])
MY=array(m_obj.variables['y'][:])



Xmax,Ymax=mapll(latmax,lonmax,1)
Xmin,Ymin=mapll(latmin,lonmin,1)

iXmax=nonzero(MX>=Xmax)[0][0]
iXmin=nonzero(MX>=Xmin)[0][0]
iYmax=nonzero(MY>=Ymax)[0][0]
iYmin=nonzero(MY>=Ymin)[0][0]

MX=MX[iXmin:iXmax]
MY=MY[iYmin:iYmax]
MP=MP[iYmin:iYmax,iXmin:iXmax]*100.


#ICE3=ICE+MP

#figure()
#imshow(ICE2)
#colorbar()
#title('ICE2')
#figure()
#imshow(ICE)
#colorbar()
#title('ICE')
#figure()
#imshow(ICE3)
#colorbar()
#title('ICE3')
#figure()
#imshow(OW)
#colorbar()
#title('OW')

############################################
#read AMSRE data
####---snow1:
AMSREfolder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/GRID_AMSRE/'
AMSREfile=AMSREfolder+'asi-n6250-20110625-v5i.nc'
AMSREfile2=AMSREfolder+'AMSR_E_L3_SeaIce12km_V12_20110625.hdf'
##ASI
obj=NetCDFFile(AMSREfile)
asi=array(obj.variables['icecon'][:], dtype=float32)/10
asi[asi==-100]=nan


###NASA-TEAM2
AMSRE_NASA=SD.SD(AMSREfile2)
nt2_key=26#'SI_12km_NH_89V_DSC'
bt_key=29#'SI_12km_NH_ICECON_DSC'
sd_ice_nt2   = AMSRE_NASA.select(nt2_key)
ice_nt2      = sd_ice_nt2.get()
ice_nt2      = array(ice_nt2,dtype=float)[::-1]
ice_nt2[ice_nt2>=101]=nan


###NASA-BOOTSTRAP
sd_ice_bt   = AMSRE_NASA.select(bt_key)
ice_bt      = sd_ice_bt.get()    # stored as difference to NT2 data (ice_bst-ice_nt2)
ice_bt      = array(ice_bt,dtype=float)[::-1]
ice_bt      = ice_bt + ice_nt2
ice_bt[ice_bt>=101]=nan




region='Arc'
cs=6.25
AX,AY=XYgrid2(region,cs)
x0,x2=AX.max(),AX.min()
y0,y2=AY.max(),AY.min()

iXmax=nonzero(AX>=Xmax)[0][0]
iXmin=nonzero(AX>=Xmin)[0][0]
iYmax=nonzero(AY>=Ymax)[0][0]
iYmin=nonzero(AY>=Ymin)[0][0]

AX=AX[iXmin:iXmax]
AY=AY[iYmin:iYmax]
asi=asi[iYmin:iYmax,iXmin:iXmax]

region='Arc'
cs=12.5
AX,AY=XYgrid2(region,cs)
x0,x2=AX.max(),AX.min()
y0,y2=AY.max(),AY.min()

iXmax=nonzero(AX>=Xmax)[0][0]
iXmin=nonzero(AX>=Xmin)[0][0]
iYmax=nonzero(AY>=Ymax)[0][0]
iYmin=nonzero(AY>=Ymin)[0][0]

AX=AX[iXmin:iXmax]
AY=AY[iYmin:iYmax]
ice_nt2=ice_nt2[iYmin:iYmax,iXmin:iXmax]
ice_bt=ice_bt[iYmin:iYmax,iXmin:iXmax]


figure()
imshow(asi)
colorbar()

figure()
imshow(ice_nt2)
colorbar()

figure()
imshow(ice_bt)
colorbar()



Tmap=Basemap(resolution='h',projection='stere', boundinglat=60.,llcrnrlon=lonmin, llcrnrlat=latmin, urcrnrlon=lonmax, urcrnrlat=latmax, lon_0=lon_0, lat_0=lat_0)


figure()
subplot(2,3,1)
Tmap.imshow(asi, interpolation='nearest',aspect='auto',cmap='bone',vmin=0,vmax=100)
colorbar(orientation='horizontal', shrink=0.8, ticks=(0,20,40,60,80,100))
Tmap.drawcoastlines(color='black')
Tmap.fillcontinents(color='gray')
Tmap.drawmeridians(np.arange(0,360,10),labels=[1,0,0,0])
Tmap.drawparallels(np.arange(-90,90,5),labels=[0,0,0,1])
text(-00000.0,-550000.0,'ice concentration [%]')
text(180000.0,1250000.0,'AMSR-E ASI',fontsize='large')

subplot(2,3,2)
Tmap.imshow(ice_nt2, interpolation='nearest',aspect='auto',cmap='bone',vmin=0,vmax=100)
colorbar(orientation='horizontal', shrink=0.8, ticks=(0,20,40,60,80,100))
Tmap.drawcoastlines(color='black')
Tmap.fillcontinents(color='gray')
Tmap.drawmeridians(np.arange(0,360,10),labels=[1,0,0,0])
Tmap.drawparallels(np.arange(-90,90,5),labels=[0,0,0,1])
text(-00000.0,-550000.0,'ice concentration [%]')
text(180000.0,1250000.0,'AMSR-E NT2',fontsize='large')
subplot(2,3,3)
Tmap.imshow(ice_bt, interpolation='nearest',aspect='auto',cmap='bone',vmin=0,vmax=100)
colorbar(orientation='horizontal', shrink=0.8, ticks=(0,20,40,60,80,100))
Tmap.drawcoastlines(color='black')
Tmap.fillcontinents(color='gray')
Tmap.drawmeridians(np.arange(0,360,10),labels=[1,0,0,0])
Tmap.drawparallels(np.arange(-90,90,5),labels=[0,0,0,1])
text(-00000.0,-550000.0,'ice concentration [%]')
text(180000.0,1250000.0,'AMSR-E BT',fontsize='large')




subplot(2,3,4)
Tmap.imshow(ICE, interpolation='nearest',cmap='bone',vmin=0,vmax=100)
colorbar(orientation='horizontal', shrink=0.8, ticks=(0,20,40,60,80,100))
Tmap.drawcoastlines(color='black')
Tmap.fillcontinents(color='gray')
Tmap.drawmeridians(np.arange(0,360,10),labels=[1,0,0,0])
Tmap.drawparallels(np.arange(-90,90,5),labels=[0,0,0,1])
text(-00000.0,-550000.0,'ice concentration [%]')
text(180000.0,1250000.0,'MODIS',fontsize='large')
subplot(2,3,5)
Tmap.imshow(MP, interpolation='nearest',cmap='jet',vmin=0,vmax=70)
colorbar(orientation='horizontal', shrink=0.8, ticks=(0,20,40,60))
Tmap.drawcoastlines(color='black')
Tmap.fillcontinents(color='gray')
Tmap.drawmeridians(np.arange(0,360,10),labels=[1,0,0,0])
Tmap.drawparallels(np.arange(-90,90,5),labels=[0,0,0,1])
text(-00000.0,-550000.0,'melt pond fraction [%]')
text(180000.0,1250000.0,'MODIS',fontsize='large')
savefig('/pf/u/u241127/plots/amsre_modis_ic_comp_24_06_2008_allAMSRE.png')

"""
#plot the image
figure()
Tmap.imshow(AMSRE, interpolation='nearest',aspect='auto',cmap='bone')
colorbar()
title('AMSRE - ASI Ice Concentration 25.6.2008')
Tmap.drawcoastlines(color='black')
Tmap.fillcontinents(color='gray')
Tmap.drawmeridians(np.arange(0,360,10),labels=[1,0,0,0])
Tmap.drawparallels(np.arange(-90,90,5),labels=[0,0,0,1])


#plot the image
figure()
Tmap.imshow(ICE, interpolation='nearest',cmap='bone')
colorbar()
title('MODIS ICE C 25.6.2008')
Tmap.drawcoastlines(color='black')
Tmap.fillcontinents(color='gray')
Tmap.drawmeridians(np.arange(0,360,10),labels=[1,0,0,0])
Tmap.drawparallels(np.arange(-90,90,5),labels=[0,0,0,1])


#plot the image
figure()
Tmap.imshow(MP, interpolation='nearest',cmap='jet')
colorbar()
title('MODIS MP Frac 25.6.2008')
Tmap.drawcoastlines(color='black')
Tmap.fillcontinents(color='gray')
Tmap.drawmeridians(np.arange(0,360,10),labels=[1,0,0,0])
Tmap.drawparallels(np.arange(-90,90,5),labels=[0,0,0,1])

########savefig('/pf/u/u241127/plots/amsre_modis_ic_comp_24_06_2008.png')

"""
