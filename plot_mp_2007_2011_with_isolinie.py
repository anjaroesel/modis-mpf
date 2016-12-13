from scipy import *
import sys,pipes,struct,os,glob
from pylab import *
from mpl_toolkits.basemap import Basemap
from polar_projection import *
import Nio
from sattools import *
from gmttools import *
import time,glob,os,os.path,struct
from mpl_toolkits.basemap import Basemap
import scipy.io as sio
import scipy.ndimage as ndi
import scipy.stats as stats

from gmt_tools import *
from pylab import *
from pyhdf.SD import *
from grids import def_regions

# run plot_mp_2007_2011_with_isolinie
# python plot_mp_2007_2011_with_isolinie

def JulDay2Date(year,n):
    import datetime
    d = datetime.date(year, 1, 1) + datetime.timedelta(n - 1)
    month=d.month
    day=d.day
    return month, day




#read MODIS mp data
####---snow1:
folder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/PRODUCTS/'
plotfolder='/pf/u/u241127/plots/'

pathfilename_ct07='/scratch/clisap/seaice/SSMI_ASI/CERSAT/Median5/Arctic_median5_2007.nc'
pathfilename_ct11='/scratch/clisap/seaice/SSMI_ASI/CERSAT/Median5/Arctic_median5_2011_0101_1110.nc'


# Ice concentration
fid07=sio.netcdf_file(pathfilename_ct07,'r')
ic07=fid07.variables['concentration'][262,:,:].astype(float32)
ic07[ic07<0]=nan
ic07[ic07>100]=nan
ic07=ic07[::-1]


fid11=sio.netcdf_file(pathfilename_ct11,'r')
ic11=fid11.variables['concentration'][262,:,:].astype(float32)
ic11[ic11<0]=nan
ic11[ic11>100]=nan
ic11=ic11[::-1]


#2007
workfile=folder+'2007_169/2007_169_mp_125masked.nc'
pic1 = NetCDFFile(workfile,'r')
day169 = array(pic1.variables['z'][:,:])
X = array(pic1.variables['x'][:])
Y = array(pic1.variables['y'][:])

####map setup
lon_0=-45. #center coordinates of the plot
lat_0=90.0

####cut plot region 
latmin,lonmin=55,-90.
latmax,lonmax=55,90.

Xmax,Ymax=mapll(latmax,lonmax,1)
Xmin,Ymin=mapll(latmin,lonmin,1)

iXmax=nonzero(X>=Xmax)[0][0]
iXmin=nonzero(X>=Xmin)[0][0]
iYmax=nonzero(Y>=Ymax)[0][0]
iYmin=nonzero(Y>=Ymin)[0][0]

X=X[iXmin:iXmax]
Y=Y[iYmin:iYmax]

lats, lons = mapxy(X,Y,1)
day169=day169[iYmin:iYmax,iXmin:iXmax]
#MP05=MP05[iYmin:iYmax,iXmin:iXmax]
ic07=ic07[iYmin:iYmax,iXmin:iXmax]

Tmap=Basemap(resolution='l',projection='stere', boundinglat=70.,llcrnrlon=lonmin, llcrnrlat=latmin, urcrnrlon=lonmax, urcrnrlat=latmax, lon_0=lon_0, lat_0=lat_0)


m=169
year=2007

m=209
year=2003

figure(figsize=(7,7))
Tmap.drawcoastlines()
Tmap.fillcontinents(color='gray')
Tmap.drawmeridians(np.arange(0,360,40),labels=[0,0,0,0])  #labels=[1,0,0,0]
Tmap.drawparallels(np.arange(-90,90,10),labels=[0,0,0,0])  #labels=[0,0,0,1]
#title('day '+m)
#conturplot
levels = [15]
x,y=Tmap(lons,lats)
CS = contour(x, y, ic07, levels, origin='lower',colors = ('k',),linewidths = (2,))

month,day=JulDay2Date(int(year),int(m))
text(2000000.0,5200000.0,str(month)+'.'+str(day)+'.'+str(year), fontsize='28', color='white')
Tmap.imshow(day169, interpolation='nearest',vmin=0.,vmax=0.5)
#Tmap.imshow(MP05, interpolation='nearest',vmin=0.,vmax=0.5)
cax = axes([0.2, 0.06, 0.6, 0.03])
colorbar(cax=cax,orientation='horizontal',ticks=[0.0,0.1,0.2,0.3,0.4,0.5],extend='max') 
savefig('/pf/u/u241127/plots/mp_oneyearplt'+str(year)+str(day)+'masked_w_ISO.png')
    
    
    

#2011
workfile=folder+'2011_169/2011_169_mp_125masked.nc'
pic1 = NetCDFFile(workfile,'r')
day169 = array(pic1.variables['z'][:,:])
X = array(pic1.variables['x'][:])
Y = array(pic1.variables['y'][:])

####map setup
lon_0=-45. #center coordinates of the plot
lat_0=90.0

####cut plot region 
latmin,lonmin=55,-90.
latmax,lonmax=55,90.

Xmax,Ymax=mapll(latmax,lonmax,1)
Xmin,Ymin=mapll(latmin,lonmin,1)

iXmax=nonzero(X>=Xmax)[0][0]
iXmin=nonzero(X>=Xmin)[0][0]
iYmax=nonzero(Y>=Ymax)[0][0]
iYmin=nonzero(Y>=Ymin)[0][0]

X=X[iXmin:iXmax]
Y=Y[iYmin:iYmax]

lats, lons = mapxy(X,Y,1)
day169=day169[iYmin:iYmax,iXmin:iXmax]
ic11=ic11[iYmin:iYmax,iXmin:iXmax]

Tmap=Basemap(resolution='l',projection='stere', boundinglat=70.,llcrnrlon=lonmin, llcrnrlat=latmin, urcrnrlon=lonmax, urcrnrlat=latmax, lon_0=lon_0, lat_0=lat_0)


m=169
year=2011
figure(figsize=(7,7))
Tmap.drawcoastlines()
Tmap.fillcontinents(color='gray')
Tmap.drawmeridians(np.arange(0,360,40),labels=[0,0,0,0])  #labels=[1,0,0,0]
Tmap.drawparallels(np.arange(-90,90,10),labels=[0,0,0,0])  #labels=[0,0,0,1]
#title('day '+m)
#conturplot
levels = [15]
x,y=Tmap(lons,lats)
CS = contour(x, y, ic11, levels, origin='lower',colors = ('k',),linewidths = (2,))

month,day=JulDay2Date(int(year),int(m))
text(2000000.0,5200000.0,str(month)+'.'+str(day)+'.'+str(year), fontsize='28', color='white')
Tmap.imshow(day169, interpolation='nearest',vmin=0.,vmax=0.5)
cax = axes([0.2, 0.06, 0.6, 0.03])
colorbar(cax=cax,orientation='horizontal',ticks=[0.0,0.1,0.2,0.3,0.4,0.5],extend='max') 
savefig('/pf/u/u241127/plots/mp_oneyearplt'+str(year)+str(day)+'masked_w_ISO.png')






