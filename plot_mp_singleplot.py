from scipy import *
import sys,pipes,struct,os,glob
from pylab import *
from mpl_toolkits.basemap import Basemap
from polar_projection import *
import Nio
from sattools import *
from gmttools import *

# run plot_one_year_mp_on_one_page.py
# python plot_one_year_mp_on_one_page.py

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


date='2003_209'
year=int(date[:4])
n=int(date[5:])

workfile=folder+date+'/'+date+'_mp_125masked.nc'
pic1 = NetCDFFile(workfile,'r')
img = array(pic1.variables['z'][:,:])
d[date]=[]
d[date].append(img)


datasets= sort(d.keys())
X = array(pic1.variables['x'][:])
Y = array(pic1.variables['y'][:])
I=d[datasets[0]][0]

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


I=I[iYmin:iYmax,iXmin:iXmax]


Tmap=Basemap(resolution='l',projection='stere', boundinglat=70.,llcrnrlon=lonmin, llcrnrlat=latmin, urcrnrlon=lonmax, urcrnrlat=latmax, lon_0=lon_0, lat_0=lat_0)



    
#####PLOT
month,day=JulDay2Date(int(year),int(n))
figure(figsize=(5,5))
Tmap.drawcoastlines()
Tmap.fillcontinents(color='gray')
Tmap.drawmeridians(np.arange(0,360,40),labels=[0,0,0,0])  #labels=[1,0,0,0]
Tmap.drawparallels(np.arange(-90,90,10),labels=[0,0,0,0])  #labels=[0,0,0,1]
Tmap.imshow(I, interpolation='nearest',vmin=0.,vmax=0.5)
text(2000000.0,5200000.0,str(month)+'.'+str(day)+'.'+str(year), fontsize='xx-large', color='white')

subplots_adjust(bottom=0.1, right=0.9, top=0.9,left=0.1)
cax = axes([0.2, 0.04, 0.6, 0.03])
colorbar(cax=cax, orientation='horizontal',ticks=[0.0,0.1,0.2,0.3,0.4,0.5],extend='max')  #ticks=[0.0,0.08,0.16,0.24,0.32,0.40]
#text(0.3,-4.5,'melt pond fraction',fontsize='large',color='black')
savefig('/pf/u/u241127/plots/mp_singleplt_masked'+str(year)+str(n)+'.png')

