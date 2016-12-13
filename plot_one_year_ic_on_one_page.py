from scipy import *
import sys,pipes,struct,os,glob
from pylab import *
from mpl_toolkits.basemap import Basemap
from polar_projection import *
import Nio
from sattools import *
from gmttools import *

#run plot_one_year_ic_on_one_page.py


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

years=['2000','2002','2003','2004','2005','2006','2007','2008','2009','2010','2011']
#years=['2011']

for year in years:
    print year
    #year='2003'
    metalist= glob.glob(folder+year+'_*')
    d={}
    for fn in metalist:
        print fn
        date=fn[-8:]
        workfile=fn+'/'+date+'_ow_125.nc'
        pic1 = NetCDFFile(workfile,'r')
        img = array(pic1.variables['z'][:,:])
        d[date]=[]
        d[date].append(img)


    datasets= sort(d.keys())
    X = array(pic1.variables['x'][:])
    Y = array(pic1.variables['y'][:])

    #day257=d[datasets[16]][0]
    day249=d[datasets[15]][0]
    day241=d[datasets[14]][0]
    day233=d[datasets[13]][0]
    day225=d[datasets[12]][0]
    day217=d[datasets[11]][0]
    day209=d[datasets[10]][0]
    day201=d[datasets[9]][0]
    day193=d[datasets[8]][0]
    day185=d[datasets[7]][0]
    day177=d[datasets[6]][0]
    day169=d[datasets[5]][0]
    day161=d[datasets[4]][0]
    day153=d[datasets[3]][0]
    day145=d[datasets[2]][0]
    day137=d[datasets[1]][0]
    day129=d[datasets[0]][0]

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

    #ice konz = (1-ow)*100
    #day257=(1-day257[iYmin:iYmax,iXmin:iXmax])*100
    day249=(1-day249[iYmin:iYmax,iXmin:iXmax])*100
    day241=(1-day241[iYmin:iYmax,iXmin:iXmax])*100
    day233=(1-day233[iYmin:iYmax,iXmin:iXmax])*100
    day225=(1-day225[iYmin:iYmax,iXmin:iXmax])*100
    day217=(1-day217[iYmin:iYmax,iXmin:iXmax])*100
    day209=(1-day209[iYmin:iYmax,iXmin:iXmax])*100
    day201=(1-day201[iYmin:iYmax,iXmin:iXmax])*100
    day193=(1-day193[iYmin:iYmax,iXmin:iXmax])*100
    day185=(1-day185[iYmin:iYmax,iXmin:iXmax])*100
    day177=(1-day177[iYmin:iYmax,iXmin:iXmax])*100
    day169=(1-day169[iYmin:iYmax,iXmin:iXmax])*100
    day161=(1-day161[iYmin:iYmax,iXmin:iXmax])*100
    day153=(1-day153[iYmin:iYmax,iXmin:iXmax])*100
    day145=(1-day145[iYmin:iYmax,iXmin:iXmax])*100
    day137=(1-day137[iYmin:iYmax,iXmin:iXmax])*100
    day129=(1-day129[iYmin:iYmax,iXmin:iXmax])*100

    Tmap=Basemap(resolution='l',projection='stere', boundinglat=70.,llcrnrlon=lonmin, llcrnrlat=latmin, urcrnrlon=lonmax, urcrnrlat=latmax, lon_0=lon_0, lat_0=lat_0)

    liste=[day129,day137,day145,day153,day161,day169,day177,day185,day193,day201,day209,day217,day225,day233,day241,day249]
    liste2=['129','137','145','153','161','169','177','185','193','201','209','217','225','233','241','249']

    year=int(year)
    figure(figsize=(10,10))
    axisNum = 0
    for n in liste:
        axisNum += 1
        #month,day=JulDay2Date(int(year),int(m))
        subplot(4, 4, axisNum)
        Tmap.drawcoastlines()
        Tmap.fillcontinents(color='gray')
        Tmap.drawmeridians(np.arange(0,360,40),labels=[0,0,0,0])  #labels=[1,0,0,0]
        Tmap.drawparallels(np.arange(-90,90,10),labels=[0,0,0,0])  #labels=[0,0,0,1]
        #title('day '+m)
        Tmap.imshow(n, interpolation='nearest',vmin=0.,vmax=100,cmap='bone')
    
    
    axisNum = 0
    for m in liste2:
        axisNum += 1
        subplot(4, 4, axisNum)
        month,day=JulDay2Date(int(year),int(m))
        #title(str(day)+'.'+str(month)+'.')
        text(1400000.0,5000000.0,str(month)+'.'+str(day)+'.'+str(year), fontsize='large', color='white')

    subplots_adjust(bottom=0.1, right=0.9, top=0.9,left=0.1)
    cax = axes([0.3, 0.06, 0.4, 0.01])
    colorbar(cax=cax, orientation='horizontal',ticks=[0,20,40,60,80,100])
    text(0.2,-4.5,'sea ice concentration in %',fontsize='large',color='black')
    savefig('/pf/u/u241127/plots/ic_oneyearplt'+str(year)+'.png')


#######################################
###########year2001 #day169 is missing#
#######################################
year='2001'
metalist= glob.glob(folder+year+'_*')
d={}
for fn in metalist:
    print fn
    date=fn[-8:]
    workfile=fn+'/'+date+'_ow_125.nc'
    pic1 = NetCDFFile(workfile,'r')
    img = array(pic1.variables['z'][:,:])
    d[date]=[]
    d[date].append(img)


datasets= sort(d.keys())
X = array(pic1.variables['x'][:])
Y = array(pic1.variables['y'][:])


day249=d[datasets[13]][0]
day241=d[datasets[12]][0]
day233=d[datasets[11]][0]
day225=d[datasets[10]][0]
day217=d[datasets[9]][0]
day209=d[datasets[8]][0]
day201=d[datasets[7]][0]
day193=d[datasets[6]][0]
day185=d[datasets[5]][0]
#day177=d[datasets[5]][0]
#day169=d[datasets[5]][0]
day161=d[datasets[4]][0]
day153=d[datasets[3]][0]
day145=d[datasets[2]][0]
day137=d[datasets[1]][0]
day129=d[datasets[0]][0]

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



#day257=day257[iYmin:iYmax,iXmin:iXmax]
day249=(1-day249[iYmin:iYmax,iXmin:iXmax])*100
day241=(1-day241[iYmin:iYmax,iXmin:iXmax])*100
day233=(1-day233[iYmin:iYmax,iXmin:iXmax])*100
day225=(1-day225[iYmin:iYmax,iXmin:iXmax])*100
day217=(1-day217[iYmin:iYmax,iXmin:iXmax])*100
day209=(1-day209[iYmin:iYmax,iXmin:iXmax])*100
day201=(1-day201[iYmin:iYmax,iXmin:iXmax])*100
day193=(1-day193[iYmin:iYmax,iXmin:iXmax])*100
day185=(1-day185[iYmin:iYmax,iXmin:iXmax])*100
#day177=(1-day177[iYmin:iYmax,iXmin:iXmax])*100
#day169=day169[iYmin:iYmax,iXmin:iXmax]
day161=(1-day161[iYmin:iYmax,iXmin:iXmax])*100
day153=(1-day153[iYmin:iYmax,iXmin:iXmax])*100
day145=(1-day145[iYmin:iYmax,iXmin:iXmax])*100
day137=(1-day137[iYmin:iYmax,iXmin:iXmax])*100
day129=(1-day129[iYmin:iYmax,iXmin:iXmax])*100


Tmap=Basemap(resolution='l',projection='stere', boundinglat=70.,llcrnrlon=lonmin, llcrnrlat=latmin, urcrnrlon=lonmax, urcrnrlat=latmax, lon_0=lon_0, lat_0=lat_0)


liste=[day129,day137,day145,day153,day161,day185,day193,day201,
day209,day217,day225,day233,day241,day249]
liste2=['129','137','145','153','161','185','193','201','209','217','225','233','241','249']

year=int(year)
figure(figsize=(10,10))
axisNum = 0
for n in liste:
    axisNum += 1
    #month,day=JulDay2Date(int(year),int(m))
    subplot(4, 4, axisNum)
    Tmap.drawcoastlines()
    Tmap.fillcontinents(color='gray')
    Tmap.drawmeridians(np.arange(0,360,40),labels=[0,0,0,0])  #labels=[1,0,0,0]
    Tmap.drawparallels(np.arange(-90,90,10),labels=[0,0,0,0])  #labels=[0,0,0,1]
    #title('day '+m)
    Tmap.imshow(n, interpolation='nearest',vmin=0.,vmax=100,cmap='bone')
    
    
axisNum = 0
for m in liste2:
    axisNum += 1
    subplot(4, 4, axisNum)
    month,day=JulDay2Date(int(year),int(m))
    #title(str(day)+'.'+str(month)+'.')
    text(1400000.0,5000000.0,str(month)+'.'+str(day)+'.'+str(year),fontsize='large',color='white')

subplots_adjust(bottom=0.1, right=0.9, top=0.9,left=0.1)
cax = axes([0.3, 0.06, 0.4, 0.01])
colorbar(cax=cax, orientation='horizontal',ticks=[0,20,40,60,80,100])
text(0.2,-4.5,'sea ice concentration in %',fontsize='large',color='black')
savefig('/pf/u/u241127/plots/ic_oneyearplt'+str(year)+'.png')



