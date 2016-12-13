from scipy import *
import sys,pipes,struct,os,glob
from pylab import *
from mpl_toolkits.basemap import Basemap
from polar_projection import *
import Nio
from sattools import *
from gmttools import *

# run plot_one_year_mp_on_one_page_masked.py
# python plot_one_year_mp_on_one_page_masked.py

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
#years=['2000']
for year in years:
    print year
    #year='2003'
    metalist= glob.glob(folder+year+'_*')
    d={}
    for fn in metalist:
        print fn
        date=fn[-8:]
        workfile=fn+'/'+date+'_mp_125masked.nc'
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



    #day257=day257[iYmin:iYmax,iXmin:iXmax]
    day249=day249[iYmin:iYmax,iXmin:iXmax]
    day241=day241[iYmin:iYmax,iXmin:iXmax]
    day233=day233[iYmin:iYmax,iXmin:iXmax]
    day225=day225[iYmin:iYmax,iXmin:iXmax]
    day217=day217[iYmin:iYmax,iXmin:iXmax]
    day209=day209[iYmin:iYmax,iXmin:iXmax]
    day201=day201[iYmin:iYmax,iXmin:iXmax]
    day193=day193[iYmin:iYmax,iXmin:iXmax]
    day185=day185[iYmin:iYmax,iXmin:iXmax]
    day177=day177[iYmin:iYmax,iXmin:iXmax]
    day169=day169[iYmin:iYmax,iXmin:iXmax]
    day161=day161[iYmin:iYmax,iXmin:iXmax]
    day153=day153[iYmin:iYmax,iXmin:iXmax]
    day145=day145[iYmin:iYmax,iXmin:iXmax]
    day137=day137[iYmin:iYmax,iXmin:iXmax]
    day129=day129[iYmin:iYmax,iXmin:iXmax]



    #Tmap=Basemap(resolution='h',projection='stere', boundinglat=70., lon_0=lon_0, lat_0=lat_0,     width=plotsize, height=plotsize, rsphere=(6378273.,6356889.))
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
        Tmap.imshow(n, interpolation='nearest',vmin=0.,vmax=0.5)
    
    
    axisNum = 0
    for m in liste2:
        axisNum += 1
        subplot(4, 4, axisNum)
        month,day=JulDay2Date(int(year),int(m))
        #title(str(day)+'.'+str(month)+'.')
        text(1400000.0,5000000.0,str(month)+'.'+str(day)+'.'+str(year), fontsize='large', color='white')

    subplots_adjust(bottom=0.1, right=0.9, top=0.9,left=0.1)
    cax = axes([0.3, 0.06, 0.4, 0.01])
    colorbar(cax=cax, orientation='horizontal',ticks=[0.0,0.1,0.2,0.3,0.4,0.5],extend='max')  #ticks=[0.0,0.08,0.16,0.24,0.32,0.40]
    text(0.3,-4.5,'melt pond fraction',fontsize='large',color='black')
    savefig('/pf/u/u241127/plots/mp_oneyearplt'+str(year)+'_masked.png')


#######################################
###########year2001 #day169 is missing#
#######################################
year='2001'
metalist= glob.glob(folder+year+'_*')
d={}
for fn in metalist:
    print fn
    date=fn[-8:]
    workfile=fn+'/'+date+'_mp_125masked.nc'
    pic1 = NetCDFFile(workfile,'r')
    img = array(pic1.variables['z'][:,:])
    d[date]=[]
    d[date].append(img)


datasets= sort(d.keys())
X = array(pic1.variables['x'][:])
Y = array(pic1.variables['y'][:])


day249=d[datasets[14]][0]
day241=d[datasets[13]][0]
day233=d[datasets[12]][0]
day225=d[datasets[11]][0]
day217=d[datasets[10]][0]
day209=d[datasets[9]][0]
day201=d[datasets[8]][0]
day193=d[datasets[7]][0]
day185=d[datasets[6]][0]
day177=d[datasets[5]][0]
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
day249=day249[iYmin:iYmax,iXmin:iXmax]
day241=day241[iYmin:iYmax,iXmin:iXmax]
day233=day233[iYmin:iYmax,iXmin:iXmax]
day225=day225[iYmin:iYmax,iXmin:iXmax]
day217=day217[iYmin:iYmax,iXmin:iXmax]
day209=day209[iYmin:iYmax,iXmin:iXmax]
day201=day201[iYmin:iYmax,iXmin:iXmax]
day193=day193[iYmin:iYmax,iXmin:iXmax]
day185=day185[iYmin:iYmax,iXmin:iXmax]
day177=day177[iYmin:iYmax,iXmin:iXmax]
#day169=day169[iYmin:iYmax,iXmin:iXmax]
day161=day161[iYmin:iYmax,iXmin:iXmax]
day153=day153[iYmin:iYmax,iXmin:iXmax]
day145=day145[iYmin:iYmax,iXmin:iXmax]
day137=day137[iYmin:iYmax,iXmin:iXmax]
day129=day129[iYmin:iYmax,iXmin:iXmax]



#Tmap=Basemap(resolution='h',projection='stere', boundinglat=70., lon_0=lon_0, lat_0=lat_0, width=plotsize, height=plotsize, rsphere=(6378273.,6356889.))
Tmap=Basemap(resolution='l',projection='stere', boundinglat=70.,llcrnrlon=lonmin, llcrnrlat=latmin, urcrnrlon=lonmax, urcrnrlat=latmax, lon_0=lon_0, lat_0=lat_0)


liste=[day129,day137,day145,day153,day161,day177,day185,day193,day201,
day209,day217,day225,day233,day241,day249]
liste2=['129','137','145','153','161','177','185','193','201','209','217','225','233','241','249']

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
    Tmap.imshow(n, interpolation='nearest',vmin=0.,vmax=0.56)
    
    
axisNum = 0
for m in liste2:
    axisNum += 1
    subplot(4, 4, axisNum)
    month,day=JulDay2Date(int(year),int(m))
    #title(str(day)+'.'+str(month)+'.')
    text(1400000.0,5000000.0,str(month)+'.'+str(day)+'.'+str(year),fontsize='large',color='white')

subplots_adjust(bottom=0.1, right=0.9, top=0.9,left=0.1)
cax = axes([0.3, 0.06, 0.4, 0.01])
colorbar(cax=cax, orientation='horizontal',ticks=[0.0,0.1,0.2,0.3,0.4,0.5],extend='max')
text(0.3,-4.5,'melt pond fraction',fontsize='large',color='black')
savefig('/pf/u/u241127/plots/mp_oneyearplt'+str(year)+'_masked.png')

####day 169 year 2007 and 2011
years=['2007','2011']
#years=['2000']
for year in years:
    print year
    #year='2003'
    metalist= glob.glob(folder+year+'_169')
    d={}
    for fn in metalist:
        print fn
        date=fn[-8:]
        workfile=fn+'/'+date+'_mp_125masked.nc'
        pic1 = NetCDFFile(workfile,'r')
        img = array(pic1.variables['z'][:,:])
        d[date]=img



    datasets= sort(d.keys())
    X = array(pic1.variables['x'][:])
    Y = array(pic1.variables['y'][:])
    day169=d[datasets[0]]
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
    day169=day169[iYmin:iYmax,iXmin:iXmax]

    #Tmap=Basemap(resolution='h',projection='stere', boundinglat=70., lon_0=lon_0, lat_0=lat_0,     width=plotsize, height=plotsize, rsphere=(6378273.,6356889.))
    Tmap=Basemap(resolution='l',projection='stere', boundinglat=70.,llcrnrlon=lonmin, llcrnrlat=latmin, urcrnrlon=lonmax, urcrnrlat=latmax, lon_0=lon_0, lat_0=lat_0)

    
    year=int(year)
    figure(figsize=(7,7))
    Tmap.drawcoastlines()
    Tmap.fillcontinents(color='gray')
    Tmap.drawmeridians(np.arange(0,360,40),labels=[0,0,0,0])  #labels=[1,0,0,0]
    Tmap.drawparallels(np.arange(-90,90,10),labels=[0,0,0,0])  #labels=[0,0,0,1]
    #title('day '+m)
    Tmap.imshow(day169, interpolation='nearest',vmin=0.,vmax=0.5)
    m=169
    month,day=JulDay2Date(int(year),int(m))
    text(2000000.0,5200000.0,str(month)+'.'+str(day)+'.'+str(year), fontsize='28', color='white')
    cax = axes([0.2, 0.06, 0.6, 0.03])
    colorbar(cax=cax,orientation='horizontal',ticks=[0.0,0.1,0.2,0.3,0.4,0.5],extend='max')  #ticks=[0.0,0.08,0.16,0.24,0.32,0.40]
    savefig('/pf/u/u241127/plots/mp_oneyearplt'+str(year)+str(day)+'masked.png')
    
    
    







"""
figure(figsize=(11,15))
#129
subplot(441)
Tmap.drawcoastlines()
Tmap.fillcontinents(color='gray')
Tmap.drawmeridians(np.arange(0,360,40),labels=[0,0,0,0])  #labels=[1,0,0,0]
Tmap.drawparallels(np.arange(-90,90,10),labels=[0,0,0,0])  #labels=[0,0,0,1]
title('day 129')
Tmap.imshow(day129, interpolation='nearest',vmin=0.,vmax=0.4)
#137
subplot(442)
Tmap.drawcoastlines()
Tmap.fillcontinents(color='gray')
Tmap.drawmeridians(np.arange(0,360,40),labels=[0,0,0,0])  #labels=[1,0,0,0]
Tmap.drawparallels(np.arange(-90,90,10),labels=[0,0,0,0])  #labels=[0,0,0,1]
title('day 137')
Tmap.imshow(day137, interpolation='nearest',vmin=0.,vmax=0.4)
#145
subplot(443)
Tmap.drawcoastlines()
Tmap.fillcontinents(color='gray')
Tmap.drawmeridians(np.arange(0,360,40),labels=[0,0,0,0])  #labels=[1,0,0,0]
Tmap.drawparallels(np.arange(-90,90,10),labels=[0,0,0,0])  #labels=[0,0,0,1]
title('day 145')
Tmap.imshow(day145, interpolation='nearest',vmin=0.,vmax=0.4)
#153
subplot(444)
Tmap.drawcoastlines()
Tmap.fillcontinents(color='gray')
Tmap.drawmeridians(np.arange(0,360,40),labels=[0,0,0,0])  #labels=[1,0,0,0]
Tmap.drawparallels(np.arange(-90,90,10),labels=[0,0,0,0])  #labels=[0,0,0,1]
title('day 153')
Tmap.imshow(day153, interpolation='nearest',vmin=0.,vmax=0.4)
#161
subplot(445)
Tmap.drawcoastlines()
Tmap.fillcontinents(color='gray')
Tmap.drawmeridians(np.arange(0,360,40),labels=[0,0,0,0])  #labels=[1,0,0,0]
Tmap.drawparallels(np.arange(-90,90,10),labels=[0,0,0,0])  #labels=[0,0,0,1]
title('day 161')
Tmap.imshow(day161, interpolation='nearest',vmin=0.,vmax=0.4)
#169
subplot(446)
Tmap.drawcoastlines()
Tmap.fillcontinents(color='gray')
Tmap.drawmeridians(np.arange(0,360,40),labels=[0,0,0,0])  #labels=[1,0,0,0]
Tmap.drawparallels(np.arange(-90,90,10),labels=[0,0,0,0])  #labels=[0,0,0,1]
title('day 169')
Tmap.imshow(day169, interpolation='nearest',vmin=0.,vmax=0.4)
#177
subplot(447)
Tmap.drawcoastlines()
Tmap.fillcontinents(color='gray')
Tmap.drawmeridians(np.arange(0,360,40),labels=[0,0,0,0])  #labels=[1,0,0,0]
Tmap.drawparallels(np.arange(-90,90,10),labels=[0,0,0,0])  #labels=[0,0,0,1]
title('day 177')
Tmap.imshow(day177, interpolation='nearest',vmin=0.,vmax=0.4)
#185
subplot(448)
Tmap.drawcoastlines()
Tmap.fillcontinents(color='gray')
Tmap.drawmeridians(np.arange(0,360,40),labels=[0,0,0,0])  #labels=[1,0,0,0]
Tmap.drawparallels(np.arange(-90,90,10),labels=[0,0,0,0])  #labels=[0,0,0,1]
title('day 185')
Tmap.imshow(day185, interpolation='nearest',vmin=0.,vmax=0.4)
#193
subplot(449)
Tmap.drawcoastlines()
Tmap.fillcontinents(color='gray')
Tmap.drawmeridians(np.arange(0,360,40),labels=[0,0,0,0])  #labels=[1,0,0,0]
Tmap.drawparallels(np.arange(-90,90,10),labels=[0,0,0,0])  #labels=[0,0,0,1]
title('day 193')
Tmap.imshow(day193, interpolation='nearest',vmin=0.,vmax=0.4)
#201
subplot(4,4,10)
Tmap.drawcoastlines()
Tmap.fillcontinents(color='gray')
Tmap.drawmeridians(np.arange(0,360,40),labels=[0,0,0,0])  #labels=[1,0,0,0]
Tmap.drawparallels(np.arange(-90,90,10),labels=[0,0,0,0])  #labels=[0,0,0,1]
title('day 201')
Tmap.imshow(day201, interpolation='nearest',vmin=0.,vmax=0.4)
#209
subplot(4,4,11)
Tmap.drawcoastlines()
Tmap.fillcontinents(color='gray')
Tmap.drawmeridians(np.arange(0,360,40),labels=[0,0,0,0])  #labels=[1,0,0,0]
Tmap.drawparallels(np.arange(-90,90,10),labels=[0,0,0,0])  #labels=[0,0,0,1]
title('day 209')
Tmap.imshow(day209, interpolation='nearest',vmin=0.,vmax=0.4)
#217
subplot(4,4,12)
Tmap.drawcoastlines()
Tmap.fillcontinents(color='gray')
Tmap.drawmeridians(np.arange(0,360,40),labels=[0,0,0,0])  #labels=[1,0,0,0]
Tmap.drawparallels(np.arange(-90,90,10),labels=[0,0,0,0])  #labels=[0,0,0,1]
title('day 217')
Tmap.imshow(day217, interpolation='nearest',vmin=0.,vmax=0.4)
#225
subplot(4,4,13)
Tmap.drawcoastlines()
Tmap.fillcontinents(color='gray')
Tmap.drawmeridians(np.arange(0,360,40),labels=[0,0,0,0])  #labels=[1,0,0,0]
Tmap.drawparallels(np.arange(-90,90,10),labels=[0,0,0,0])  #labels=[0,0,0,1]
title('day 225')
Tmap.imshow(day225, interpolation='nearest',vmin=0.,vmax=0.4)
#233
subplot(4,4,14)
Tmap.drawcoastlines()
Tmap.fillcontinents(color='gray')
Tmap.drawmeridians(np.arange(0,360,40),labels=[0,0,0,0])  #labels=[1,0,0,0]
Tmap.drawparallels(np.arange(-90,90,10),labels=[0,0,0,0])  #labels=[0,0,0,1]
title('day 233')
Tmap.imshow(day233, interpolation='nearest',vmin=0.,vmax=0.4)
#241
subplot(4,4,15)
Tmap.drawcoastlines()
Tmap.fillcontinents(color='gray')
Tmap.drawmeridians(np.arange(0,360,40),labels=[0,0,0,0])  #labels=[1,0,0,0]
Tmap.drawparallels(np.arange(-90,90,10),labels=[0,0,0,0])  #labels=[0,0,0,1]
title('day 241')
Tmap.imshow(day241, interpolation='nearest',vmin=0.,vmax=0.4)
#249
subplot(4,4,16)
Tmap.drawcoastlines()
Tmap.fillcontinents(color='gray')
Tmap.drawmeridians(np.arange(0,360,40),labels=[0,0,0,0])  #labels=[1,0,0,0]
Tmap.drawparallels(np.arange(-90,90,10),labels=[0,0,0,0])  #labels=[0,0,0,1]
title('day 249')
Tmap.imshow(day249, interpolation='nearest',vmin=0.,vmax=0.4)



subplots_adjust(bottom=0.1, right=0.9, top=0.9,left=0.1)
cax = axes([0.11, 0.05, 0.78, 0.02])
colorbar(cax=cax, orientation='horizontal')
#savefig('/pf/u/u241127/plots/mp_oneyearplt2008.png')



"""






