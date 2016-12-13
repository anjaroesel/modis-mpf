#! /usr/bin/env python
#sattools.py


from pylab import *
from scipy import *
from pyhdf.SD import *
from scipy import ndimage
from osgeo.gdalconst import *
from osgeo import gdal 
import pipes,struct,os,glob,fnmatch
from scipy.interpolate import *
import time as tm
from mpl_toolkits.basemap import NetCDFFile
from polar_projection import *

def modis_granules(localdir):
    '''returns available modis hdf files and their available data'''
    sats=['Aqua','Terra']
    md={'Aqua':'MYD','Terra':'MOD'}
    granules={}
    for sat in sats:
        mtypes=[md[sat]+'02HKM',md[sat]+'03',md[sat]+'02QKM',md[sat]+'021KM',md[sat]+'01SS']
        #,md[sat]+'04_L2',md[sat]+'05_L2',md[sat]+'06_L2',md[sat]+'07_L2']
        files=glob.glob(localdir+'*.*.*.*.*.hdf')
        files=glob.glob(localdir+'*.hdf')
        mtdic={}
        for mt in mtypes:
            pattern='*'+mt+'*'
            ffile=fnmatch.filter(files,pattern)
            mtdic[mt]=ffile
        for f1 in mtdic[md[sat]+'03']:
            basefile=os.path.basename(f1)
            blist=basefile.split('.')
            gran=blist[1]+'.'+blist[2]
            granules[gran]={'geo': basefile}
            for f2 in mtdic[md[sat]+'021KM']:
                basefile2=os.path.basename(f2)
                blist2=basefile2.split('.')
                gran2=blist2[1]+'.'+blist2[2]
                if gran==gran2:
                    granules[gran]['1km']=basefile2
            for f2 in mtdic[md[sat]+'02HKM']:
                basefile2=os.path.basename(f2)
                blist2=basefile2.split('.')
                gran2=blist2[1]+'.'+blist2[2]
                if gran==gran2:
                    granules[gran]['500m']=basefile2
            for f2 in mtdic[md[sat]+'02QKM']:
                basefile2=os.path.basename(f2)
                blist2=basefile2.split('.')
                gran2=blist2[1]+'.'+blist2[2]
                if gran==gran2:
                    granules[gran]['250m']=basefile2
            for f2 in mtdic[md[sat]+'01SS']:
                basefile2=os.path.basename(f2)
                blist2=basefile2.split('.')
                gran2=blist2[1]+'.'+blist2[2]
                if gran==gran2:
                    granules[gran]['01SS']=basefile2
    return granules

def read_band(fileObj,key,b1,b2,d):
    '''reads MODIS data,scale them with the Solar Zenith and returns one band'''
    dsObj=fileObj.select(key)
    B=array(dsObj.get(),float32) 
    attr=dsObj.attributes()
    scale=attr['reflectance_scales'][b1]
    offset=attr['reflectance_offsets'][b1]
    dn1=B[b1,:,:]
    dn2=B[b2,:,:]
    bad_ind1=(dn1>=32767).nonzero() # Bit16 zeigt fehlerhafte Messungen an
    bad_ind2=(dn2>=32767).nonzero()
    band1=(scale*(dn1-offset))/cos(d['SolarZenith'])
    band2=(scale*(dn2-offset))/cos(d['SolarZenith'])
    band1[bad_ind1]=0
    band2[bad_ind2]=0
    return band1,band2

def modis_read(folder='20070220/',resolution='1km',band1=0, band2=2):
    d={}
    localdir=folder
    resfactor={'1km':1,'500m':2,'250m':4}
    rf=resfactor[resolution]
    dataset={'1km':'EV_1KM_RefSB','500m':'EV_500_RefSB','250m':'EV_250_RefSB'}
    #Find Modis Files and recognize available Data
    granules=modis_granules(localdir)
    #Chose one File
    k=granules.keys()[0]
    #Read SolarZenith, Latitude and Longitude
    filename_1km=granules[k]['1km']
    print 'read '+localdir+filename_1km
    f_1km=SD(localdir+filename_1km)
    for key in ['SolarZenith','Latitude','Longitude']:
        d[key]=(array(f_1km.select(key).get(),float32))
        if key=='SolarZenith':
            d[key]=d[key]/100.0/360*2*pi # Umrechnung in rad
        d[key]=ndimage.interpolation.zoom(d[key],5*rf,order=1)[:,0:1354*rf]
        
    #Read Chanel
    filename=granules[k][resolution]
    print 'read '+localdir+filename
    f=SD(localdir+filename)
    d['img']=read_band(f,dataset[resolution],band1,band2,d)#d.kd.
    print 'band index '+str(band1) + str(d['img'][0].shape)
    lat=d['Latitude'][:-4*rf,:-4*rf]
    lon=d['Longitude'][:-4*rf,:-4*rf]
    img1=d['img'][0][:-4*rf,:-4*rf]
    img2=d['img'][1][:-4*rf,:-4*rf]
    return lon,lat,img1,img2,k

def select_band(satdata):
    img=satdata[2]
    #print img
    return img

def sonnendistanz(julday):
    day=array([1.,15.,32.,46,60,74,91,106,121,135,152,166,182,196,213,227,242,258,274,288,305,319,335,349,365])
    dist=array  ([0.9832,0.9836,0.9853,0.9878,0.9909,0.9945,0.9993,1.0033,1.0076,1.0109,1.0140,1.0158,1.0167,1.0165,1.0149,1.0128,1.0092,1.0057,1.0011,0.9972,0.9925,0.9892,0.9860,0.9843,0.9833])
    fit=splrep(day,dist)
    dist_sun=splev(julday,fit)
    return dist_sun

def read_dataset(filename): #read in all bands and calculates at sensor's radiance
    files=glob.glob(filename+'*_warp.TIF')
    print files
    dataset={}
    for fn in files:
        band=fn[-11:-10]
        dataset[band] = gdal.Open(fn)
    return dataset

def read_landsat(dataset,satellite,jahr,JulTag,elev_angle,L_min,L_max): 
    #read in all bands and calculates radiance
    #Baender werden eingelsesen
    B={}
    B={'1':array(dataset['1'].ReadAsArray(),dtype=float32),\
       '2':array(dataset['2'].ReadAsArray(),dtype=float32),\
       '3':array(dataset['3'].ReadAsArray(),dtype=float32),\
       '4':array(dataset['4'].ReadAsArray(),dtype=float32),\
       '5':array(dataset['5'].ReadAsArray(),dtype=float32)}#,\
       #'7':array(dataset['7'].ReadAsArray(),dtype=float32)}
    #Radianzen berechnen
    B_rad={}
    band=B.keys
    for i in band():
        #for LS7 - Before July 1, 2000 use this:
        if (satellite=='Landsat7'):
            print 'calc radiance for band '+i+' - Landsat 7'
            #L_min={'1':-6.2,'2':-6.4,'3':-5.0,'4':-5.1,'5':-1.}
            #L_max={'1':293.7,'2':300.9,'3':234.4,'4':241.1,'5':31.06}
            QcalMax = 255.0
            QcalMin = 1.0
            B_rad[i]=((L_max[i]-L_min[i])/(QcalMax-QcalMin))*(B[i]-QcalMin)+L_min[i]
        if (satellite=='Landsat5'):
            print 'calc radiance for band '+i+' - Landsat 5'
            #calibration coefi nach Chander,2003
            #L_min={'1':-1.52,'2':-2.84,'3':-1.17,'4':-1.51,'5':-0.37}
            #L_max={'1':193,'2':365,'3':264,'4':221,'5':30.2}
            QcalMax = 255.0
            QcalMin = 1.0
            B_rad[i]=((L_max[i]-L_min[i])/(QcalMax))*B[i]+L_min[i]
    #return B_rad
    #Reflektanzen berechnen
    B_refl={}
    band=B.keys
    for i in band():
        sun=sonnendistanz(JulTag)
        #print sun
        if (satellite=='Landsat7'):  
            print 'calc reflectance for band '+i+' - Landsat 7'      
            ESUN={'1':1969,'2':1840,'3':1551,'4':1044,'5':225.7,'7':82.07}
        if (satellite=='Landsat5'):
            print 'calc reflectance for band '+i+' - Landsat 5'
            ESUN={'1':1957,'2':1826,'3':1554,'4':1036,'5':215.0,'7':80.67}
        B_refl[i]=pi*B_rad[i]*sun**2/ESUN[i]/cos((90-elev_angle)*pi/180)
    #return B_refl
    #Daten einlesen #(filename,Y1,Y2,X1,X2): #(I=Y,J=X)
    print 'Load data'
    band1=array(B_refl.get('1'),dtype=float32)
    band2=array(B_refl.get('2'),dtype=float32)
    band3=array(B_refl.get('3'),dtype=float32)
    band4=array(B_refl.get('4'),dtype=float32)
    band5=array(B_refl.get('5'),dtype=float32)
    #band7=array(B_refl.get('7'),dtype=float32)
    print '...done'
    print 'image shape: '+ str(band1.shape)
    return band1,band2,band3,band4,band5#,band7

def warping(datadir):
    files=glob.glob(datadir+'*.TIF')
    for fn in files:
        #if not os.path.exists(datadir+'B10_warp.tif'):
            cmd('gdalwarp -s_srs "+proj=utm  +zone='+UTMzone+'N +datum=WGS84" -t_srs "+proj=stere +lat_ts=70 +lat_0=90 +lon_0=0.0 +k_0=1.0 +x_0=0.0 +y_0=0.0" '+fn+' '+datadir+fn[-7:-4]+'_warp.tif')

def cmd(cmdline):
    print cmdline
    os.system(cmdline)
    return

def extractXY(b1,b2,y1,y2,x1,x2):
    X=b1[y1:y2,x1:x2].flatten()
    Y=b2[y1:y2,x1:x2].flatten()
    return X,Y

def plotsat(img,ueberschrift='',cmap=''):
    figure()
    imshow(img,cmap,origin='lower',interpolation='nearest')
    title(ueberschrift)
    colorbar()

def plotimg(mymap,x,y,img,src):
    figure(figsize=[6,6])
    mymap.fillcontinents(color='green')
    mymap.drawmeridians(arange(0,360,2),linewidth=0.5,dashes=[1,0],labels=[0,0,0,1])
    mymap.drawparallels(arange(60,90,0.5),linewidth=0.5,dashes=[1,0],labels=[0,1,0,0])
    mymap.pcolormesh(x,y,img,vmin=0.,vmax=0.8,cmap=cm.gray)
    datestring='%(jahr)04d.%(monat)02d.%(tag)02d' %  \
                       {'jahr':jahr,'monat':monat,'tag':tag}
    title(src+' '+datestring)

def find_pos(lons,lats,lon_0,lat_0,f=4):
    """find the best fitting position of a coordinate pair in an image
    lon: longitude image
    lat: latitude image
    lon_0, lat_0: coordinates of desired position
    f: scale factor to speed up the calculation. Degrades the accurance of the output position
    returns a tuple x,y, which gives the best position in lon,lat of lon_0,lat_0"""
    best_lat_i=lats[::f,::f]
    best_lat_i[best_lat_i>lat_0]=2*lat_0-best_lat_i[best_lat_i>lat_0]
    best_lon_i=lons[::f,::f]
    best_lon_i[best_lon_i>lon_0]=2*lon_0-best_lon_i[best_lon_i>lon_0]
    best_position=best_lat_i/lat_0.__abs__()+best_lon_i/lon_0.__abs__()
    best_y,best_x=ndimage.measurements.maximum_position(best_position)
    best_y,best_x=f*best_y,f*best_x
    return best_x,best_y


def display_rgb(B3,B4,B5):
    RGB=array([B5,B4,B3]).swapaxes(0,1).swapaxes(1,2)
    imshow(RGB,origin='lower',interpolation='nearest')
    show()
    return

def display_mono(B5,cmap=''):
    #imshow(B5,cmap,origin='lower',interpolation='nearest',vmin=0,vmax=1)
    imshow(B5,cmap,origin='lower')
    show()
    return


def get_parameters(metalist,workfolder):
    for case in metalist:
        File = open(case)
        lat = zeros(4)
        lon = zeros(4)
        ppx = zeros(4)
        ppy = zeros(4)
    #set the date and parameters
        for i in range(173):
            line = File.readline()
            if line[4:20] == 'ACQUISITION_DATE':
                jahr, monat, tag = line[23:34].split('-')
                jahr = int(jahr)
                monat = int(monat)
                tag = int(tag)
            if line[4:17] == 'SUN_ELEVATION':
                elev_angle = float(line[20:29])
            if line[4:15] == 'ZONE_NUMBER':
                UTMzone = int(line[18:20])
            if line[4:17] == 'SPACECRAFT_ID':
                satellite = line[21:29] ### andere Bezeichnung als in sattools.read_landsat
            if line[4:25] == 'PRODUCT_UL_CORNER_LAT' :
                lat[0] = line[28:39]
            if line[4:25] == 'PRODUCT_UR_CORNER_LAT' :
                lat[1] = line[28:39]
            if line[4:25] == 'PRODUCT_LL_CORNER_LAT' :
                lat[2] = line[28:39]
            if line[4:25] == 'PRODUCT_LR_CORNER_LAT' :
                lat[3] = line[28:39]    
            if line[4:25] == 'PRODUCT_UL_CORNER_LON' :
                lon[0] = line[28:39]
            if line[4:25] == 'PRODUCT_UR_CORNER_LON' :
                lon[1] = line[28:39]
            if line[4:25] == 'PRODUCT_LL_CORNER_LON' :
                lon[2] = line[28:39]
            if line[4:25] == 'PRODUCT_LR_CORNER_LON' :
                lon[3] = line[28:39] 

        for k in range(4):
            ppx[k], ppy[k] = mapll(float(lat[k]),float(lon[k]),1)
            x2 = min(ppx)
            x0 = max(ppx)
            y2 = min(ppy)
            y0 = max(ppy)
        
        JulTag = Julday2(jahr,monat,tag)
        landsatpath = case[0:79]
        workfile=workfolder+case[48:79]+'ls_'+str(jahr)+'_'+str(monat)+'_'+str(tag)
        tmpfile=workfile+'.gz'
        mapfile=workfile+'.ps'
        ctabfile=workfile+'.ctab'
    return workfile,tmpfile,mapfile,ctabfile,satellite,landsatpath,JulTag,x2,x0,y2,y0,jahr,monat,tag

def Julday(year,month,day):  #liefert fuer 2009 falsche juldays -Julday2 verwenden!!!
    hr = 12  
    t = tm.mktime((year, month, day, hr, 0, 0, 0, 0, -1))
    JulDay = tm.localtime(t)[7]
    return JulDay



def Julday2(year,month,day):
    #Schaltjahre 1996,2000,2004,2008
    if year==2000 or 1992 or 1988 or 1996 or 2004 or 2008 or 2012 or 1984 or 2016:
        if month==1:
            mm=0
        if month==2:
            mm=31
        if month==3:
            mm=31+29
        if month==4:
            mm=31+29+30
        if month==5:
            mm=31+29+30+31
        if month==6:
            mm=31+29+30+31+30
        if month==7:
            mm=31+29+30+31+30+31
        if month==8:
            mm=31+29+30+31+30+31+30
        if month==9:
            mm=31+29+30+31+30+31+30+31
        if month==10:
            mm=31+29+30+31+30+31+30+31+30
        if month==11:
            mm=31+29+30+31+30+31+30+31+30+31
        if month==12:
            mm=31+29+30+31+30+31+30+31+30+31+30
        julday=mm+day
    else:
        if month==1:
            mm=0
        if month==2:
            mm=31
        if month==3:
            mm=31+28
        if month==4:
            mm=31+28+30
        if month==5:
            mm=31+28+30+31
        if month==6:
            mm=31+28+30+31+30
        if month==7:
            mm=31+28+30+31+30+31
        if month==8:
            mm=31+28+30+31+30+31+30
        if month==9:
            mm=31+28+30+31+30+31+30+31
        if month==10:
            mm=31+28+30+31+30+31+30+31+30
        if month==11:
            mm=31+28+30+31+30+31+30+31+30+31
        if month==12:
            mm=31+28+30+31+30+31+30+31+30+31+30
        julday=mm+day    
    return julday



def JulDay2Date(year,n):
    import datetime
    d = datetime.date(year, 1, 1) + datetime.timedelta(n - 1)
    month=d.month
    day=d.day
    return month, day





def find_cut_area_coord(MY,MX,lon_y,lat_x,f):
    """find the best fitting position of a coordinate pair in an image
    MY: longitude image array
    MX: latitude image
    lon_0, lat_0: coordinates of desired position
    f: +- amount of pixel around the centercoordinate"""
    for index, item in enumerate(MX):
        if item > lat_x and item < lat_x+1:
            ix=index
            ix_value=item
            #print(index,item)
    for index, item in enumerate(MY):
        if item > lon_y and item < lon_y+1:
            iy=index
            iy_value=item
            #print(index,item)
    X0=ix-f
    X2=ix+f
    Y0=iy-f
    Y2=iy+f
    return X0,X2,Y0,Y2
    
