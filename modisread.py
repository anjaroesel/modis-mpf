#modisread.py
#Python Module for reading modis Level1B image(s), geolocations and removing of the bow-tie effect.
#
#written by Johannes Roehrs and Lars Kaleschke, IfM, University of Hamburg

import scipy as sp
import scipy.ndimage as nd
from scipy.interpolate import interp1d
from pyhdf.SD import *
import pipes,struct,os,glob,fnmatch
import pylab as lab

def convert_decimal_to_8bit(x):
    bstr=str(1-(sp.sign(x)+1)/2)
    for n in range(6,-1,-1):
        bstr=bstr+str(abs(x)/(2**n))
        x=abs(x)%(2**n)
    return bstr

ovp={'250m':[2.72813387e-06,-1.48673672e-02,2.06285669e+01],'500m':[5.36825540e-06,-1.46392798e-02,1.02198638e+01],'1km':[  8.60640972e-06,-1.17260231e-02,4.79268727e+00]}

def modisread(folder, resolution='1km', bandlist=[0], calc_overlap='True'):
    """reads modis Level1B image(s), geolocations and removes the bow-tie effect.

    folder:       location of the image files
    resolution:   '1km', '500m' or '250m'
    bandlist:     list of the bands to read,
                  where the bands are the i_th band of the given resolution starting from 0
    calc_overlap: set this to 'True', if the bowtie-overlap pattern
                  shall be calculated for the specific image

    returns lon,lat,image_list
    lon:          array of Longitude coordinates
    lat:          array of Latitude coordinates
    image_list:   List of 2d Arrays which are the images corresponding to bandlist

    written by Johannes Roehrs and Lars Kaleschke, IfM, University of Hamburg

    """
    print ' '
    print '____________Read Modis Data_____________'
    print ' '
    #os.system('rm bowtie.so') # Comment out for rebuilding Fortran module

    if not(os.path.exists('bowtie.so')):
        os.system('f2py --fcompiler=gnu95  -m bowtie -c ./bowtie.f95 ')

    import bowtie #reload(bowtie)

    d={}
    cs_dict={'1km':1000,'500m':500,'250m':250}
    cs=cs_dict[resolution]
    dataset={'1km':'EV_1KM_RefSB','500m':'EV_500_RefSB','250m':'EV_250_RefSB'}
    #Find Modis Files and recognize available Data
    granules=modis_granules(folder)
    #Chose one File
    k=granules.keys()[0]
    #
    # Read SolarZenith, Latitude and Longitude
    #
    filename_1km=granules[k]['1km']
    print 'read '+folder+filename_1km
    f_1km=SD(folder+filename_1km)
    print ' '   
    for key in ['SensorAzimuth','SensorZenith','SolarAzimuth','SolarZenith','Latitude','Longitude']:
        d[key]=(sp.array(f_1km.select(key).get(),sp.float32))
        if key=='SolarZenith':
            #d[key]=d[key]/100.0
            d[key]=d[key]/100.0/360*2*sp.pi
        if key=='SolarAzimuth':
            d[key]=d[key]/100.0   #scalefactor 0.01
        if key=='SensorZenith':
            d[key]=d[key]/100.0   #scalefactor 0.01
        if key=='SensorAzimuth':
            d[key]=d[key]/100.0   #scalefactor 0.01
        print 'interpolating ' +key+' onto a grid...'
        # remove bowtie effect:
        # this is quick and dirty but works quite well in the test region
        # Averaging should not be done with lon/lat values!
        #avfilter=sp.array([[0,1,0],[0,2,0],[0,1,0]])/4.
        #d[key]=nd.filters.convolve(d[key],avfilter) #weighted average between neighbored location
        # resample to image resolution:
        d[key]=nd.interpolation.zoom(d[key],5000/cs,order=1,mode='constant')[:,:-1000/cs]
        print ' done.'

 
    #Read Chanel
    filename=granules[k][resolution]
    print 'load '+folder+filename+'... '
    f=SD(folder+filename)
    print ' done'
    d['img']=[]
    for band in bandlist:
        IMG=read_band(f,dataset[resolution],band,d)#d.kd.
        # get the bowtie overlap pattern
        if calc_overlap=='True':
            p=sp.array(bowtie_polynom(IMG,cs,folder),dtype=float,order='Fortran')
        else:
            p=sp.array(ovp[resolution],dtype=float,order='Fortran')
        # remove bowtie-effect from image
        modis_img=sp.array(IMG.transpose(),dtype=float,order='Fortran')
        bowtie.bowtie(modis_img,p,int(10000/cs),modis_img.shape[0],modis_img.shape[1])
        IMG=sp.array(modis_img.transpose(),order='C')
        d['img'].append(IMG[:,1000/cs:]) #remove the left-hand-side edge pixel, so that geolocations fit
    #Make the image fit to the geolocations:
    d['Latitude']=d['Latitude'][:,:-1000/cs] #remove the right-hand-side edge pixels
    d['Longitude']=d['Longitude'][:,:-1000/cs]
    d['SolarZenith']=d['SolarZenith'][:,:-1000/cs]
    d['SolarAzimuth']=d['SolarAzimuth'][:,:-1000/cs]
    d['SensorAzimuth']=d['SensorAzimuth'][:,:-1000/cs]
    d['SensorZenith']=d['SensorZenith'][:,:-1000/cs]
    print 'image shape: '+ str(d['img'][0].shape)
    print ' '
    return d['Longitude'],d['Latitude'],d['img'],d['SolarZenith'],d['SolarAzimuth'],d['SensorAzimuth'],d['SensorZenith']


def modisread_day_of_year(folder):
    print ' '
    print '____________Read Modis Day_of_Year_____________'
    print ' '
    d={}
    cs_clouds=1000
    granules=modis_granules(folder)
    #Chose one File
    k=granules.keys()[0]
    filename_doy=granules[k]['day_of_year']
    print 'read day_of_year: '+folder+filename_doy
    f_doy=SD(folder+filename_doy)
    print ' '
    for key in ['day_of_year','Latitude','Longitude']:
        d[key]=sp.array(f_doy.select(key).get(),sp.int16)
        if key=='day_of_year':
            day_of_year=[]
            for n in d[key][0].flatten():
                #flag=convert_decimal_to_8bit(n)
                flag=sp.binary_repr(n,width=16)
                day_of_year.append(int(flag)) #0->cloudy, 1->clear
            doy=sp.array(day_of_year).reshape(d[key].shape[1:3])
            d[key]=doy

    print ' done.'
    d['day_of_year']=d['day_of_year'][:,:-1000/cs_clouds]
    d['Latitude']=d['Latitude'][:,:-1000/cs_clouds] #remove the right-hand-side edge pixels
    d['Longitude']=d['Longitude'][:,:-1000/cs_clouds]
    print 'image shape: '+ str(d['day_of_year'].shape)
    print ' '
    return d['day_of_year'],d['Longitude'],d['Latitude']

def modisread_clouds(folder):
    print ' '
    print '____________Read Modis Cloudmask - kann etwas dauern_____________'
    print ' '
    d={}
    cs_clouds=1000
    granules=modis_granules(folder)
    #Chose one File
    k=granules.keys()[0]
    #
    # Read Cloudmask
    #
    filename_clouds=granules[k]['clouds']
    print 'read clouds: '+folder+filename_clouds
    f_clouds=SD(folder+filename_clouds)
    print ' '
    for key in ['Cloud_Mask','Latitude','Longitude']:
        d[key]=sp.array(f_clouds.select(key).get(),sp.int8)
        if key=='Cloud_Mask':
            cloud_mask=[]
            land_mask=[]
            for n in d[key][0].flatten():
                #flag=convert_decimal_to_8bit(n)
                flag=sp.binary_repr(n,width=8)
                cloud_mask.append(int(flag[5])) #0->cloudy, 1->clear
                land_mask.append(int(flag[0])) #0->water 1->coast,desert,land
            clouds=sp.array(cloud_mask,sp.int8).reshape(d[key].shape[1:3])
            land=abs(sp.array(land_mask,sp.int8).reshape(d[key].shape[1:3])-1)#1->water 0->coast,desert,land
            d[key]=land*clouds

    print ' done.'
    d['Cloud_Mask']=d['Cloud_Mask'][:,:-1000/cs_clouds]
    d['Latitude']=d['Latitude'][:,:-1000/cs_clouds] #remove the right-hand-side edge pixels
    d['Longitude']=d['Longitude'][:,:-1000/cs_clouds]
    print 'image shape: '+ str(d['Cloud_Mask'].shape)
    print ' '
    return d['Cloud_Mask'],d['Longitude'],d['Latitude']





def modis_granules(localdir):
    '''returns available modis hdf files and their available data'''
    sats=['Terra','Aqua']
    md={'Terra':'MOD','Aqua':'MYD'}
    granules={}
    for sat in sats:
        mtypes=[md[sat]+'02HKM',md[sat]+'03',md[sat]+'02QKM',md[sat]+'021KM',md[sat]+'01SS',md[sat]+'35_L2']
        #,md[sat]+'05_L2',md[sat]+'06_L2',md[sat]+'07_L2']
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
            for f3 in mtdic[md[sat]+'35_L2']:
                basefile3=os.path.basename(f3)
                blist3=basefile3.split('.')
                gran3=blist3[1]+'.'+blist3[2]
                if gran==gran3:
                    granules[gran]['clouds']=basefile3
    return granules

def read_band(fileObj,key,b,d):
    '''reads MODIS data,scale them with the Solar Zenith and returns one band'''
    print 'read image '+str(b)+' from dataset '+str(key) +'... '
    dsObj=fileObj.select(key)
    B=sp.array(dsObj.get(),sp.float32) # float32
    attr=dsObj.attributes()
    scale=attr['reflectance_scales'][b]
    offset=attr['reflectance_offsets'][b]
    dn=B[b,:,:]
    bad_ind=(dn>=32767).nonzero() # Bit16 zeigt fehlerhafte Messungen an
    band=(scale*(dn-offset))/sp.cos(d['SolarZenith'])
    band[bad_ind]=0
    print ' done'
    return band

def bowtie_py(modis_img,p,cs):
    print 'removie bowtie effect from image... '
    stripwidth=10000/cs
    #Loop over every x coordinate of the image
    for x in sp.arange(modis_img.shape[1]):
        #Loop over every sanning strip
        overlap=sp.polyval(p,x).round() #get the overlap from the polynom
        if overlap > 0:
            for y in sp.arange(stripwidth,modis_img.shape[0],stripwidth):
                #cut out the current part:
                strippart=modis_img[y:y+stripwidth,x]
                #delete the upper and lower few pixels of the strippart:
                strippart=strippart[int(overlap/2.):-int(round(overlap/2.))]
                #Interpolat to stipwidth length
                f=interp1d(sp.arange(0,len(strippart)),strippart)
                strippart=f(sp.linspace(0,len(strippart)-1,stripwidth))
                #replace the current sick part in the image by the new healthy one
                modis_img[y:y+stripwidth,x]=strippart
    print 'done'
    return modis_img

def bowtie_polynom(modis_img,cs,folder):
    print 'Determine overlap pattern... '
    sw=10000/cs #stripwidth
    overlaplist=[]#define list to store number of overlapped lines
    #devide in parts with a width of 40 pixel
    for i in sp.arange(0,modis_img.shape[1]-40,40):
        part=modis_img[:,i:i+39]
        #search in every scanning strip
        samples=[]
        for j in sp.arange(sw-2,part.shape[0]-sw,sw):
            target=part[j-1:j+1,:] #cut out a target, which overlapped counter-part shall be found
            searchwindow=part[j+2:j+sw+2] #,: cut out the window, where the overlapped counter part might be located
            #start the search
            c=[] #calculate correlation coefficients of every given offset from 3 to 11
            for offset in sp.arange(3,sw/2+1):
                imgpart=searchwindow[offset-3:offset-1] #,: cut out image, which has to be compared with the target
                c.append(sp.corrcoef(imgpart.flatten(),target.flatten())[0,1])#calculate correlatoin coefficient
            c=sp.array(c)
            overl=sp.ndimage.measurements.maximum_position(c)[0]+3 #find the overlap with the highes correlation coefficient
            samples.append([overl,c.max()]) #attach overlap and correlation coefficient to the sample list
        samples=sp.array(samples)
        #print i, samples[:,1].mean()
        if samples[:,1].mean() > 0.9: #chek the mean correlation coefficient:
            #print('Bowtie Correlation high - removing effect')
            overlaplist.append([i+20,samples[:,0].mean()]) #save result, if correlation coefficient is high
            #print(overlaplist)
            o=sp.array(overlaplist)
            X=o[:,0]
            overlap=o[:,1]
            #Calculate a second order Polynom to describe the overlap
            p=sp.polyfit(X,overlap,2)
            #print 'done, Overlap polynom: '+str(p)
        else:
            #print('low Bowtie correlation')
            p = [1.,  1.,  1.]
            #overlaplist.append([i+20,1])
            #os.system('rm -r '+folder)
            #print('scene deleted') 
    return p



def find_pos(lon,lat,lon_0,lat_0,f=4):
    '''find the best fitting position of a coordinate pair in an image
    lon: longitude image
    lat: latitude image
    lon_0, lat_0: coordinates of desired position
    f: scale factor to speed up the calculation. Degrades the accurance of the output position
    returns a tuple x,y, which gives the best position in lon,lat of lon_0,lat_0'''
    best_lat_i=lat[::f,::f]
    best_lat_i[best_lat_i>lat_0]=2*lat_0-best_lat_i[best_lat_i>lat_0]
    best_lon_i=lon[::f,::f]
    best_lon_i[best_lon_i>lon_0]=2*lon_0-best_lon_i[best_lon_i>lon_0]
    best_position=best_lat_i/lat_0.__abs__()+best_lon_i/lon_0.__abs__()
    best_y,best_x=nd.measurements.maximum_position(best_position)
    best_y,best_x=f*best_y,f*best_x
    return best_x,best_y





    

