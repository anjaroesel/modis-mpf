# run read_MOD09GA_2grd_snow.py
# python read_MOD09GA_2grd_snow.py
from pylab import *
from scipy import *
import pipes,struct,os,glob,fnmatch
import os
from osgeo.gdalconst import *
from osgeo import gdal 
from gmttools import *
import scipy.ndimage as nd
from sattools import *

def read_dataset(folder): 
    files=glob.glob(folder+'/*.tif') 
    dataset={}
    for fn in files:   
        band=fn[-28:-15]
        dataset[band] = gdal.Open(fn)
    return dataset

def read_MOD(dataset): 
    print('read bands...')
    B={}
    B={'1':array(dataset['_sur_refl_b01'].ReadAsArray(),dtype=float32),\
       '2':array(dataset['_sur_refl_b02'].ReadAsArray(),dtype=float32),\
       '3':array(dataset['_sur_refl_b03'].ReadAsArray(),dtype=float32),\
       '4':array(dataset['_sur_refl_b04'].ReadAsArray(),dtype=float32)}
    band=B.keys
    #scaling of reflectances
    for i in band():
        print('scale band '+i)
        B[i]=B[i]*0.0001
        B[i][B[i]==-28672]=NaN
        B[i][B[i]>=1.0]=nan
        B[i][B[i]<=0.0]=nan
    bands={}
    for i in band():
        bands[i]=B[i]
    return bands

def read_MOD_FLAG(dataset): 
    print('read flags...')
    B={}
    B={'L':array(dataset['_2D_state_1km'].ReadAsArray())}
    band=B.keys
    bands={}
    for i in band():
        bands[i]=B[i]
    return bands




#run on snow
folder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/MOD09/daily/processed/'
metalist = glob.glob(folder+'A2011*/')
print metalist




for fn in metalist:
    foldername1=fn[-16:]
    folderyear=foldername1[1:5]
    folderday=foldername1[5:8]
    foldertile=foldername1[9:15]
    foldername=folder+folderyear+'_'+folderday+'_'+foldertile
    if not os.path.exists(foldername):
        print('grd doesnot exist - start processing of '+fn)
        year=fn[-15:-11]
        julday=fn[-11:-8]
        tile=fn[-7:-1]
        directory=year+'_'+julday+'_'+tile+'/'
        dataset = read_dataset(fn)
        MOD=read_MOD(dataset)
        B1,B2,B3,B4=MOD['1'],MOD['2'],MOD['3'],MOD['4']
        B1[isnan(B1)==True]=0
        B2[isnan(B2)==True]=0
        B3[isnan(B3)==True]=0
        B4[isnan(B4)==True]=0
        FL=read_MOD_FLAG(dataset)
        bflag=FL['L']
        cloud_mask=[]
        land_mask=[]
        print('convert cloud and land flags....')
        for n in bflag.flatten():
            flag=binary_repr(n,width=16)
            #cloudmask #flags von hinten zaehlen!!!!-->little endian
            if flag[15]=='0' and flag[14]=='0' and flag[13]=='0':   
                f=0 #clear
            else:
                f=1 #cloudy    if flag[0]=='1':
            cloud_mask.append(f) #1->clear, 0->cloudy
            #landmask
            if flag[12]=='0' and flag[11]=='0' and flag[10]=='1' or flag[12]=='0' and flag[11]=='1' and flag[10]=='0' or flag[12]=='1' and flag[11]=='0' and flag[10]=='1' or flag[12]=='1' and flag[11]=='0' and flag[10]=='0':
                l=1 #land
            else:
                l=0 #sea
            land_mask.append(l) #1->sea, 0->land
   
        clouds=array(cloud_mask).reshape(bflag.shape)
        land=array(land_mask).reshape(bflag.shape)
        maske=clouds+land 
        maske[maske==2]=1
        maske=abs(maske-1)  
        m=nd.interpolation.zoom(maske,2)
        #dtype muss float32 sein!!!!!!!!!:
        m=array(m,dtype=float32)
        m[m>1.]=0.
        m[m<0.]=0.
        if B1.shape[0]==m.shape[0] and B1.shape[1]==m.shape[1]:
            print('shape fits')
        if B1.shape[0]>m.shape[0] and B1.shape[1]>m.shape[1]:
            print('B1>m')
            B1=B1[:m.shape[0]-B1.shape[0],:m.shape[1]-B1.shape[1]]
            B2=B2[:m.shape[0]-B2.shape[0],:m.shape[1]-B2.shape[1]]
            B3=B3[:m.shape[0]-B3.shape[0],:m.shape[1]-B3.shape[1]]
            B4=B4[:m.shape[0]-B4.shape[0],:m.shape[1]-B4.shape[1]]
        if B1.shape[0]>m.shape[0] and B1.shape[1]==m.shape[1]:
            print('only 1st B1>m')
            B1=B1[:m.shape[0]-B1.shape[0],:]
            B2=B2[:m.shape[0]-B2.shape[0],:]
            B3=B3[:m.shape[0]-B3.shape[0],:]
            B4=B4[:m.shape[0]-B4.shape[0],:]
        if B1.shape[0]==m.shape[0] and B1.shape[1]>m.shape[1]:
            print('only 2nd B1>m')
            B1=B1[:,:m.shape[1]-B1.shape[1]]
            B2=B2[:,:m.shape[1]-B2.shape[1]]
            B3=B3[:,:m.shape[1]-B3.shape[1]]
            B4=B4[:,:m.shape[1]-B4.shape[1]]
        if B1.shape[0]<m.shape[0] and B1.shape[1]<m.shape[1]:
            print('B1<m')
            m=m[:B1.shape[0]-m.shape[0],:B1.shape[1]-m.shape[1]]
        if B1.shape[0]<m.shape[0] and B1.shape[1]==m.shape[1]:
            print('only 1st B1<m')
            m=m[:B1.shape[0]-m.shape[0],:]
        if B1.shape[0]==m.shape[0] and B1.shape[1]<m.shape[1]:
            print('only 2nd B1<m')
            m=m[:,:B1.shape[1]-m.shape[1]]
	if B1.shape[0]>m.shape[0] and B1.shape[1]<m.shape[1]:
            print('first B1>m, second B1<m')
            m=m[:,:B1.shape[1]-m.shape[1]]
	    B1=B1[:m.shape[0]-B1.shape[0],:]
            B2=B2[:m.shape[0]-B2.shape[0],:]
            B3=B3[:m.shape[0]-B3.shape[0],:]
            B4=B4[:m.shape[0]-B4.shape[0],:]
	
	
	if B1.shape[0]<m.shape[0] and B1.shape[1]>m.shape[1]:
            print('first B1<m, second B1>m')
	    m=m[:B1.shape[0]-m.shape[0],:]
            B1=B1[:,:m.shape[1]-B1.shape[1]]
            B2=B2[:,:m.shape[1]-B2.shape[1]]
	    B3=B3[:,:m.shape[1]-B3.shape[1]]
	    B4=B4[:,:m.shape[1]-B4.shape[1]]


        
        ###
        print 'masking'
        B1=B1*m
        B2=B2*m
        B3=B3*m
        B4=B4*m  

        #alle Nullen muessen zu nans werden:
        print '0 zu nans'
        B1[B1==0.]=nan
        B2[B2==0.]=nan
        B3[B3==0.]=nan
        B4[B4==0.]=nan
        """
        figure()
        subplot(221)
        title('cloudmask')
        imshow(clouds)
        colorbar()
        subplot(222)
        title('landmask')
        imshow(land)
        colorbar()
        subplot(223)
        title('total mask')
        imshow(m)
        colorbar()
        show()
        subplot(224)
        title('band with mask')
        imshow(B1)
        colorbar()
        
        """



        #Koordinaten und Pixelgroesse zum spaeteren gridden werden gelesen
        print 'get coordinates'
        geotransform=dataset['_sur_refl_b01'].GetGeoTransform()
        [X, Y]= meshgrid(arange(B1.shape[1],dtype=float32), arange(B1.shape[0],dtype=float32))
        X_polar = (geotransform[0]+geotransform[1]*X+geotransform[2]*Y)/1000#[:,::-1]
        Y_polar = (geotransform[3]+geotransform[4]*X+geotransform[5]*Y)/1000#[:,::-1]
        print '...done'

        #creating filnames and dirs:
        #folder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/MOD09/MOD09GA/'
        os.system('mkdir -p '+folder+directory)
        workfile=folder+'processed/'+directory+year+'_'+julday+'_'+tile
        plotfile=workfile
        out_file1=workfile+'_1.nc'
        out_file2=workfile+'_2.nc'
        out_file3=workfile+'_3.nc'
        out_file4=workfile+'_4.nc'
        tmpfile=workfile+'.gz'
        mapfile=plotfile+'.ps'
        ctabfile=plotfile+'.ctab'
        print(' ')
        print('_____GRIDDING MOD______')
        G1={}
        cs1=0.5#geotransform[1]/1000     #0.5 #500m
        x0,y0,x2,y2=X_polar.max(),Y_polar.max(),X_polar.min(),Y_polar.min()
        grid_par_neu(1,x0,y0,x2,y2,cs1,G1)
        G1=makecpt(ctabfile,'gray',0.0,0.9,0.2,G1)
        #check if grdfiles exist - if not process them
        grdfiles=glob.glob(workfile+'*.nc')
        filelist1=[out_file1,out_file2,out_file3,out_file4]
        print(filelist1)
        for i in range(len(filelist1)):
            if not os.path.exists(filelist1[i]):  
                print('nc gibts noch nicht')
                nearneighbor(X_polar,Y_polar,B1,out_file1,G1)
                nearneighbor(X_polar,Y_polar,B2,out_file2,G1)
                nearneighbor(X_polar,Y_polar,B3,out_file3,G1)
                nearneighbor(X_polar,Y_polar,B4,out_file4,G1)
        #prepare headline:
        titlestr='MOD09 - surface refl: '+year+' '+julday+' tile: '+tile
        pos_x,pos_y=str(x2), str(y0+50)
    
        cmd('gmtset ELLIPSOID WGS-84 OBLIQUE_ANOTATION 1 PAPER_MEDIA A4+ PAGE_ORIENTATION portrait PSIMAGE_FORMAT hex PLOT_DEGREE_FORMAT ddd:mm:ss COLOR_NAN 255/255/255')
        cmd('grdimage '+out_file1+' '+opts(G1,['Rx','Jx','Bx','C'])+' -Y3.5 -K -V> '+mapfile)
        cmd('pscoast '+opts(G1,['Rll','Bll','coast'])+' -O -V -K>> '+mapfile)
        cmd('echo '+pos_x+' '+pos_y+' 18 0 0 0 "'+titlestr+'" | pstext '+opts(G1,['Rx','Jx','N'])+' -V -O >> '+mapfile)
        #cmd('gv '+mapfile+'&')
        
    else:
        print(fn+'was already processed - skip processing')
     
####clean memory
    #del X,Y,X_polar,Y_polar,dataset,MOD,B1,B2,B3,B4,FL,bflag,cloud_mask,land_mask,clouds,land,maske,m
    #print('Schleife zuende')   


