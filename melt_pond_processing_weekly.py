#! /usr/bin/env python
# run melt_pond_processing_weekly.py
# python melt_pond_processing_weekly.py
from pylab import *
from scipy import *
from osgeo.gdalconst import *
from osgeo import gdal 
from gmt_tools import *
from sattools import *
import sys,pipes,struct,os,glob,fnmatch,pickle

#from mpl_toolkits.basemap import NetCDFFile  #laeuft so auf so snow
from netCDF4 import Dataset as NetCDFFile    # for local runs

from scipy.interpolate import *
from scipy.ndimage import *
import time
from ffnet import ffnet,mlgraph
import networkx
from scipy.interpolate import *
import scipy.io as io
import scipy.optimize as opti
from mpl_toolkits.basemap import Basemap
import scipy.ndimage as nd


#READ FIRST:
#before you start, you need to go through the following steps:

# 1. download MODIS weekly *.hdf files in a specific directory
# 2. cp the files 'execute_reproj_via_shell.sh' and 'modis_hdf2tif_polster_wgs84.sh' in the same directory and adjust pathes in the 'execute_reproj_via_shell.sh' file
# 3. apply 'execute_reproj_via_shell.sh' with the comand sh execute_reproj_via_shell.sh
# 4. Change Directories in this program


#run on snow
####DEFINE FOLDERS####
#satellite data folder
#folder='/scratch/clisap/seaice/MODIS/2012_weekly/'

#workfolder=folder+'processed/'
#tmpfolder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/tmp/'

#Productfolder
#prodfolder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/PRODUCTS/2012_weekly/'

####run local
####DEFINE FOLDERS####
#satellite data folder
folder='/home/u241127/data/modis/weekly/'

workfolder=folder+'processed/'
tmpfolder='/home/u241127/tmp/'
NNfolder='/home/u241127/data/NN/'

#Productfolder
prodfolder='/home/u241127/MELT_PONDS/PRODUCTS/2011_weekly/'
year='2011'

####FUNCTIONS#####
def read_dataset(folder): 
    files=glob.glob(folder+'/*.tif') 
    dataset={}
    for fn in files:   
        band=fn[-26:-13]
        dataset[band] = gdal.Open(fn)
    return dataset

def read_MOD(dataset): 
    print('read bands 3412...')
    B={}
    B={'3':array(dataset['_sur_refl_b01'].ReadAsArray(),dtype=float32),\
       '4':array(dataset['_sur_refl_b02'].ReadAsArray(),dtype=float32),\
       '1':array(dataset['_sur_refl_b03'].ReadAsArray(),dtype=float32),\
       '2':array(dataset['_sur_refl_b04'].ReadAsArray(),dtype=float32)}
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
    B={'L':array(dataset['fl_state_500m'].ReadAsArray())}
    band=B.keys
    bands={}
    for i in band():
        bands[i]=B[i]
    return bands
    
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


def JulDay2Date(year,n):
    import datetime
    d = datetime.date(year, 1, 1) + datetime.timedelta(n - 1)
    month=d.month
    day=d.day
    return month, day




########START
#######read_MODIS2grd


metalist = glob.glob(folder+'A20*/')
print metalist

cmd('rm -r '+folder+'*v03')


for fn in metalist:
    foldername1=fn[-16:]
    folderyear=foldername1[1:5]
    folderday=foldername1[5:8]
    foldertile=foldername1[9:15]
    foldername2=folder+'processed/'+folderyear+'_'+folderday+'_'+foldertile
    if not os.path.exists(foldername2):
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
        #dtype muss float32 sein!!!!!!!!!:
        m=array(maske,dtype=float32)
  
        ###masking
        B1=B1*m
        B2=B2*m
        B3=B3*m
        B4=B4*m  

        #all zeros to nans:
        B1[B1==0.]=nan
        B2[B2==0.]=nan
        B3[B3==0.]=nan
        B4[B4==0.]=nan

	




        #Koordinaten und Pixelgroesse zum spaeteren gridden werden gelesen
        print 'get coordinates'
        geotransform=dataset['_sur_refl_b01'].GetGeoTransform()
        [X, Y]= meshgrid(arange(B1.shape[1],dtype=float32), arange(B1.shape[0],dtype=float32))
        X_polar = (geotransform[0]+geotransform[1]*X+geotransform[2]*Y)/1000#[:,::-1]
        Y_polar = (geotransform[3]+geotransform[4]*X+geotransform[5]*Y)/1000#[:,::-1]
        print '...done'

        #creating filnames and dirs:
        #folder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/MOD09/MOD09GA/'
        os.system('mkdir -p '+folder+'processed/'+directory)
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
                near_neighbor(X_polar,Y_polar,B1,out_file1,G1)
                near_neighbor(X_polar,Y_polar,B2,out_file2,G1)
                near_neighbor(X_polar,Y_polar,B3,out_file3,G1)
                near_neighbor(X_polar,Y_polar,B4,out_file4,G1)
        #prepare headline:
        titlestr='MOD09 - surface refl: '+year+' '+julday+' tile: '+tile
        pos_x,pos_y=str(x2), str(y0+50)
    
        cmd('GMT gmtset ELLIPSOID WGS-84 OBLIQUE_ANOTATION 1 PAPER_MEDIA A4+ PAGE_ORIENTATION portrait PSIMAGE_FORMAT hex PLOT_DEGREE_FORMAT ddd:mm:ss COLOR_NAN 255/255/255')
        cmd('GMT grdimage '+out_file1+' '+opts(G1,['Rx','Jx','Bx','C'])+' -Y3.5 -K -V> '+mapfile)
        cmd('GMT pscoast '+opts(G1,['Rll','Bll','coast'])+' -O -V -K>> '+mapfile)
        cmd('GMT echo '+pos_x+' '+pos_y+' 18 0 0 0 "'+titlestr+'" | pstext '+opts(G1,['Rx','Jx','N'])+' -V -O >> '+mapfile)
        #cmd('gv '+mapfile+'&')
        
    else:
        print(fn+'was already processed - skip processing')
     
####clean memory
    #del X,Y,X_polar,Y_polar,dataset,MOD,B1,B2,B3,B4,FL,bflag,cloud_mask,land_mask,clouds,land,maske,m
    #print('Schleife zuende')   
    
    
  
for fn in metalist:
    foldername1=fn[-16:]
    folderyear=foldername1[1:5]
    folderday=foldername1[5:8]
    foldertile=foldername1[9:15]
    foldername=folder+'processed/'+folderyear+'_'+folderday+'_'+foldertile
    year=fn[-15:-11]
    julday=fn[-11:-8]
    tile=fn[-7:-1]

#####MOSAIC
print('MOSAIKE WERDEN GEBAUT')

metalist_m=glob.glob(workfolder+year+'*/')
print(metalist_m)
daylist=[]

for day in metalist_m:
    d=day[-11:-8]
    daylist.append(d)

days=list(sorted(set(daylist)))

#####
#####
#days=['147', '148', '149', '150', '151', '152', '153', '154', '155', '156', '157', '158']
#days=['249']
print(days)
#############
#BLENDFILE
###############


#for fn in metalist_m:
    #if not os.path.exists(foldername):
        #print('grd doesnot exist - start processing of '+fn)



for b in range(4):
    band=b+1
    print('working on band '+str(band))
    for d in days:
        os.system('mkdir '+workfolder+year+'_'+str(d)+'_mosaic/')
        workfolder_mosaic=workfolder+year+'_'+str(d)+'_mosaic/'
        workfile=workfolder_mosaic+year+'_'+str(d)+'_mosaic'
        blendfile=workfile+str(band)+'.job'
        blendoutfile=workfile+str(band)+'.nc'    
        ctabfile=workfile+str(band)+'.ctab'
        mapfile=workfile+str(band)+'.ps'
        piclist = glob.glob(workfolder+year+'_'+str(d)+'*/*'+str(band)+'.nc') 
        print(piclist)
        mosaic = file(blendfile,'w')
        for fn in piclist:
            print('now working on '+fn)
            file1=fn
            #read Netcdfs and extract Bands
            mb1=NetCDFFile(file1)
            print('read nc') 
            MX=array(mb1.variables['x'][:])
            MY=array(mb1.variables['y'][:])
            xmin = min(MX)
            xmax = max(MX)
            ymin = min(MY)
            ymax = max(MY)
            mosaic.write(fn+' -R'+str(xmin)+'/'+str(xmax)+'/'+str(ymin)+'/'+str(ymax)+' 1'+'\n') 
        mosaic.close()
        #print('######')
        print(file(blendfile).read())	
        #print('######')
        G={}
        cs=0.5 #500m
        region='Arc'
        grid_par_reg(region,cs,G)
        G=makecpt(ctabfile,'gray',0.0,0.9,0.2,G)
        cmd('GMT gmtset OBLIQUE_ANOTATION 1 PAPER_MEDIA A4+ PAGE_ORIENTATION portrait PLOT_DEGREE_FORMAT ddd:mm:ss   COLOR_NAN 255/255/255')
        cmd('GMT grdblend '+blendfile+' -G'+blendoutfile+' '+opts(G,['Rx','I'])+' -V > '+mapfile)
        cmd('GMT grdimage '+blendoutfile+' '+opts(G,['Rx','Jx','Bx','C'])+' -P -K -V>> '+mapfile)
        cmd('GMT psscale -D8c/-0.5c/12c/0.4ch -C'+ctabfile+' -O -B0.1 -V -K >> '+mapfile)
        cmd('GMT pscoast '+opts(G,['Rll','Bll','coast'])+' -G128/128/128 -O -V -K >> '+mapfile)
        cmd('GMT echo -2000 6500 18 0 0 0 "Overview B'+str(band)+' Day '+str(d)+'" | pstext '+opts(G,['Rx','Jx','N'])+' -V -O >> '+mapfile)
        #cmd('gv '+mapfile+'&')
	
#aufraeumen
cmd('mkdir '+workfolder+'tiles')
cmd('mv '+workfolder+' 20*_h* '+workfolder+'tiles')

####Neural Net apply
print('NEURONALES NETZ WIRD AKTIVIERT')

############
##METALIST## 
############
metalist_n = glob.glob(workfolder+'*_mosaic/')
print(metalist_n)


for folder in metalist_n:
    print folder
    jahr=folder[-16:-12]
    julday=folder[-11:-8]
    date=jahr+'_'+julday
    workfile=folder+date+'_mosaic'
    print('work on '+workfile)
    if not os.path.exists(prodfolder+date+'/'+date+'_ow.data'):
        print date+' does not exist'
        file1,file2,file4=workfile+'1.nc',workfile+'2.nc',workfile+'4.nc'
        #read Netcdfs and extract Bands
        mb1,mb2,mb4=NetCDFFile(file1),NetCDFFile(file2),NetCDFFile(file4)
        print('read nc') 
        #MX=array(mb1.variables['x'][:])
        #MY=array(mb1.variables['y'][:])
        MB1=array(mb1.variables['z'][:,:])
        MB2=array(mb2.variables['z'][:,:])
        MB4=array(mb4.variables['z'][:,:])
        titlestr='MODIS mosaic- '+date
        MB1[MB1==0.]=NaN
        MB2[MB2==0.]=NaN
        MB4[MB4==0.]=NaN
        BAND1=MB1.flatten()
        BAND2=MB2.flatten()
        BAND4=MB4.flatten()
        R=array([BAND1,BAND2,BAND4])
        print('...done')
        print 'loading NN'
        #net=pickle.load(open(tmpfolder+'neural_net_TD.data','r'))
        #net=pickle.load(open(tmpfolder+'neural_net_TDneu3_only_nomp.data','r'))
        net=pickle.load(open(tmpfolder+'neural_net_TDneu.data','r'))   
        net=pickle.load(open(NNfolder+'net5neu.data','r'))
        #net=pickle.load(open(tmpfolder+'neural_net_my_sib_can_fram_MP_valuesTSCHUDI','r')) #ergebnisse schlecht
	    ###net=pickle.load(open(tmpfolder+'neural_net_my_sib_can_fram','r'))  #GUTE RESULTS
        #zeitnahme start
        start = time.clock()
        print('...done')
        print 'Applying NN'
        Z_neural=net(R.transpose()).transpose()
        #zeitnahme ende
        stend = time.clock()
        print('dauer: '+str(stend-start))

        mp=Z_neural[0].reshape((MB1.shape))
        snow=Z_neural[1].reshape((MB1.shape))
        ow=Z_neural[2].reshape((MB1.shape))
        print('...done')
        print('pickle data')
        cmd('mkdir '+prodfolder+date)
        pickle.dump(mp,open(prodfolder+date+'/'+date+'_mp.data','w'))
        print('mp pickled')
        pickle.dump(snow,open(prodfolder+date+'/'+date+'_snow.data','w'))
        print('snow pickled')
        pickle.dump(ow,open(prodfolder+date+'/'+date+'_ow.data','w'))
        print('ow pickeld......done')


########regrid
#prodfolder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/PRODUCTS/2012_June_daily/'
#prodfolder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/PRODUCTS/KARA/'
#prodfolder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/PRODUCTS/GREENLAND/'


############
print('REGRID')
metalist = glob.glob(prodfolder+'20*/*mp.data')
print(metalist)

for dataset in metalist:
    print dataset
    #####select date:
    date=dataset[-16:-8]#'2008_169'
    year=int(date[0:4])
    n=int(date[5:])
    print date, n, year
    month, day = JulDay2Date(year,n)
    month=str(month)
    day=str(day)
    year=str(year)
    DATE=day+'.'+month+'.'+year
    print DATE
    ####---snow1:
    workfolder=prodfolder+date+'/'
    tmpfolder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/tmp/'
    #plotfolder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/PRODUCTS/daily/plots/'
    plotfolder=folder+'plots/'
    jahr=date[0:4]
    julday=date[5:8]
    titlestr='MODIS MP - '+DATE+'('+str(n)+')'
    #check this:!!!
    if not os.path.exists(workfolder+date+'_mp_125.nc' or workfolder+date+'_mp_25.nc' or workfolder+date+'_mp_05.nc'):
        print('need to load all objects')
        print('load mp object')
        mp=pickle.load(open(workfolder+date+'_mp.data','r'))
        print('load ice object')
        sn=pickle.load(open(workfolder+date+'_snow.data','r'))
        print('load ow object')
        ow=pickle.load(open(workfolder+date+'_ow.data','r'))
        print('...done')
        print('mpset ('+workfolder+date+'_mp_*.nc) doesnot exist - start processing of '+dataset)
        mp=array(mp,dtype=float32)
        print('mp shape:'+str(mp.shape))
        cs_sm=0.5
        region='Arc'
        XYs=XYgrid(region,cs_sm)
        X=array(XYs[0],dtype=float32)
        Y=array(XYs[1],dtype=float32)
        i=isfinite(mp)
        mpi=mp[i]
	owi=ow[i]
	sni=sn[i]
	Xi=X[i]
        Yi=Y[i]
        FILENAMEmp=workfolder+date+'_mp.xyz'
        if not os.path.exists(workfolder+date+'_mp.xyz'):
            print('write table (xyz-file)')
            write_table(FILENAMEmp,Xi,Yi,mpi)
	FILENAMEow=workfolder+date+'_ow.xyz'
        if not os.path.exists(workfolder+date+'_ow.xyz'):
            print('write table (xyz-file)')
            write_table(FILENAMEow,Xi,Yi,owi)
	FILENAMEsn=workfolder+date+'_sn.xyz'
        if not os.path.exists(workfolder+date+'_sn.xyz'):
            print('write table (xyz-file)')
            write_table(FILENAMEsn,Xi,Yi,sni)

        ctabfilemp=folder+'mp.ctab'   # makecpt -Cjet -D -Z -T0.0/0.4/0.1 > mp.ctab
        ctabfilestdev=folder+'stdev.ctab'
        pos_x,pos_y=str(-3850), str(6350)

	
	if not os.path.exists(workfolder+date+'_mp_05.nc'):
		print 'start processind 0.5 km grid'
		#create 0.5 grid
		#plotfolder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/PRODUCTS/daily/plots/'
		plotfolder=folder+'plots/'
		outfile_mp_05=workfolder+date+'_mp_05.nc'
		mapfile05=plotfolder+date+'_mp_05.ps'
		G={}
		cs_sm=0.5
		region='Arc'
		grid_par_reg(region,cs_sm,G)
		titlestr05=titlestr+' mp 0.5 km grid'
		xyz2grd(Xi,Yi,mpi,outfile_mp_05,G)
		#check plot
	
	#scale to 12.5 km
	titlestr125=titlestr+' 12.5 km'
	if not os.path.exists(workfolder+date+'_mp_125.nc'):
		print 'start processing 12.5 km mp grid'
		G={}	
		cs=12.5
		region='Arc'
		grid_par_reg(region,cs,G)
		titlestr125=titlestr+' mp 12.5 km grid'
		titlestr125SD=titlestr+' SD 12.5 km grid'
		titlestr125W=titlestr+' Weight 12.5 km grid'
		outfile=workfolder+date+'_mp_125.xyz'
		MEANoutfile_mp125=workfolder+date+'_mp_125.nc'
		
		os.system('GMT blockmean '+FILENAMEmp+opts(G,['Rx','I'])+' -Wo -Eslh -V > '+outfile)
		os.system('GMT xyz2grd '+outfile+opts(G,['Rx','I'])+' -V -G'+MEANoutfile_mp125)
		
		MEANfile=workfolder+date+'_mp_125MEAN.xyz'
		SDfile=workfolder+date+'_mp_125SD.xyz'
	
		SDoutfile_mp125=workfolder+date+'_mp_125SD.nc'
		MEANmapfile=plotfolder+date+'_mp_125.ps'
		SDmapfile=plotfolder+date+'_mp_125SD.ps'
		WEIGHTfile=workfolder+date+'_mp_125WEIGHT.xyz'
		WEIGHToutfile_mp125=workfolder+date+'_mp_125WEIGHT.nc'
		WEIGHTmapfile=plotfolder+date+'_mp_125WEIGHT.ps'
		ctabfile25weight=folder+'25weight.ctab'
	
		
	
		cmd('GMT gmtconvert '+outfile+' -F0,1,3 -V > '+SDfile)
		cmd('GMT xyz2grd '+SDfile+opts(G,['Rx','I'])+' -V -G'+SDoutfile_mp125)
		
		cmd('GMT gmtconvert '+outfile+' -F0,1,6 -V > '+WEIGHTfile)
		cmd('GMT xyz2grd '+WEIGHTfile+opts(G,['Rx','I'])+' -V -G'+WEIGHToutfile_mp125)
	
		
	#################ow concentration####################
	if not os.path.exists(workfolder+date+'_ow_05.nc'):
		print 'start processind 0.5 km grid'
		#create 0.5 grid
		#plotfolder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/PRODUCTS/daily/plots/'
		plotfolder=folder+'plots/'
		outfile_ow_05=workfolder+date+'_ow_05.nc'
		mapfile05=plotfolder+date+'_ow_05.ps'
		G={}
		cs_sm=0.5
		region='Arc'
		grid_par_reg(region,cs_sm,G)
		titlestr05=titlestr+' ow 0.5 km grid'
		xyz2grd(Xi,Yi,owi,outfile_ow_05,G)
		#check plot
	
	if not os.path.exists(workfolder+date+'_ow_125.nc'):
		print 'start processing 12.5 km ow grid'
		#scale to 12.5 km
		outfile=workfolder+date+'_ow_125k.xyz'
		outfile_ow125=workfolder+date+'_ow_125.nc'
		mapfile=plotfolder+date+'_ow_125.ps'
		cs=12.5
		G={}
		grid_par_reg(region,cs,G)
		titlestr125=titlestr+' 12.5 km'
		cmd('GMT gmtset OBLIQUE_ANOTATION 1 PAPER_MEDIA A4+ PAGE_ORIENTATION portrait PLOT_DEGREE_FORMAT ddd:mm:ss -V')
		cmd('GMT gmtset COLOR_NAN 200')
		os.system('GMT blockmean '+FILENAMEow+opts(G,['Rx','I'])+' -Eslh -V > '+outfile)
		os.system('GMT xyz2grd '+outfile+opts(G,['Rx','I'])+' -V -G'+outfile_ow125) 
	
	#################sn concentration####################
	if not os.path.exists(workfolder+date+'_sn_05.nc'):
		print 'start processind 0.5 km grid'
		#create 0.5 grid
		#plotfolder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/PRODUCTS/daily/plots/'
		plotfolder=folder+'plots/'
		outfile_sn_05=workfolder+date+'_sn_05.nc'
		mapfile05=plotfolder+date+'_sn_05.ps'
		G={}
		cs_sm=0.5
		region='Arc'
		grid_par_reg(region,cs_sm,G)
		titlestr05=titlestr+' ow 0.5 km grid'
		xyz2grd(Xi,Yi,sni,outfile_sn_05,G)
		#check plot
	
	if not os.path.exists(workfolder+date+'_sn_125.nc'):
		print 'start processing 12.5 km sn grid'
		#scale to 12.5 km
		outfile=workfolder+date+'_sn_125k.xyz'
		outfile_ow125=workfolder+date+'_sn_125.nc'
		mapfile=plotfolder+date+'_sn_125.ps'
		cs=12.5
		G={}
		grid_par_reg(region,cs,G)
		titlestr125=titlestr+' 12.5 km'
		cmd('GMT gmtset OBLIQUE_ANOTATION 1 PAPER_MEDIA A4+ PAGE_ORIENTATION portrait PLOT_DEGREE_FORMAT ddd:mm:ss -V')
		cmd('GMT gmtset COLOR_NAN 200')
		os.system('GMT blockmean '+FILENAMEsn+opts(G,['Rx','I'])+' -Eslh -V > '+outfile)
		os.system('GMT xyz2grd '+outfile+opts(G,['Rx','I'])+' -V -G'+outfile_ow125) 
	
	#############ICE CONCENTRATION###################
	#####  ice=1-ow
	##run ice_konz_regrid.py
	##python ice_konz_regrid.py
	##check here:
	#if not os.path.exists(workfolder+date+'_ice_05.nc'):
		#print('ice concentration set doesnot exist - start processing of '+dataset)
		#ice=1-array(ow,dtype=float32)    #  ice=1-ow
		#print('ice shape:'+str(ice.shape))
		#region='Arc'
		#cs_sm=0.5
		#XYs=XYgrid(region,cs_sm)
		#X=array(XYs[0],dtype=float32)
		#Y=array(XYs[1],dtype=float32)
		#i=isfinite(ice)
		#icei=ice[i]
		#Xi=X[i]
		#Yi=Y[i]
		#pos_x,pos_y=str(-3850), str(6350)
		#ctabfileice=folder+'ice.ctab'
		#titlestr='MODIS ICE CONCENTRATION - '+DATE+'('+str(n)+')'
		#FILENAME=workfolder+date+'_ice.xyz'
		#write_table(FILENAME,Xi,Yi,icei)
	
	
	#if not os.path.exists(workfolder+date+'_ice_25.nc'):
		#print 'start processing 25 km ice grid'
		##scale to 25 km
		#outfile=workfolder+date+'_ice_25k.xyz'
		#outfile_ice25=workfolder+date+'_ice_25.nc'
		#mapfile=plotfolder+date+'_ice_25.ps'
		#cs=25
		#G={}
		#grid_par_reg(region,cs,G)
		##prepare headline:
		#titlestr25=titlestr+' 25 km'
		#cmd('gmtset OBLIQUE_ANOTATION 1 PAPER_MEDIA A4+ PAGE_ORIENTATION portrait PLOT_DEGREE_FORMAT ddd:mm:ss -V')
		#cmd('gmtset COLOR_NAN 200')
		#os.system('blockmean '+FILENAME+opts(G,['Rx','I'])+' -Eslh -V > '+outfile)
		#os.system('xyz2grd '+outfile+opts(G,['Rx','I'])+' -V -G'+outfile_ice25) 
		
	#if not os.path.exists(workfolder+date+'_ice_125.nc'):
		#print 'start processing 12.5 km ice grid'
		##scale to 12.5 km
		#outfile=workfolder+date+'_ice_125k.xyz'
		#outfile_ice125=workfolder+date+'_ice_125.nc'
		#mapfile=plotfolder+date+'_ice_125.ps'
		#cs=12.5
		#G={}
		#grid_par_reg(region,cs,G)
		#titlestr125=titlestr+' 12.5 km'
		#cmd('gmtset OBLIQUE_ANOTATION 1 PAPER_MEDIA A4+ PAGE_ORIENTATION portrait PLOT_DEGREE_FORMAT ddd:mm:ss -V')
		#cmd('gmtset COLOR_NAN 200')
		#os.system('blockmean '+FILENAME+opts(G,['Rx','I'])+' -Eslh -V > '+outfile)
		#os.system('xyz2grd '+outfile+opts(G,['Rx','I'])+' -V -G'+outfile_ice125) 

