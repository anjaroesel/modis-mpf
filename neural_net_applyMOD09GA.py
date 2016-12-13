# applay neural net from programm neural_net_train on full modis scene 
# based onclassifcation Tschudi 2008 for MODIS - haupt_mod_tschudi.py

# run neural_net_applyMOD09GA.py
# python neural_net_applyMOD09GA.py

import sys,pipes,struct,os,glob,fnmatch,pickle
from mpl_toolkits.basemap import NetCDFFile
from scipy.interpolate import *
from sattools import *
from scipy.ndimage import *
import time
from ffnet import ffnet,mlgraph
import networkx
from pylab import *
from scipy.interpolate import *
import scipy.io as io
import scipy.optimize as opti
from mpl_toolkits.basemap import Basemap



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





#####select date:
#date='2010_097'
####---snow1:
workfolder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/GRID_MOSAIC/daily/'
tmpfolder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/tmp/'
############
##METALIST## 
############
metalist = glob.glob(workfolder+'2011_*/')
print(metalist)


for folder in metalist:
    print folder
    jahr=folder[-16:-12]
    julday=folder[-11:-8]
    date=jahr+'_'+julday
    workfile=folder+date+'_mosaic'
    print('work on '+workfile)
    if not os.path.exists('/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/PRODUCTS/daily/'+date+'/'+date+'_ow.data'):
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
        #net=pickle.load(open(tmpfolder+'neural_net_TDneu.data','r'))   ###ergebnisse ok
        net=pickle.load(open(tmpfolder+'neural_net_my_sib_can_fram_MP_valuesTSCHUDI','r'))
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
        cmd('mkdir /scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/PRODUCTS/daily/'+date)
        pickle.dump(mp,open('/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/PRODUCTS/daily/'+date+'/'+date+'_mp.data','w'))
        print('mp pickled')
        pickle.dump(snow,open('/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/PRODUCTS/daily/'+date+'/'+date+'_snow.data','w'))
        print('snow pickled')
        pickle.dump(ow,open('/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/PRODUCTS/daily/'+date+'/'+date+'_ow.data','w'))
        print('ow pickeld......done')


