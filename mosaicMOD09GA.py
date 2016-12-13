#run mosaicMOD09GA.py
#python mosaicMOD09GA.py


import sys,pipes,struct,os,glob,fnmatch,pickle
from mpl_toolkits.basemap import NetCDFFile
from sattools import *
from scipy.ndimage import *
from gmttools import *

from pylab import *




####---snow1:
year='2011'

workfolder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/MOD09/daily/processed/'
metalist=glob.glob(workfolder+year+'*/')
print(metalist)

daylist=[]

for day in metalist:
    d=day[-11:-8]
    daylist.append(d)

days=list(sorted(set(daylist)))
#####
#####
#days=['147', '148', '149', '150', '151', '152', '153', '154', '155', '156', '157', '158']
print(days)
#############
#BLENDFILE
###############


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
        cmd('gmtset OBLIQUE_ANOTATION 1 PAPER_MEDIA A4+ PAGE_ORIENTATION portrait PLOT_DEGREE_FORMAT ddd:mm:ss   COLOR_NAN 255/255/255')
        cmd('grdblend '+blendfile+' -G'+blendoutfile+' '+opts(G,['Rx','I'])+' -V > '+mapfile)
        cmd('grdimage '+blendoutfile+' '+opts(G,['Rx','Jx','Bx','C'])+' -P -K -V>> '+mapfile)
        cmd('psscale -D8c/-0.5c/12c/0.4ch -C'+ctabfile+' -O -B0.1 -V -K >> '+mapfile)
        cmd('pscoast '+opts(G,['Rll','Bll','coast'])+' -G128/128/128 -O -V -K >> '+mapfile)
        cmd('echo -2000 6500 18 0 0 0 "Overview B'+str(band)+' Day '+str(d)+'" | pstext '+opts(G,['Rx','Jx','N'])+' -V -O >> '+mapfile)
        #cmd('gv '+mapfile+'&')


