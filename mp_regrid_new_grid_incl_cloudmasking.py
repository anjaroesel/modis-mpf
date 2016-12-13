#run mp_regrid_new_grid_incl_cloudmasking.py
#python mp_and_ice_regrid_only_rmp_incl_masking.py

import sys,pipes,struct,os,glob,fnmatch,pickle
from scipy.interpolate import *
from sattools import *
#import scipy.ndimage as nd
#import time
from gmttools import *
from pylab import *
import scipy.io as io
from mpl_toolkits.basemap import Basemap
import Nio

#folder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/PRODUCTS/daily/'
folder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/PRODUCTS/'

def JulDay2Date(year,n):
    import datetime
    d = datetime.date(year, 1, 1) + datetime.timedelta(n - 1)
    month=d.month
    day=d.day
    return month, day

############
############
#orig
#metalist = glob.glob(folder+'2007_233*')
#metalist = glob.glob(folder+'2005_169/*mp.data')
metalist = glob.glob(folder+'20*/*_mp_125.nc')
#print(metalist)

#only for testing
#metalist = glob.glob(folder+'2008_169/*mp.data')
print(metalist)

for dataset in metalist:
	print dataset
	#####select date:
	date=dataset[-18:-10]#'2008_169'
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
	workfolder=folder+date+'/'
	workfile=workfolder+date
	tmpfolder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/tmp/'
	plotfolder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/PRODUCTS/plots/'
	jahr=date[0:4]
	julday=date[5:8]
	mpfile=workfile+'_mp_125.nc'
	weightfile=workfile+'_mp_125WEIGHT.nc'
	mp,weight=NetCDFFile(mpfile),NetCDFFile(weightfile)
	X=array(mp.variables['x'][:])
	Y=array(mp.variables['y'][:])
	X,Y=meshgrid(X,Y)
	cs=12.5
	region='Arc'
	XYs=XYgrid(region,cs)
	X=array(XYs[0],dtype=float32)
	Y=array(XYs[1],dtype=float32)
	MP=array(mp.variables['z'][:,:])
	MP=MP[:-1,:-1]
	WEIGHT=array(weight.variables['z'][:,:])
	WEIGHT=WEIGHT[:-1,:-1]
	#alle Pixel mit weniger als 50% gewicht auf 0 setzen
	WEIGHT[WEIGHT<nanmax(WEIGHT)*.5]=0
	WEIGHT[WEIGHT>=nanmax(WEIGHT)*.5]=1
	MP=MP*WEIGHT
	MP[MP==0]=NAN
	FILENAME=workfolder+date+'_mp_masked.xyz'
	ctabfilemp=folder+'mp.ctab'   # makecpt -Cjet -D -Z -T0.0/0.4/0.1 > mp.ctab
	pos_x,pos_y=str(-3850), str(6350)
	titlestr='MODIS MP masked - '+DATE+'('+str(n)+')'
	#G={}
	#grid_par_reg(region,cs,G)
	
	#os.system('mkdir '+plotfolder+year)
	#pltfolder=plotfolder+year+'/'
	
	#figure()
	#imshow(MP,interpolation='Nearest',origin='lower',vmin=0.,vmax=1)
	#colorbar()
	#title('Masked MP '+DATE)
	#savefig(pltfolder+'Masked_MP'+jahr+julday+'.png')
	ctabfilemp=folder+'mp.ctab'   # makecpt -Cjet -D -Z -T0.0/0.4/0.1 > mp.ctab
	
	
	
	print 'start processing masked'
	outfile_mp_125masked=workfolder+date+'_mp_125masked.nc'
	mapfileM=plotfolder+date+'_mp_125masked66.ps'
	G={}
	grid_par_reg(region,cs,G)
	xyz2grd(X,Y,MP,outfile_mp_125masked,G)
	#check plot
	cmd('gmtset OBLIQUE_ANOTATION 1 PAPER_MEDIA A4+ PAGE_ORIENTATION portrait PLOT_DEGREE_FORMAT ddd:mm:ss -V')
	cmd('gmtset COLOR_NAN 200')
	cmd('grdimage '+outfile_mp_125masked+' '+opts(G,['Rx','Jx','Bx'])+' -C'+ctabfilemp+' -P -K -V> '+mapfileM)
	cmd('psscale -D8c/-0.5c/12c/0.4ch -C'+ctabfilemp+' -Ef -O -B0.1 -V -K >> '+mapfileM)
	cmd('pscoast '+opts(G,['Rll','Bll','coast'])+' -G128/128/128 -O -V -K >> '+mapfileM)    
	cmd('echo '+pos_x+' '+pos_y+' 18 0 0 0 "'+titlestr+'" | pstext '+opts(G,['Rx','Jx','N'])+' -V -O >> '+mapfileM)
	cmd('gv '+mapfileM+'&')
	
	
	
