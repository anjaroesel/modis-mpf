#run plot_AMSRE_MODIS_IC_timeseries.py

from scipy import *
import sys,pipes,struct,os,glob,fnmatch,pickle
from pylab import *
from scipy.stats import nanmean, nanstd
from mpl_toolkits.basemap import Basemap
from polar_projection import *
import scipy.io as io
from mpl_toolkits.basemap import Basemap
import Nio
from sattools import *
from gmttools import *
import scipy.ndimage as ndi
from pyhdf import SD
import datetime as dt



def JulDay2Date(year,n):
    import datetime
    d = dt.date(year, 1, 1) + dt.timedelta(n - 1)
    month=d.month
    day=d.day
    return month, day

def XYgrid2(region,cs):
    x0,y0,x2,y2=x0y0x2y2(region)
    y=linspace(y2+cs/2,y0-cs/2,(y0-y2)/cs)
    x=linspace(x2+cs/2,x0-cs/2,(x0-x2)/cs)
    return x,y





####map setup
lon_0=-45. #center coordinates of the plot
lat_0=90.0

####cut plot region #canadian achipelago
latmin,lonmin=65,-100.
latmax,lonmax=80,-110.

##Meltex
#latmin,lonmin=69.,-135.
#latmax,lonmax=73.,-146.


Xmax,Ymax=mapll(latmax,lonmax,1)
Xmin,Ymin=mapll(latmin,lonmin,1)


############################################
#read MODIS data
####---snow1:
####---snow1:
workfolder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/PRODUCTS/'
tmpfolder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/tmp/'
AMSRE1folder='/scratch/clisap/seaice/AMSR-E/L3polargrids/89_6250/'  #ASI
AMSRE2folder='/scratch/clisap/seaice/AMSR-E/L3polargrids/all_12500/'  #NSIDC
############
##METALIST##
############/
metalist = glob.glob(workfolder+'2011_169*')
print metalist

#folder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/PRODUCTS/'
for m in metalist:
	year=m[-8:-4]
	doy=m[-3:]
	month,day=JulDay2Date(int(year),int(doy))
	month=str(0)+str(month)
	if len(str(day))<2:
		day=str(0)+str(day)
		print 'day ist einstellig'
		
	if len(str(day))==2:
		day=str(day)
		print 'day ist zweistellig'
		
	print year, month, day
	##read modis ice and mp data
	owdataset=m+'/'+year+'_'+doy+'_ow_125.nc'
	OW_obj=NetCDFFile(owdataset)
	OW=array(OW_obj.variables['z'][:])
	ICE=1-array(OW_obj.variables['z'][:])
	ICEX=array(OW_obj.variables['x'][:])
	ICEY=array(OW_obj.variables['y'][:])
	iXmax=nonzero(ICEX>=Xmax)[0][0]
	iXmin=nonzero(ICEX>=Xmin)[0][0]
	iYmax=nonzero(ICEY>=Ymax)[0][0]
	iYmin=nonzero(ICEY>=Ymin)[0][0]
	ICEX=ICEX[iXmin:iXmax]
	ICEY=ICEY[iYmin:iYmax]
	#OW=OW[iYmin:iYmax,iXmin:iXmax]
	ICE=ICE[iYmin:iYmax,iXmin:iXmax]*100
	mpdataset=m+'/'+year+'_'+doy+'_mp_125.nc'
	print mpdataset
	m_obj=NetCDFFile(mpdataset)
	MP=array(m_obj.variables['z'][:])
	MX=array(m_obj.variables['x'][:])
	MY=array(m_obj.variables['y'][:])
	iXmax=nonzero(MX>=Xmax)[0][0]
	iXmin=nonzero(MX>=Xmin)[0][0]
	iYmax=nonzero(MY>=Ymax)[0][0]
	iYmin=nonzero(MY>=Ymin)[0][0]
	MX=MX[iXmin:iXmax]
	MY=MY[iYmin:iYmax]
	MP=MP[iYmin:iYmax,iXmin:iXmax]*100.
	#ICE=ICE+MP
	############################################
	#read AMSRE data
	asi_grid=6.25
	nasa_grid=12.5
	asi_path = '/net/eiszeit/qfs2/icdc/DATA/ice_and_snow/asi_iceconc/DATA/'
	#asi_path  = '/data/icdc/ice_and_snow/asi_iceconc/DATA/'
	AMSREfile=asi_path+year+'/'+'asi-n6250-'+year+month+day+'-v5i.nc'
	AMSREfile2gz=AMSRE2folder+year+'.'+month+'.'+day+'/AMSR_E_L3_SeaIce12km_V13_'+year+month+day+'.hdf.gz'
	AMSREfile2=workfolder+'AMSR_E_L3_SeaIce12km_V13_'+year+month+day+'.hdf'
	###unpacking
	os.system('cp '+AMSREfile2gz+' '+workfolder)
	os.system('uncompress '+AMSREfile2+'.gz '+AMSREfile2)
	##ASI
	obj=NetCDFFile(AMSREfile)
	asi=array(obj.variables['icecon'][:], dtype=float32)/10
	asi[asi==-100]=nan
	
	###NASA-TEAM2
	AMSRE_NASA=SD.SD(AMSREfile2)
	nt2_key=26#'SI_12km_NH_89V_DSC'
	bt_key=29#'SI_12km_NH_ICECON_DSC'
	sd_ice_nt2   = AMSRE_NASA.select(nt2_key)
	ice_nt2      = sd_ice_nt2.get()
	ice_nt2      = array(ice_nt2,dtype=float)[::-1]
	ice_nt2[ice_nt2>=101]=nan


	###NASA-BOOTSTRAP
	sd_ice_bt   = AMSRE_NASA.select(bt_key)
	ice_bt      = sd_ice_bt.get()    # stored as difference to NT2 data (ice_bst-ice_nt2)
	ice_bt      = array(ice_bt,dtype=float)[::-1]
	ice_bt      = ice_bt + ice_nt2
	ice_bt[ice_bt>=101]=nan

	#figure()
	#imshow(ice_nt2)
	#figure()
	#imshow(ice_bt)
	#figure()
	#imshow(asi)


	region='Arc'
	cs_asi=asi_grid
	AX,AY=XYgrid2(region,cs_asi)
	x0,x2=AX.max(),AX.min()
	y0,y2=AY.max(),AY.min()

	iXmax=nonzero(AX>=Xmax)[0][0]
	iXmin=nonzero(AX>=Xmin)[0][0]
	iYmax=nonzero(AY>=Ymax)[0][0]
	iYmin=nonzero(AY>=Ymin)[0][0]

	AX=AX[iXmin:iXmax]
	AY=AY[iYmin:iYmax]
	asi=asi[iYmin:iYmax,iXmin:iXmax]
	cs=nasa_grid
	AX,AY=XYgrid2(region,cs)
	x0,x2=AX.max(),AX.min()
	y0,y2=AY.max(),AY.min()

	iXmax=nonzero(AX>=Xmax)[0][0]
	iXmin=nonzero(AX>=Xmin)[0][0]
	iYmax=nonzero(AY>=Ymax)[0][0]
	iYmin=nonzero(AY>=Ymin)[0][0]
	
	AX=AX[iXmin:iXmax]
	AY=AY[iYmin:iYmax]
	ice_nt2=ice_nt2[iYmin:iYmax,iXmin:iXmax]
	ice_bt=ice_bt[iYmin:iYmax,iXmin:iXmax]
	#OWmask=ice_bt.copy()    ####wassermaske
	#OWmask[OWmask<15.0]=0
	#OWmask[OWmask>15.0]=1
	#ICEmasked=ICE*OWmask

	
	
	Tmap=Basemap(resolution='h',projection='stere', boundinglat=60.,llcrnrlon=lonmin, llcrnrlat=latmin, urcrnrlon=lonmax, urcrnrlat=latmax, lon_0=lon_0, lat_0=lat_0)
	
	
	figure(figsize=(15,10))
	title(day+month+year)
	subplot(2,3,1)
	Tmap.imshow(asi, interpolation='nearest',aspect='auto',cmap='bone',vmin=0,vmax=100)
	colorbar(orientation='vertical', shrink=0.8, ticks=(0,20,40,60,80,100))
	Tmap.drawcoastlines(color='black')
	Tmap.fillcontinents(color='gray')
	Tmap.drawmeridians(np.arange(0,360,10),labels=[1,0,0,0])
	Tmap.drawparallels(np.arange(-90,90,5),labels=[0,0,0,1])
	text(300000.0,-550000.0,'ice concentration [%]')
	text(400000.0,1250000.0,'AMSR-E ASI - '+year+month+day)

	subplot(2,3,2)
	Tmap.imshow(ice_nt2, interpolation='nearest',aspect='auto',cmap='bone',vmin=0,vmax=100)
	colorbar(orientation='vertical', shrink=0.8, ticks=(0,20,40,60,80,100))
	Tmap.drawcoastlines(color='black')
	Tmap.fillcontinents(color='gray')
	Tmap.drawmeridians(np.arange(0,360,10),labels=[1,0,0,0])
	Tmap.drawparallels(np.arange(-90,90,5),labels=[0,0,0,1])
	#text(-00000.0,-550000.0,'ice concentration [%]')
	text(400000.0,1250000.0,'AMSR-E NT2 ')#+day+'.'+month+'.'+year)
	subplot(2,3,3)
	Tmap.imshow(ice_bt, interpolation='nearest',aspect='auto',cmap='bone',vmin=0,vmax=100)
	colorbar(orientation='vertical', shrink=0.8, ticks=(0,20,40,60,80,100))
	Tmap.drawcoastlines(color='black')
	Tmap.fillcontinents(color='gray')
	Tmap.drawmeridians(np.arange(0,360,10),labels=[1,0,0,0])
	Tmap.drawparallels(np.arange(-90,90,5),labels=[0,0,0,1])
	#text(-00000.0,-550000.0,'ice concentration [%]')
	text(400000.0,1250000.0,'AMSR-E BT ')#+day+'.'+month+'.'+year)
	
	
	
	
	subplot(2,3,4)
	Tmap.imshow(ICE, interpolation='nearest',cmap='bone',vmin=0,vmax=100)
	colorbar(orientation='vertical', shrink=0.8, ticks=(0,20,40,60,80,100))
	Tmap.drawcoastlines(color='black')
	Tmap.fillcontinents(color='gray')
	Tmap.drawmeridians(np.arange(0,360,10),labels=[1,0,0,0])
	Tmap.drawparallels(np.arange(-90,90,5),labels=[0,0,0,1])
	#text(-00000.0,-550000.0,'ice concentration [%]')
	text(300000.0,1250000.0,'MODIS SEA ICE ')#+day+'.'+month+'.'+year)
	subplot(2,3,5)
	Tmap.imshow(MP, interpolation='nearest',cmap='jet',vmin=0,vmax=70)
	colorbar(orientation='vertical', shrink=0.8, ticks=(0,20,40,60))
	Tmap.drawcoastlines(color='black')
	Tmap.fillcontinents(color='gray')
	Tmap.drawmeridians(np.arange(0,360,10),labels=[1,0,0,0])
	Tmap.drawparallels(np.arange(-90,90,5),labels=[0,0,0,1])
	text(200000.0,-550000.0,'melt pond fraction [%]')
	text(300000.0,1250000.0,'MODIS MP ')#+day+'.'+month+'.'+year)
	savefig('/pf/u/u241127/plots/amsre_modis_ic_comp_'+year+'_'+month+'_'+day+'_allAMSRE.png')
	###aufraeumen
	os.system('rm '+AMSREfile2)


NT2_mean=mean(nanmean(ice_nt2))
BT_mean=mean(nanmean(ice_bt))
ASI_mean=mean(nanmean(asi))
MODICE=mean(nanmean(ICE))
M_P=mean(nanmean(MP))

print( 'nt2_mean,bt_mean, asi_mean,modis_ice,mp: '+ str(NT2_mean)+','+str(BT_mean)+','+str(ASI_mean)+','+str(MODICE)+','+str(M_P) )
#94.1277752667 93.0629555122 51.2498602873

M=ICE[50:70,40:60]
a=asi[100:140,80:120]
nt2=ice_nt2[50:70,40:60]
bt=ice_bt[50:70,40:60]
mp=MP[50:70,40:60]

figure()
imshow(M, vmin=0,vmax=100,origin='lower',cmap='bone')
title('MODIS')
colorbar()
figure()
imshow(a,vmin=0,vmax=100,origin='lower',cmap='bone')
title('ASI')
colorbar()
figure()
imshow(nt2,vmin=0,vmax=100,origin='lower',cmap='bone')
title('NT2')
colorbar()
figure()
imshow(bt,vmin=0,vmax=100,origin='lower',cmap='bone')
title('BT')
colorbar()
figure()
imshow(mp,vmin=0,vmax=70,origin='lower',cmap='jet')
title('Mp')
colorbar()





nt2_mean=mean(nanmean(nt2))
bt_mean=mean(nanmean(bt))
asi_mean=mean(nanmean(a))
modis=mean(nanmean(M))
M_P=mean(nanmean(mp))

print( 'nt2_mean,bt_mean, asi_mean,modis,mp: '+ str(nt2_mean)+','+str(bt_mean)+','+str(asi_mean)+','+str(modis)+','+str(M_P))
#94.1277752667 93.0629555122 51.2498602873
