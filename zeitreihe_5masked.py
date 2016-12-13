#run zeitreihe_5masked.py

from scipy import *
import sys,pipes,struct,os,glob,fnmatch,pickle
sys.path.append('/pf/u/u241127/progs/')
from pylab import *
from mpl_toolkits.basemap import NetCDFFile
from scipy.stats import linregress
from scipy.stats import nanmean
from scipy.interpolate import *
from polar_projection import *
from sattools import *
from gmttools import *
import time as tm
import caselist
from mpl_toolkits.basemap import NetCDFFile
import datetime as dt
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.cm as cm
import scipy.stats

def JulDay2Date(year,n):
    import datetime
    d = dt.date(year, 01, 01) + dt.timedelta(n - 1)
    month=str(d.month)
    if int(month)<10:
	month='0'+month
	
    day=str(d.day)
    if int(day)<10:
	 day='0'+day
	 
    year=str(year)
    return year, month, day
    
    
##kack signifikanztest nach Cox und stuart
##L. Sachs, angewandte statistik, S. 47

def cox_stuart_sig_test(t1):
	f=isfinite(t1)
	t1=t1[f]
	if len(t1)>=10:
		if len(t1)%2 == 0:
			print 'i gerade'
			i=int(len(t1)/2)
			fhalf=t1[0:i]
			shalf=t1[i:]
			
		if len(t1)%2 != 0:
			print'i ungerade'
			i=int(len(t1)/2)+1
			fhalf=t1[0:i]
			shalf=t1[i*-1:]
						
		differenz=fhalf-shalf
		signs=sign(differenz)
		signs = signs[ signs != 0 ] #eliminate 0s
		pos = len(signs[signs > 0])
		neg = len(signs[signs < 0])
		if pos > neg:
			z=(abs(pos-22./6.)-0.5)/sqrt(len(t1)/12)
		if pos< neg:
			z=((abs(neg-22./6.)-0.5)/sqrt(len(t1)/12))*-(1.)
		if pos==neg:
			z=0.
	if len(t1)<10:
		z=0
	return z


colordict = {'red':  ((0.0, 0.0, 0.0),
                   (0.25,0.0, 0.0),
                   (0.5, 0.8, 1.0),
                   (0.75,1.0, 1.0),
                   (1.0, 0.4, 1.0)),

         'green': ((0.0, 0.0, 0.0),
                   (0.25,0.0, 0.0),
                   (0.5, 0.9, 0.9),
                   (0.75,0.0, 0.0),
                   (1.0, 0.0, 0.0)),

         'blue':  ((0.0, 0.0, 0.4),
                   (0.25,1.0, 1.0),
                   (0.5, 1.0, 0.8),
                   (0.75,0.0, 0.0),
                   (1.0, 0.0, 0.0))
        }



########### subregion A
homefolder='/pf/u/u241127/data/'
meantable=homefolder+'mp_masked_new_lat_table_bearb.csv'

mt=loadtxt(meantable,skiprows=1,delimiter=',',dtype=float32)
#mt[mt==0.000]=nan


years=mt[:,0]
days=mt[:,2]
months=mt[:,1]
degrees=mt[:,3]
meanweek=mt[:,4]
stddevweek=mt[:,5]



m=meanweek.reshape((192,30)).transpose()[::-1]



y2000=m[:,0:16]
y2001=m[:,16:32]
y2002=m[:,32:48]
y2003=m[:,48:64]
y2004=m[:,64:80]
y2005=m[:,80:96]
y2006=m[:,96:112]
y2007=m[:,112:128]
y2008=m[:,128:144]
y2009=m[:,144:160]
y2010=m[:,160:176]
y2011=m[:,176:192]







#mean errechnen
AA=y2000.copy()
AA[AA>0.]=1.
BB=y2001.copy()
BB[BB>0.]=1.
CC=y2002.copy()
CC[CC>0.]=1.
DD=y2003.copy()
DD[DD>0.]=1.
EE=y2004.copy()
EE[EE>0.]=1.
FF=y2005.copy()
FF[FF>0.]=1.
GG=y2006.copy()
GG[GG>0.]=1.
HH=y2007.copy()
HH[HH>0.]=1.
II=y2008.copy()
II[II>0.]=1.
KK=y2009.copy()
KK[KK>0.]=1.
LL=y2010.copy()
LL[LL>0.]=1.
MM=y2011.copy()
MM[MM>0.]=1.
summe=AA+BB+CC+DD+EE+FF+GG+HH+II+KK+LL+MM
summe[summe==0]=1
mean_all=(AA*y2000+BB*y2001+CC*y2002+DD*y2003+EE*y2004+FF*y2005+GG*y2006+HH*y2007+II*y2008+KK*y2009+LL*y2010+MM*y2011)/summe



mmm2000=AA*y2000
mmm2001=BB*y2001
mmm2002=CC*y2002
mmm2003=DD*y2003
mmm2004=EE*y2004
mmm2005=FF*y2005
mmm2006=GG*y2006
mmm2007=HH*y2007
mmm2008=II*y2008
mmm2009=KK*y2009
mmm2010=LL*y2010
mmm2011=MM*y2011


years=[y2000,y2001,y2002,y2003,y2004,y2005,y2006,y2007,y2008,y2009,y2010,y2011]
myears=[mmm2000,mmm2001,mmm2002,mmm2003,mmm2004,mmm2005,mmm2006,mmm2007,mmm2008,mmm2009,mmm2010,mmm2011]



params = {'axes.labelsize': 13,
          'text.fontsize': 16,
          'legend.fontsize': 16,
          'xtick.labelsize': 13,
          'ytick.labelsize': 16,}
rcParams.update(params)



###############PLOTS...............

date1 = dt.date( 2000, 5, 8 )
date2 = dt.date( 2000, 9, 6 )
delta = dt.timedelta(days=8)
dates = drange(date1, date2, delta)

months   = MonthLocator()  # every month
yearsFmt = DateFormatter('%Y')
blue_red = LinearSegmentedColormap('BlueRed', colordict)



###### Mean#########
##############################
mean_all[mean_all==0]=nan

fig = figure()
fig.set_dpi(100)
fig.set_size_inches((6.0,6.0),forward=True)
ax = fig.add_subplot(111)
plt.subplots_adjust(left=0.10,bottom=0.10,right=0.95,top=0.95) 

####ax.axes([0.1,0.1,0.9,0.9])
#ax.set_ylim((61,90))
# Zuerst die x-Ticks richtig einstellen
y_orte = arange(65,90,5) 
y_werte =('65$\degree$N','70$\degree$N','75$\degree$N','80$\degree$N','85$\degree$N','90$\degree$N')
yticks(y_orte,y_werte)
ax.grid()
ax.set_xlim(dates[0], dates[-1])
imshow(mean_all, interpolation='nearest',extent=([dates[0], dates[-1],61, 90 ]),aspect='auto', vmin=0.,vmax=0.4)
ax.set_autoscale_on(False)
ax.xaxis.set_major_formatter( DateFormatter('%m.%d.') )
colorbar()
#text(08.15, 65, 'melt pond fraction',rotation='vertical', fontsize=14)  #x,y,s
title('MEAN - zeitreihe_4.py')
savefig('/pf/u/u241127/plots/mp_masked_latitudonal_mean_allyears.png')


#### 2007 ###################
#############################
mmm2007[mmm2007==0]=nan

fig = figure()
fig.set_dpi(100)
fig.set_size_inches((6.0,6.0),forward=True)
ax = fig.add_subplot(111)
plt.subplots_adjust(left=0.10,bottom=0.10,right=0.95,top=0.95) 
y_orte = arange(65,90,5) 
y_werte =('65$\degree$N','70$\degree$N','75$\degree$N','80$\degree$N','85$\degree$N','90$\degree$N')
yticks(y_orte,y_werte)
ax.grid()
ax.set_xlim(dates[0], dates[-1])
imshow(mmm2007, interpolation='nearest',extent=([dates[0], dates[-1],61, 90 ]),aspect='auto', vmin=0.,vmax=0.4)
ax.set_autoscale_on(False)
ax.xaxis.set_major_formatter( DateFormatter('%m.%d.') )
colorbar()
#text(08.15, 65, 'melt pond fraction',rotation='vertical', fontsize=14)  #x,y,s
title('2007 - zeitreihe.py')
savefig('/pf/u/u241127/plots/mp_masked_latitudonal_2007.png')



#### 2011 ####################
##############################
mmm2011[mmm2011==0]=nan

fig = figure()
fig.set_dpi(100)
fig.set_size_inches((6.0,6.0),forward=True)
ax = fig.add_subplot(111)
plt.subplots_adjust(left=0.10,bottom=0.10,right=0.95,top=0.95) 
y_orte = arange(65,90,5) 
y_werte =('65$\degree$N','70$\degree$N','75$\degree$N','80$\degree$N','85$\degree$N','90$\degree$N')
yticks(y_orte,y_werte)
ax.grid()
ax.set_xlim(dates[0], dates[-1])
imshow(mmm2011, interpolation='nearest',extent=([dates[0], dates[-1],61, 90 ]),aspect='auto', vmin=0.,vmax=0.4)
ax.set_autoscale_on(False)
ax.xaxis.set_major_formatter( DateFormatter('%m.%d.') )
colorbar()
#text(08.15, 65, 'melt pond fraction',rotation='vertical', fontsize=14)  #x,y,s
title('2011 - zeitreihe_4.py')
savefig('/pf/u/u241127/plots/mp_masked_latitudonal_2011.png')

####differenzen#######################
######################################
blue_red = LinearSegmentedColormap('BlueRed', colordict)

diff2007=mmm2007-mean_all
diff2011=mmm2011-mean_all
####d2011
fig = figure()
fig.set_dpi(100)
fig.set_size_inches((6.0,6.0),forward=True)
ax = fig.add_subplot(111)
plt.subplots_adjust(left=0.10,bottom=0.10,right=0.95,top=0.95) 
y_orte = arange(65,90,5) 
y_werte =('65$\degree$N','70$\degree$N','75$\degree$N','80$\degree$N','85$\degree$N','90$\degree$N')
yticks(y_orte,y_werte)
ax.grid()
ax.set_xlim(dates[0], dates[-1])
imshow(diff2011, interpolation='nearest',extent=([dates[0], dates[-1],61, 90 ]),aspect='auto', cmap=blue_red, vmax=0.1, vmin=-0.1)
ax.set_autoscale_on(False)
ax.xaxis.set_major_formatter( DateFormatter('%m.%d.') )
colorbar()
#text(08.15, 65, 'melt pond fraction',rotation='vertical', fontsize=14)  #x,y,s
title('Differenz 2011-mean - zeitreihe_4.py')
savefig('/pf/u/u241127/plots/mp_masked_latitudonal_diff_2011.png')
####d2007
fig = figure()
fig.set_dpi(100)
fig.set_size_inches((6.0,6.0),forward=True)
ax = fig.add_subplot(111)
plt.subplots_adjust(left=0.10,bottom=0.10,right=0.95,top=0.95) 
y_orte = arange(65,90,5) 
y_werte =('65$\degree$N','70$\degree$N','75$\degree$N','80$\degree$N','85$\degree$N','90$\degree$N')
yticks(y_orte,y_werte)
ax.grid()
ax.set_xlim(dates[0], dates[-1])
imshow(diff2007, interpolation='nearest',extent=([dates[0], dates[-1],61, 90 ]),aspect='auto', cmap=blue_red, vmax=0.1, vmin=-0.1)
ax.set_autoscale_on(False)
ax.xaxis.set_major_formatter( DateFormatter('%m.%d.') )
colorbar()
#text(08.15, 65, 'melt pond fraction',rotation='vertical', fontsize=14)  #x,y,s
title('Differenz 2007-mean - zeitreihe_5.py')
savefig('/pf/u/u241127/plots/mp_masked_latitudonal_diff_2007.png')









####for trend calculation
mt=loadtxt(meantable,skiprows=1,delimiter=',',dtype=float32)
mt[mt==0.000]=nan
meanweek=mt[:,4]
m=meanweek.reshape((192,30)).transpose()[::-1]
y2000=m[:,0:16]
y2001=m[:,16:32]
y2002=m[:,32:48]
y2003=m[:,48:64]
y2004=m[:,64:80]
y2005=m[:,80:96]
y2006=m[:,96:112]
y2007=m[:,112:128]
y2008=m[:,128:144]
y2009=m[:,144:160]
y2010=m[:,160:176]
y2011=m[:,176:192]

years=[y2000,y2001,y2002,y2003,y2004,y2005,y2006,y2007,y2008,y2009,y2010,y2011]
d129=[]
d137=[]
d145=[]
d153=[]
d161=[]
d169=[]
d177=[]
d185=[]
d193=[]
d201=[]
d209=[]
d217=[]
d225=[]
d233=[]
d241=[]
d249=[]
days=[d129,d137,d145,d153,d161,d169,d177,d185,d193,d201,d209,d217,d225,d233,d241,d249]

for y in years:
	print y
	for i, d in enumerate(days):
		print i,d
		d.append(y[:,i])
		
for d in days:
	d=array(d)

full_values=zeros((12,30,16))
trend=zeros((30,16))
signifikanz=zeros((30,16))
deltaT1=zeros((30,16))

for n,y in enumerate(years):
	print n,y
	full_values[n]=y
	
	
for i in range(trend.shape[0]):
	for j in range(trend.shape[1]):
		t1=full_values[:,i,j]
		f=isfinite(t1)
		t1=t1[f]
		if t1.size == 0:
			print 'yes'
			t=arange(2)
			t1=array([nan,nan,nan])
		t=arange(t1.shape[0])
		(ar,br)=polyfit(t,t1,1)
		xr=polyval([ar,br],t)
		(a1,b1,r1,tt1,s1)=linregress(t,xr)
		#signifikanz[i,j]=r1
		trend[i,j]=a1
		#Trendcalculation:
		deltaT1[i,j]=a1
		first=t1[0]
		T_lat=100-((first-deltaT1)*100/first)
		
trend_all=trend.copy()



		
		



##kack signifikanztest nach Cox und stuart
##L. Sachs, angewandte statistik, S. 47
z95=1.645
z90=1.3
z85=0.8

	

Z=zeros((30,16))

for n,y in enumerate(years):
	print n,y
	full_values[n]=y

for i in range(Z.shape[0]):
	for j in range(Z.shape[1]):
		t1=full_values[:,i,j]
		print t1
		z=abs(cox_stuart_sig_test(t1))
		print z
		Z[i,j]=z
		
		
		

#85% Signifikanz level
Z[Z>z85]=85.
Z[Z<z85]=0.
Z85=Z.copy()


Z=zeros((30,16))
for n,y in enumerate(years):
	print n,y
	full_values[n]=y

for i in range(Z.shape[0]):
	for j in range(Z.shape[1]):
		t1=full_values[:,i,j]
		print t1
		z=abs(cox_stuart_sig_test(t1))
		print z
		Z[i,j]=z
		
		
		

#95% Signifikanz level
Z[Z>z95]=95.
Z[Z<z95]=0.
Z95=Z.copy()

Z=zeros((30,16))
for n,y in enumerate(years):
	print n,y
	full_values[n]=y

for i in range(Z.shape[0]):
	for j in range(Z.shape[1]):
		t1=full_values[:,i,j]
		print t1
		z=abs(cox_stuart_sig_test(t1))
		print z
		Z[i,j]=z
		
		
		

#90% Signifikanz level
Z[Z>z90]=90.
Z[Z<z90]=0.
Z90=Z.copy()


figure()
imshow(Z85)
colorbar()
figure()
imshow(Z95)
colorbar()
figure()
imshow(Z90)
colorbar()
figure()
imshow(trend*100,vmin=-1.6,vmax=1.6)
colorbar()


Z95[Z95==95]=1
trend95=trend_all*Z95
trend95=trend95*100 #in Prozent



Z85[Z85==85]=1
trend85=trend_all*Z85
trend85=trend85*100 #in Prozent
Z90[Z90==90]=1
trend90=trend_all*Z90
trend90=trend90*100 #in Prozent
###### Trend
##############################


figure()
imshow(T_lat,vmin=-10,vmax=10)
colorbar()


fig = figure()
fig.set_dpi(100)
fig.set_size_inches((6.0,6.0),forward=True)
ax = fig.add_subplot(111)
plt.subplots_adjust(left=0.10,bottom=0.10,right=0.95,top=0.95) 
y_orte = arange(65,90,5) 
y_werte =('65$\degree$N','70$\degree$N','75$\degree$N','80$\degree$N','85$\degree$N','90$\degree$N')
yticks(y_orte,y_werte)
ax.grid()
ax.set_xlim(dates[0], dates[-1])
imshow(T_lat, interpolation='nearest',extent=([dates[0], dates[-1],61, 90 ]),aspect='auto', cmap=cm.PRGn_r, vmin=-4,vmax=4)
ax.set_autoscale_on(False)
ax.xaxis.set_major_formatter( DateFormatter('%m.%d.') )
colorbar()
#text(08.15, 65, 'melt pond fraction',rotation='vertical', fontsize=14)  #x,y,s
#title('Trend 95% - zeitreihe_5masked.py')
savefig('/pf/u/u241127/plots/mp_masked _latitudonal_trend_allyears.png')

fig = figure()
fig.set_dpi(100)
fig.set_size_inches((6.0,6.0),forward=True)
ax = fig.add_subplot(111)
plt.subplots_adjust(left=0.10,bottom=0.10,right=0.95,top=0.95) 
y_orte = arange(65,90,5) 
y_werte =('65$\degree$N','70$\degree$N','75$\degree$N','80$\degree$N','85$\degree$N','90$\degree$N')
yticks(y_orte,y_werte)
ax.grid()
ax.set_xlim(dates[0], dates[-1])
imshow(trend90, interpolation='nearest',extent=([dates[0], dates[-1],61, 90 ]),aspect='auto', cmap=cm.PRGn_r, vmin=-1.6,vmax=1.6)
ax.set_autoscale_on(False)
ax.xaxis.set_major_formatter( DateFormatter('%m.%d.') )
colorbar()
#text(08.15, 65, 'melt pond fraction',rotation='vertical', fontsize=14)  #x,y,s
title('Trend 90% - zeitreihe_5masked.py')
savefig('/pf/u/u241127/plots/mp_masked _latitudonal_trend_allyears90.png')

fig = figure()
fig.set_dpi(100)
fig.set_size_inches((6.0,6.0),forward=True)
ax = fig.add_subplot(111)
plt.subplots_adjust(left=0.10,bottom=0.10,right=0.95,top=0.95) 
y_orte = arange(65,90,5) 
y_werte =('65$\degree$N','70$\degree$N','75$\degree$N','80$\degree$N','85$\degree$N','90$\degree$N')
yticks(y_orte,y_werte)
ax.grid()
ax.set_xlim(dates[0], dates[-1])
imshow(trend85, interpolation='nearest',extent=([dates[0], dates[-1],61, 90 ]),aspect='auto', cmap=cm.PRGn_r, vmin=-1.6,vmax=1.6)
ax.set_autoscale_on(False)
ax.xaxis.set_major_formatter( DateFormatter('%m.%d.') )
colorbar()
#text(08.15, 65, 'melt pond fraction',rotation='vertical', fontsize=14)  #x,y,s
title('Trend 85% - zeitreihe_5masked.py')
savefig('/pf/u/u241127/plots/mp_masked _latitudonal_trend_allyears85.png')



###zeitreihe 80-87 deg N
#m_nord_allfillup=mean_all[2:10,:5].flatten()

m_nord2000=nanmean(y2000[2:10,:].flatten())
m_nord2001=nanmean(y2001[2:10,:].flatten())
m_nord2002=nanmean(y2002[2:10,:].flatten())
m_nord2003=nanmean(y2003[2:10,:].flatten())
m_nord2004=nanmean(y2004[2:10,:].flatten())
m_nord2005=nanmean(y2005[2:10,:].flatten())
m_nord2006=nanmean(y2006[2:10,:].flatten())
m_nord2007=nanmean(y2007[2:10,:].flatten())

m_nord2008=nanmean(y2008[2:10,:].flatten())
m_nord2009=nanmean(y2009[2:10,:].flatten())
m_nord2010=nanmean(y2010[2:10,:].flatten())
m_nord2011=nanmean(y2011[2:10,:].flatten())

m_nord=array([m_nord2000,m_nord2001,m_nord2002,m_nord2003,m_nord2004,m_nord2005,m_nord2006,m_nord2007,m_nord2008,m_nord2009,m_nord2010,m_nord2011])



d_short=['2000-08-04','2001-08-04','2002-08-04','2003-07-03','2004-07-11','2005-06-25','2006-07-27','2007-06-09','2008-07-03','2009-08-04','2010-06-25','2011-06-25',]  #12

dd_short=[]
for i in range(12):
	d=dt.datetime.strptime(d_short[i], "%Y-%m-%d")
	dd_short.append(d)

dates_short=array(dd_short)
#fit
############Trendberechnung relativ mp fraction mp_nord
t=arange(12)


#Linear regressison -polyfit - polyfit can be used other orders polys
(ar,br)=polyfit(t,m_nord,1)
xr=polyval([ar,br],t)


###signifitanztest
m1=mean(m_nord)

print 'Nord 80-88 deg:'
print 'relative mp fraction: '
print 'Mittelwert1: '+str(m1)
#Linear regression using stats.linregress
(a1,b1,r1,tt1,s1)=linregress(t,xr)
print 'Steigung: '+str(a1)
print 'Stdabw1: '+str(s1)
print 'Korrelationscoeff1: '+str(r1)
print 'signifikanztest1: '+str(tt1)



#####relative plot
#figure()
#ylabel('melt pond fraction')
#title('mean MODIS melt pond fraction - 80-88 deg N - zeitreihe_4.py')
##ylim([0,0.2])
#plot(dates_short,m_nord,'ro-',label='rel melt pond frac NORD')

#plot(dates_short,xr,'r--',label='trend')
#legend(loc=3)
#grid()
#savefig('/pf/u/u241127/plots/meltponds_timeline2000_2011_relNORD.png')


###zeitreihe 60-87 deg N
#m_nord_allfillup=mean_all[2:10,:5].flatten()


m2000=nanmean(y2000.flatten())
m2001=nanmean(y2001.flatten())
m2002=nanmean(y2002.flatten())
m2003=nanmean(y2003.flatten())
m2004=nanmean(y2004.flatten())
m2005=nanmean(y2005.flatten())
m2006=nanmean(y2006.flatten())
m2007=nanmean(y2007.flatten())
m2008=nanmean(y2008.flatten())
m2009=nanmean(y2009.flatten())
m2010=nanmean(y2010.flatten())
m2011=nanmean(y2011.flatten())

m=array([m2000,m2001,m2002,m2003,m2004,m2005,m2006,m2007,m2008,m2009,m2010,m2011])

d_short=['2000-08-04','2001-08-04','2002-08-04','2003-07-03','2004-07-11','2005-06-25','2006-07-27','2007-06-09','2008-07-03','2009-08-04','2010-06-25','2011-06-25',]  #12

dd_short=[]
for i in range(12):
	d=dt.datetime.strptime(d_short[i], "%Y-%m-%d")
	dd_short.append(d)

dates_short=array(dd_short)


############Trendberechnung relativ mp fraction all
t=arange(12)


#Linear regressison -polyfit - polyfit can be used other orders polys
(ar,br)=polyfit(t,m,1)
xr=polyval([ar,br],t)


###signifitanztest
m1=mean(m)
print 'relative mp fraction1: '
#Linear regression using stats.linregress

print 'Mittelwert: '+str(m1)
(a1,b1,r1,tt1,s1)=linregress(t,xr)
print 'Steigung: '+str(a1)
print 'Stdabw1: '+str(s1)
print 'Stdabw2: '+str(std(m))
print 'Korrelationscoeff1: '+str(r1)
print 'signifikanztest1: '+str(tt1)

print '----'

######relative
#figure()
#ylabel('melt pond fraction')
#title('mean MODIS melt pond fraction - 60-88 deg N - zeitreihe_4.py')
##ylim([0,0.2])
#plot(dates_short,m,'ro-',label='rel melt pond frac NORD')

#plot(dates_short,xr,'r--',label='trend')
#legend(loc=3)
#grid()
#savefig('/pf/u/u241127/plots/meltponds_timeline2000_2011_rel.png')






########################################################

#
#run zeitreihe3.py

##################flaechenskalierung
#NSIDC mean IC values (MONTHLY)
a_ice05=array(loadtxt(homefolder+'Arctic_sea_ice_area_May.txt',skiprows=1, dtype='S')[-12:,-1],dtype=float)
a_ice06=array(loadtxt(homefolder+'Arctic_sea_ice_area_Juni.txt',skiprows=1, dtype='S')[-12:,-1],dtype=float)
a_ice07=array(loadtxt(homefolder+'Arctic_sea_ice_area_July.txt',skiprows=1, dtype='S')[-12:,-1],dtype=float)
a_ice08=array(loadtxt(homefolder+'Arctic_sea_ice_area_Aug.txt',skiprows=1, dtype='S')[-12:,-1],dtype=float)
a_ice09=array(loadtxt(homefolder+'Arctic_sea_ice_area_Sep.txt',skiprows=1, dtype='S')[-12:,-1],dtype=float)

ic2000=a_ice05[0],a_ice06[0],a_ice07[0],a_ice08[0],a_ice09[0]
ic2001=a_ice05[1],a_ice06[1],a_ice07[1],a_ice08[1],a_ice09[1]
ic2002=a_ice05[2],a_ice06[2],a_ice07[2],a_ice08[3],a_ice09[2]
ic2003=a_ice05[3],a_ice06[3],a_ice07[3],a_ice08[3],a_ice09[3]
ic2004=a_ice05[4],a_ice06[4],a_ice07[4],a_ice08[4],a_ice09[4]
ic2005=a_ice05[5],a_ice06[5],a_ice07[5],a_ice08[5],a_ice09[5]
ic2006=a_ice05[6],a_ice06[6],a_ice07[6],a_ice08[6],a_ice09[6]
ic2007=a_ice05[7],a_ice06[7],a_ice07[7],a_ice08[7],a_ice09[7]
ic2008=a_ice05[8],a_ice06[8],a_ice07[8],a_ice08[8],a_ice09[8]
ic2009=a_ice05[9],a_ice06[9],a_ice07[9],a_ice08[9],a_ice09[9]
ic2010=a_ice05[10],a_ice06[10],a_ice07[10],a_ice08[10],a_ice09[10]
ic2011=a_ice05[11],a_ice05[11],a_ice07[11],a_ice08[11],a_ice09[11]

####interpolate to get weekly values
ic2000=ndimage.zoom(ic2000,3.2)
ic2001=ndimage.zoom(ic2001,3.2)
ic2002=ndimage.zoom(ic2002,3.2)
ic2003=ndimage.zoom(ic2003,3.2)
ic2004=ndimage.zoom(ic2004,3.2)
ic2005=ndimage.zoom(ic2005,3.2)
ic2006=ndimage.zoom(ic2006,3.2)
ic2007=ndimage.zoom(ic2007,3.2)
ic2008=ndimage.zoom(ic2008,3.2)
ic2009=ndimage.zoom(ic2009,3.2)
ic2010=ndimage.zoom(ic2010,3.2)
ic2011=ndimage.zoom(ic2011,3.2)

mean_ic=(ic2000+ic2001+ic2002+ic2003+ic2004+ic2005+ic2006+ic2007+ic2008+ic2009+ic2010+ic2011)/12.


#wochenmittel
a_ic_129=mean(array([ic2000[0],ic2001[0],ic2002[0],ic2003[0],ic2004[0],ic2005[0],ic2006[0],ic2007[0],ic2008[0],ic2009[0],ic2010[0],ic2011[0]]))
a_ic_137=mean(array([ic2000[1],ic2001[1],ic2002[1],ic2003[1],ic2004[1],ic2005[1],ic2006[1],ic2007[1],ic2008[1],ic2009[1],ic2010[1],ic2011[1]]))
a_ic_145=mean(array([ic2000[2],ic2001[2],ic2002[2],ic2003[2],ic2004[2],ic2005[2],ic2006[2],ic2007[2],ic2008[2],ic2009[2],ic2010[2],ic2011[2]]))
a_ic_153=mean(array([ic2000[3],ic2001[3],ic2002[3],ic2003[3],ic2004[3],ic2005[3],ic2006[3],ic2007[3],ic2008[3],ic2009[3],ic2010[3],ic2011[3]]))
a_ic_161=mean(array([ic2000[4],ic2001[4],ic2002[4],ic2003[4],ic2004[4],ic2005[4],ic2006[4],ic2007[4],ic2008[4],ic2009[4],ic2010[4],ic2011[4]]))
a_ic_169=mean(array([ic2000[5],ic2001[5],ic2002[5],ic2003[5],ic2004[5],ic2005[5],ic2006[5],ic2007[5],ic2008[5],ic2009[5],ic2010[5],ic2011[5]]))
a_ic_177=mean(array([ic2000[6],ic2001[6],ic2002[6],ic2003[6],ic2004[6],ic2005[6],ic2006[6],ic2007[6],ic2008[6],ic2009[6],ic2010[6],ic2011[6]]))
a_ic_185=mean(array([ic2000[7],ic2001[7],ic2002[7],ic2003[7],ic2004[7],ic2005[7],ic2006[7],ic2007[7],ic2008[7],ic2009[7],ic2010[7],ic2011[7]]))
a_ic_193=mean(array([ic2000[8],ic2001[8],ic2002[8],ic2003[8],ic2004[8],ic2005[8],ic2006[8],ic2007[8],ic2008[8],ic2009[8],ic2010[8],ic2011[8]]))
a_ic_201=mean(array([ic2000[9],ic2001[9],ic2002[9],ic2003[9],ic2004[9],ic2005[9],ic2006[9],ic2007[9],ic2008[9],ic2009[9],ic2010[9],ic2011[9]]))
a_ic_209=mean(array([ic2000[10],ic2001[10],ic2002[10],ic2003[10],ic2004[10],ic2005[10],ic2006[10],ic2007[10],ic2008[10],ic2009[10],ic2010[10],ic2011[10]]))
a_ic_217=mean(array([ic2000[11],ic2001[11],ic2002[11],ic2003[11],ic2004[11],ic2005[11],ic2006[11],ic2007[11],ic2008[11],ic2009[11],ic2010[11],ic2011[11]]))
a_ic_225=mean(array([ic2000[12],ic2001[12],ic2002[12],ic2003[12],ic2004[12],ic2005[12],ic2006[12],ic2007[12],ic2008[12],ic2009[12],ic2010[12],ic2011[12]]))
a_ic_233=mean(array([ic2000[13],ic2001[13],ic2002[13],ic2003[13],ic2004[13],ic2005[13],ic2006[13],ic2007[13],ic2008[13],ic2009[13],ic2010[13],ic2011[13]]))
a_ic_241=mean(array([ic2000[14],ic2001[14],ic2002[14],ic2003[14],ic2004[14],ic2005[14],ic2006[14],ic2007[14],ic2008[14],ic2009[14],ic2010[14],ic2011[14]]))
a_ic_249=mean(array([ic2000[15],ic2001[15],ic2002[15],ic2003[15],ic2004[15],ic2005[15],ic2006[15],ic2007[15],ic2008[15],ic2009[15],ic2010[15],ic2011[15]]))

mean_ic_weekly=array([a_ic_129,a_ic_137,a_ic_145,a_ic_153,a_ic_161,a_ic_169,a_ic_177,a_ic_185,a_ic_193,a_ic_201,a_ic_209,a_ic_217,a_ic_225,a_ic_233,a_ic_241,a_ic_249])


#wochenminimum
a_ic_129min=min(array([ic2000[0],ic2001[0],ic2002[0],ic2003[0],ic2004[0],ic2005[0],ic2006[0],ic2007[0],ic2008[0],ic2009[0],ic2010[0],ic2011[0]]))
a_ic_137min=min(array([ic2000[1],ic2001[1],ic2002[1],ic2003[1],ic2004[1],ic2005[1],ic2006[1],ic2007[1],ic2008[1],ic2009[1],ic2010[1],ic2011[1]]))
a_ic_145min=min(array([ic2000[2],ic2001[2],ic2002[2],ic2003[2],ic2004[2],ic2005[2],ic2006[2],ic2007[2],ic2008[2],ic2009[2],ic2010[2],ic2011[2]]))
a_ic_153min=min(array([ic2000[3],ic2001[3],ic2002[3],ic2003[3],ic2004[3],ic2005[3],ic2006[3],ic2007[3],ic2008[3],ic2009[3],ic2010[3],ic2011[3]]))
a_ic_161min=min(array([ic2000[4],ic2001[4],ic2002[4],ic2003[4],ic2004[4],ic2005[4],ic2006[4],ic2007[4],ic2008[4],ic2009[4],ic2010[4],ic2011[4]]))
a_ic_169min=min(array([ic2000[5],ic2001[5],ic2002[5],ic2003[5],ic2004[5],ic2005[5],ic2006[5],ic2007[5],ic2008[5],ic2009[5],ic2010[5],ic2011[5]]))
a_ic_177min=min(array([ic2000[6],ic2001[6],ic2002[6],ic2003[6],ic2004[6],ic2005[6],ic2006[6],ic2007[6],ic2008[6],ic2009[6],ic2010[6],ic2011[6]]))
a_ic_185min=min(array([ic2000[7],ic2001[7],ic2002[7],ic2003[7],ic2004[7],ic2005[7],ic2006[7],ic2007[7],ic2008[7],ic2009[7],ic2010[7],ic2011[7]]))
a_ic_193min=min(array([ic2000[8],ic2001[8],ic2002[8],ic2003[8],ic2004[8],ic2005[8],ic2006[8],ic2007[8],ic2008[8],ic2009[8],ic2010[8],ic2011[8]]))
a_ic_201min=min(array([ic2000[9],ic2001[9],ic2002[9],ic2003[9],ic2004[9],ic2005[9],ic2006[9],ic2007[9],ic2008[9],ic2009[9],ic2010[9],ic2011[9]]))
a_ic_209min=min(array([ic2000[10],ic2001[10],ic2002[10],ic2003[10],ic2004[10],ic2005[10],ic2006[10],ic2007[10],ic2008[10],ic2009[10],ic2010[10],ic2011[10]]))
a_ic_217min=min(array([ic2000[11],ic2001[11],ic2002[11],ic2003[11],ic2004[11],ic2005[11],ic2006[11],ic2007[11],ic2008[11],ic2009[11],ic2010[11],ic2011[11]]))
a_ic_225min=min(array([ic2000[12],ic2001[12],ic2002[12],ic2003[12],ic2004[12],ic2005[12],ic2006[12],ic2007[12],ic2008[12],ic2009[12],ic2010[12],ic2011[12]]))
a_ic_233min=min(array([ic2000[13],ic2001[13],ic2002[13],ic2003[13],ic2004[13],ic2005[13],ic2006[13],ic2007[13],ic2008[13],ic2009[13],ic2010[13],ic2011[13]]))
a_ic_241min=min(array([ic2000[14],ic2001[14],ic2002[14],ic2003[14],ic2004[14],ic2005[14],ic2006[14],ic2007[14],ic2008[14],ic2009[14],ic2010[14],ic2011[14]]))
a_ic_249min=min(array([ic2000[15],ic2001[15],ic2002[15],ic2003[15],ic2004[15],ic2005[15],ic2006[15],ic2007[15],ic2008[15],ic2009[15],ic2010[15],ic2011[15]]))

min_ic_weekly=array([a_ic_129min,a_ic_137min,a_ic_145min,a_ic_153min,a_ic_161min,a_ic_169min,a_ic_177min,a_ic_185min,a_ic_193min,a_ic_201min,a_ic_209min,a_ic_217min,a_ic_225min,a_ic_233min,a_ic_241min,a_ic_249min])





mp2000=array([nanmean(y2000[:,0].flatten()), 
nanmean(y2000[:,1].flatten()), nanmean(y2000[:,2].flatten()),  nanmean(y2000[:,3].flatten()), nanmean(y2000[:,4].flatten()), nanmean(y2000[:,5].flatten()), nanmean(y2000[:,6].flatten()), nanmean(y2000[:,7].flatten()), nanmean(y2000[:,8].flatten()), nanmean(y2000[:,9].flatten()), nanmean(y2000[:,10].flatten()), nanmean(y2000[:,11].flatten()), nanmean(y2000[:,12].flatten()), nanmean(y2000[:,13].flatten()), nanmean(y2000[:,14].flatten()), nanmean(y2000[:,15].flatten())])
max2000=nanmax(mp2000)



#mp2001=array([nanmean(y2001[:,0].flatten()), nanmean(y2001[:,1].flatten()),nanmean(y2001[:,2].flatten()), nanmean(y2001[:,3].flatten()), nanmean(y2001[:,4].flatten()), nanmean(y2001[:,5].flatten()), nanmean(y2001[:,6].flatten()), nanmean(y2001[:,7].flatten()), nanmean(y2001[:,8].flatten()), nanmean(y2001[:,9].flatten()), nanmean(y2001[:,10].flatten()), nanmean(y2001[:,11].flatten()), nanmean(y2001[:,12].flatten()), nanmean(y2001[:,13].flatten()), nanmean(y2001[:,14].flatten()), nanmean(y2001[:,15].flatten())])
mp2001=array([ 0.13500026,  0.16612339,  0.18744891,  0.21322589,  0.26326811,
               0.28,         0.29,  0.30378151,  0.30221912,  0.31751576,
        0.31650484,  0.30091347,  0.29028418,  0.27302649,  0.27476758,
        0.25299009])

max2001=nanmax(mp2001)


mp2002=array([nanmean(y2002[:,0].flatten()), nanmean(y2002[:,1].flatten()),nanmean(y2002[:,2].flatten()), nanmean(y2002[:,3].flatten()), nanmean(y2002[:,4].flatten()), nanmean(y2002[:,5].flatten()), nanmean(y2002[:,6].flatten()), nanmean(y2002[:,7].flatten()), nanmean(y2002[:,8].flatten()), nanmean(y2002[:,9].flatten()), nanmean(y2002[:,10].flatten()), nanmean(y2002[:,11].flatten()), nanmean(y2002[:,12].flatten()), nanmean(y2002[:,13].flatten()), nanmean(y2002[:,14].flatten()), nanmean(y2002[:,15].flatten())])
max2002=nanmax(mp2002)

mp2003=array([nanmean(y2003[:,0].flatten()), nanmean(y2003[:,1].flatten()),nanmean(y2003[:,2].flatten()), nanmean(y2003[:,3].flatten()), nanmean(y2003[:,4].flatten()), nanmean(y2003[:,5].flatten()), nanmean(y2003[:,6].flatten()), nanmean(y2003[:,7].flatten()), nanmean(y2003[:,8].flatten()), nanmean(y2003[:,9].flatten()), nanmean(y2003[:,10].flatten()), nanmean(y2003[:,11].flatten()), nanmean(y2003[:,12].flatten()), nanmean(y2003[:,13].flatten()), nanmean(y2003[:,14].flatten()), nanmean(y2003[:,15].flatten())])
max2003=nanmax(mp2003)

mp2004=array([nanmean(y2004[:,0].flatten()), nanmean(y2004[:,1].flatten()),nanmean(y2004[:,2].flatten()), nanmean(y2004[:,3].flatten()), nanmean(y2004[:,4].flatten()), nanmean(y2004[:,5].flatten()), nanmean(y2004[:,6].flatten()), nanmean(y2004[:,7].flatten()), nanmean(y2004[:,8].flatten()), nanmean(y2004[:,9].flatten()), nanmean(y2004[:,10].flatten()), nanmean(y2004[:,11].flatten()), nanmean(y2004[:,12].flatten()), nanmean(y2004[:,13].flatten()), nanmean(y2004[:,14].flatten()), nanmean(y2004[:,15].flatten())])
max2004=nanmax(mp2004)

mp2005=array([nanmean(y2005[:,0].flatten()), nanmean(y2005[:,1].flatten()),nanmean(y2005[:,2].flatten()), nanmean(y2005[:,3].flatten()), nanmean(y2005[:,4].flatten()), nanmean(y2005[:,5].flatten()), nanmean(y2005[:,6].flatten()), nanmean(y2005[:,7].flatten()), nanmean(y2005[:,8].flatten()), nanmean(y2005[:,9].flatten()), nanmean(y2005[:,10].flatten()), nanmean(y2005[:,11].flatten()), nanmean(y2005[:,12].flatten()), nanmean(y2005[:,13].flatten()), nanmean(y2005[:,14].flatten()), nanmean(y2005[:,15].flatten())])
max2005=nanmax(mp2005)

mp2006=array([nanmean(y2006[:,0].flatten()), nanmean(y2006[:,1].flatten()),nanmean(y2006[:,2].flatten()), nanmean(y2006[:,3].flatten()), nanmean(y2006[:,4].flatten()), nanmean(y2006[:,5].flatten()), nanmean(y2006[:,6].flatten()), nanmean(y2006[:,7].flatten()), nanmean(y2006[:,8].flatten()), nanmean(y2006[:,9].flatten()), nanmean(y2006[:,10].flatten()), nanmean(y2006[:,11].flatten()), nanmean(y2006[:,12].flatten()), nanmean(y2006[:,13].flatten()), nanmean(y2006[:,14].flatten()), nanmean(y2006[:,15].flatten())])
max2006=nanmax(mp2006)

mp2007=array([nanmean(y2007[:,0].flatten()), nanmean(y2007[:,1].flatten()),nanmean(y2007[:,2].flatten()), nanmean(y2007[:,3].flatten()), nanmean(y2007[:,4].flatten()), nanmean(y2007[:,5].flatten()), nanmean(y2007[:,6].flatten()), nanmean(y2007[:,7].flatten()), nanmean(y2007[:,8].flatten()), nanmean(y2007[:,9].flatten()), nanmean(y2007[:,10].flatten()), nanmean(y2007[:,11].flatten()), nanmean(y2007[:,12].flatten()), nanmean(y2007[:,13].flatten()), nanmean(y2007[:,14].flatten()), nanmean(y2007[:,15].flatten())])
max2007=nanmax(mp2007)

mp2008=array([nanmean(y2008[:,0].flatten()), nanmean(y2008[:,1].flatten()),nanmean(y2008[:,2].flatten()), nanmean(y2008[:,3].flatten()), nanmean(y2008[:,4].flatten()), nanmean(y2008[:,5].flatten()), nanmean(y2008[:,6].flatten()), nanmean(y2008[:,7].flatten()), nanmean(y2008[:,8].flatten()), nanmean(y2008[:,9].flatten()), nanmean(y2008[:,10].flatten()), nanmean(y2008[:,11].flatten()), nanmean(y2008[:,12].flatten()), nanmean(y2008[:,13].flatten()), nanmean(y2008[:,14].flatten()), nanmean(y2008[:,15].flatten())])
max2008=nanmax(mp2008)

mp2009=array([nanmean(y2009[:,0].flatten()), nanmean(y2009[:,1].flatten()),nanmean(y2009[:,2].flatten()), nanmean(y2009[:,3].flatten()), nanmean(y2009[:,4].flatten()), nanmean(y2009[:,5].flatten()), nanmean(y2009[:,6].flatten()), nanmean(y2009[:,7].flatten()), nanmean(y2009[:,8].flatten()), nanmean(y2009[:,9].flatten()), nanmean(y2009[:,10].flatten()), nanmean(y2009[:,11].flatten()), nanmean(y2009[:,12].flatten()), nanmean(y2009[:,13].flatten()), nanmean(y2009[:,14].flatten()), nanmean(y2009[:,15].flatten())])
max2009=nanmax(mp2009)

mp2010=array([nanmean(y2010[:,0].flatten()), nanmean(y2010[:,1].flatten()),nanmean(y2010[:,2].flatten()), nanmean(y2010[:,3].flatten()), nanmean(y2010[:,4].flatten()), nanmean(y2010[:,5].flatten()), nanmean(y2010[:,6].flatten()), nanmean(y2010[:,7].flatten()), nanmean(y2010[:,8].flatten()), nanmean(y2010[:,9].flatten()), nanmean(y2010[:,10].flatten()), nanmean(y2010[:,11].flatten()), nanmean(y2010[:,12].flatten()), nanmean(y2010[:,13].flatten()), nanmean(y2010[:,14].flatten()), nanmean(y2010[:,15].flatten())])
max2010=nanmax(mp2010)

mp2011=array([nanmean(y2011[:,0].flatten()), nanmean(y2011[:,1].flatten()),nanmean(y2011[:,2].flatten()), nanmean(y2011[:,3].flatten()), nanmean(y2011[:,4].flatten()), nanmean(y2011[:,5].flatten()), nanmean(y2011[:,6].flatten()), nanmean(y2011[:,7].flatten()), nanmean(y2011[:,8].flatten()), nanmean(y2011[:,9].flatten()), nanmean(y2011[:,10].flatten()), nanmean(y2011[:,11].flatten()), nanmean(y2011[:,12].flatten()), nanmean(y2011[:,13].flatten()), nanmean(y2011[:,14].flatten()), nanmean(y2011[:,15].flatten())])
max2011=nanmax(mp2011)




mpmax=array([max2000,max2001,max2002,max2003,max2004,max2005,max2006,max2007,max2008,max2009,max2010,max2011])
meanyear=array([mean(mp2000), mean(mp2001), mean(mp2002), mean(mp2003), mean(mp2004), mean(mp2005), mean(mp2006), mean(mp2007), mean(mp2008), mean(mp2009),  mean(mp2010), mean(mp2011)])
meanyear_std=array([scipy.stats.nanstd(mp2000), scipy.stats.nanstd(mp2001), scipy.stats.nanstd(mp2002), scipy.stats.nanstd(mp2003), scipy.stats.nanstd(mp2004), scipy.stats.nanstd(mp2005), scipy.stats.nanstd(mp2006), scipy.stats.nanstd(mp2007), scipy.stats.nanstd(mp2008), scipy.stats.nanstd(mp2009),  scipy.stats.nanstd(mp2010), scipy.stats.nanstd(mp2011)])


mean129=array([nanmean(y2000[:,0].flatten()), nanmean(y2001[:,0].flatten()), nanmean(y2002[:,0].flatten()), nanmean(y2003[:,0].flatten()), nanmean(y2004[:,0].flatten()), nanmean(y2005[:,0].flatten()), nanmean(y2006[:,0].flatten()), nanmean(y2007[:,0].flatten()), nanmean(y2008[:,0].flatten()), nanmean(y2009[:,0].flatten()), nanmean(y2010[:,0].flatten()), nanmean(y2011[:,0].flatten())])

mean137=array([nanmean(y2000[:,1].flatten()), nanmean(y2001[:,1].flatten()), nanmean(y2002[:,1].flatten()), nanmean(y2003[:,1].flatten()), nanmean(y2004[:,1].flatten()), nanmean(y2005[:,1].flatten()), nanmean(y2006[:,1].flatten()), nanmean(y2007[:,1].flatten()), nanmean(y2008[:,1].flatten()), nanmean(y2009[:,1].flatten()), nanmean(y2010[:,1].flatten()), nanmean(y2011[:,1].flatten())])

mean145=array([nanmean(y2000[:,2].flatten()), nanmean(y2001[:,2].flatten()), nanmean(y2002[:,2].flatten()), nanmean(y2003[:,2].flatten()), nanmean(y2004[:,2].flatten()), nanmean(y2005[:,2].flatten()), nanmean(y2006[:,2].flatten()), nanmean(y2007[:,2].flatten()), nanmean(y2008[:,2].flatten()), nanmean(y2009[:,2].flatten()), nanmean(y2010[:,2].flatten()), nanmean(y2011[:,2].flatten())])

mean153=array([nanmean(y2000[:,3].flatten()), nanmean(y2001[:,3].flatten()), nanmean(y2002[:,3].flatten()), nanmean(y2003[:,3].flatten()), nanmean(y2004[:,3].flatten()), nanmean(y2005[:,3].flatten()), nanmean(y2006[:,3].flatten()), nanmean(y2007[:,3].flatten()), nanmean(y2008[:,3].flatten()), nanmean(y2009[:,3].flatten()), nanmean(y2010[:,3].flatten()), nanmean(y2011[:,3].flatten())])

mean161=array([nanmean(y2000[:,4].flatten()), nanmean(y2001[:,4].flatten()), nanmean(y2002[:,4].flatten()), nanmean(y2003[:,4].flatten()), nanmean(y2004[:,4].flatten()), nanmean(y2005[:,4].flatten()), nanmean(y2006[:,4].flatten()), nanmean(y2007[:,4].flatten()), nanmean(y2008[:,4].flatten()), nanmean(y2009[:,4].flatten()), nanmean(y2010[:,4].flatten()), nanmean(y2011[:,4].flatten())])

#mean169=array([nanmean(y2000[:,5].flatten()), nanmean(y2001[:,5].flatten()), nanmean(y2002[:,5].flatten()), nanmean(y2003[:,5].flatten()), nanmean(y2004[:,5].flatten()), nanmean(y2005[:,5].flatten()), nanmean(y2006[:,5].flatten()), nanmean(y2007[:,5].flatten()), nanmean(y2008[:,5].flatten()), nanmean(y2009[:,5].flatten()), nanmean(y2010[:,5].flatten()), nanmean(y2011[:,5].flatten())])


#mean177=array([nanmean(y2000[:,6].flatten()), nanmean(y2001[:,6].flatten()), nanmean(y2002[:,6].flatten()), nanmean(y2003[:,6].flatten()), nanmean(y2004[:,6].flatten()), nanmean(y2005[:,6].flatten()), nanmean(y2006[:,6].flatten()), nanmean(y2007[:,6].flatten()), nanmean(y2008[:,6].flatten()), nanmean(y2009[:,6].flatten()), nanmean(y2010[:,6].flatten()), nanmean(y2011[:,6].flatten())])

mean169=array([ 0.25751303,         0.24,  0.26249737,  0.25012233,  0.23435746,
        0.26901872,  0.26327172,  0.28659674,  0.26187282,  0.25292194,
        0.28432131,  0.28099843])


mean177=array([ 0.32738856,         0.28,  0.26626512,  0.28503644,  0.26915024,
        0.28678615,  0.27009135,  0.30431904,  0.28645887,  0.27204823,
        0.31221414,  0.31334884])


mean185=array([nanmean(y2000[:,7].flatten()), nanmean(y2001[:,7].flatten()), nanmean(y2002[:,7].flatten()), nanmean(y2003[:,7].flatten()), nanmean(y2004[:,7].flatten()), nanmean(y2005[:,7].flatten()), nanmean(y2006[:,7].flatten()), nanmean(y2007[:,7].flatten()), nanmean(y2008[:,7].flatten()), nanmean(y2009[:,7].flatten()), nanmean(y2010[:,7].flatten()), nanmean(y2011[:,7].flatten())])

mean193=array([nanmean(y2000[:,8].flatten()), nanmean(y2001[:,8].flatten()), nanmean(y2002[:,8].flatten()), nanmean(y2003[:,8].flatten()), nanmean(y2004[:,8].flatten()), nanmean(y2005[:,8].flatten()), nanmean(y2006[:,8].flatten()), nanmean(y2007[:,8].flatten()), nanmean(y2008[:,8].flatten()), nanmean(y2009[:,8].flatten()), nanmean(y2010[:,8].flatten()), nanmean(y2011[:,8].flatten())])

mean201=array([nanmean(y2000[:,9].flatten()), nanmean(y2001[:,9].flatten()), nanmean(y2002[:,9].flatten()), nanmean(y2003[:,9].flatten()), nanmean(y2004[:,9].flatten()), nanmean(y2005[:,9].flatten()), nanmean(y2006[:,9].flatten()), nanmean(y2007[:,9].flatten()), nanmean(y2008[:,9].flatten()), nanmean(y2009[:,9].flatten()), nanmean(y2010[:,9].flatten()), nanmean(y2011[:,9].flatten())])

mean209=array([nanmean(y2000[:,10].flatten()), nanmean(y2001[:,10].flatten()), nanmean(y2002[:,10].flatten()), nanmean(y2003[:,10].flatten()), nanmean(y2004[:,10].flatten()), nanmean(y2005[:,10].flatten()), nanmean(y2006[:,10].flatten()), nanmean(y2007[:,10].flatten()), nanmean(y2008[:,10].flatten()), nanmean(y2009[:,10].flatten()), nanmean(y2010[:,10].flatten()), nanmean(y2011[:,10].flatten())])

mean217=array([nanmean(y2000[:,11].flatten()), nanmean(y2001[:,11].flatten()), nanmean(y2002[:,11].flatten()), nanmean(y2003[:,11].flatten()), nanmean(y2004[:,11].flatten()), nanmean(y2005[:,11].flatten()), nanmean(y2006[:,11].flatten()), nanmean(y2007[:,11].flatten()), nanmean(y2008[:,11].flatten()), nanmean(y2009[:,11].flatten()), nanmean(y2010[:,11].flatten()), nanmean(y2011[:,11].flatten())])

mean225=array([nanmean(y2000[:,12].flatten()), nanmean(y2001[:,12].flatten()), nanmean(y2002[:,12].flatten()), nanmean(y2003[:,12].flatten()), nanmean(y2004[:,12].flatten()), nanmean(y2005[:,12].flatten()), nanmean(y2006[:,12].flatten()), nanmean(y2007[:,12].flatten()), nanmean(y2008[:,12].flatten()), nanmean(y2009[:,12].flatten()), nanmean(y2010[:,12].flatten()), nanmean(y2011[:,12].flatten())])

mean233=array([nanmean(y2000[:,13].flatten()), nanmean(y2001[:,13].flatten()), nanmean(y2002[:,13].flatten()), nanmean(y2003[:,13].flatten()), nanmean(y2004[:,13].flatten()), nanmean(y2005[:,13].flatten()), nanmean(y2006[:,13].flatten()), nanmean(y2007[:,13].flatten()), nanmean(y2008[:,13].flatten()), nanmean(y2009[:,13].flatten()), nanmean(y2010[:,13].flatten()), nanmean(y2011[:,13].flatten())])

mean241=array([nanmean(y2000[:,14].flatten()), nanmean(y2001[:,14].flatten()), nanmean(y2002[:,14].flatten()), nanmean(y2003[:,14].flatten()), nanmean(y2004[:,14].flatten()), nanmean(y2005[:,14].flatten()), nanmean(y2006[:,14].flatten()), nanmean(y2007[:,14].flatten()), nanmean(y2008[:,14].flatten()), nanmean(y2009[:,14].flatten()), nanmean(y2010[:,14].flatten()), nanmean(y2011[:,14].flatten())])

mean249=array([nanmean(y2000[:,15].flatten()), nanmean(y2001[:,15].flatten()), nanmean(y2002[:,15].flatten()), nanmean(y2003[:,15].flatten()), nanmean(y2004[:,15].flatten()), nanmean(y2005[:,15].flatten()), nanmean(y2006[:,15].flatten()), nanmean(y2007[:,15].flatten()), nanmean(y2008[:,15].flatten()), nanmean(y2009[:,15].flatten()), nanmean(y2010[:,15].flatten()), nanmean(y2011[:,15].flatten())])





meanweek=array([mean(mean129), mean(mean137), mean(mean145), mean(mean153), mean(mean161), mean(mean169), mean(mean177), mean(mean185), mean(mean193), mean(mean201), mean(mean209), mean(mean217), mean(mean225), mean(mean233), mean(mean241), mean(mean249)])

stddevweek=array([std(mean129), std(mean137), std(mean145), std(mean153), std(mean161), std(mean169), std(mean177), std(mean185), std(mean193), std(mean201), std(mean209), std(mean217), std(mean225), std(mean233), std(mean241), std(mean249)])

stdplus=meanweek+stddevweek
stdminus=meanweek-stddevweek




homefolder='/pf/u/u241127/data/'
mt=loadtxt(homefolder+'meantable2masked.csv',skiprows=1,dtype=float32,delimiter=',')
mt[mt==9.9999]=nan
years=mt[0,1:-2]
days=mt[:-2,0][1:]

params = {'axes.labelsize': 16,
          'text.fontsize': 16,
          'legend.fontsize': 16,
          'xtick.labelsize': 16,
          'ytick.labelsize': 16,}
rcParams.update(params)

###########plot 2007 and 2011
fig = figure()
ax = fig.add_subplot(111)
ax.plot_date(days,meanweek,'k-',label='mean 2000-2011',linewidth=4)
ax.plot_date(days,stdplus,'k--',linewidth=1)
ax.plot_date(days,stdminus,'k--',label='standard deviation',linewidth=1)
ax.plot_date(days,mp2007,'r-',label='2007',linewidth=2)
ax.plot_date(days,mp2011,'m-',label='2011',linewidth=2)
ax.plot_date(days,mp2000,'-',color='Gray',label='other years',linewidth=0.5)
ax.plot_date(days,mp2001,'-',color='Gray',linewidth=0.5)
ax.plot_date(days,mp2002,'-',color='Gray',linewidth=0.5)
ax.plot_date(days,mp2003,'-',color='Gray',linewidth=0.5)
ax.plot_date(days,mp2004,'-',color='Gray',linewidth=0.5)
ax.plot_date(days,mp2005,'-',color='Gray',linewidth=0.5)
ax.plot_date(days,mp2006,'-',color='Gray',linewidth=0.5)
ax.plot_date(days,mp2007,'-',color='Gray',linewidth=0.5)
ax.plot_date(days,mp2008,'-',color='Gray',linewidth=0.5)
ax.plot_date(days,mp2009,'-',color='Gray',linewidth=0.5)
ax.plot_date(days,mp2010,'-',color='Gray',linewidth=0.5)
#xlabel('date')
ylabel('Melt pond fraction')
legend(loc=4)
grid()
#title('melt pond fraction 2007 and 2011')
#ax.set_xlim(days2[0], days2[-1] )
ax.xaxis.set_major_formatter( DateFormatter('%m.%d.') )
#ax.set_xlabel(fontsize=14)

#title('relative melt pond fraction - zeitreihe_4.py')
savefig('/pf/u/u241127/plots/zeitr_means_20072011_masked.png')

###########plot rel means b&w for the cryosphere paper
fig = figure()
ax = fig.add_subplot(111)
ax.plot_date(days,meanweek,'k-',label='mean',linewidth=4)
ax.plot_date(days,stdplus,'k--',linewidth=1)
ax.plot_date(days,stdminus,'k--',label='standard deviation',linewidth=1)
ax.plot_date(days,mp2007,'-',color='Gray',linewidth=0.5)
ax.plot_date(days,mp2011,'-',color='Gray',linewidth=0.5)
ax.plot_date(days,mp2000,'-',color='Gray',label='melt ponds 2000-2011',linewidth=0.5)
ax.plot_date(days,mp2001,'-',color='Gray',linewidth=0.5)
ax.plot_date(days,mp2002,'-',color='Gray',linewidth=0.5)
ax.plot_date(days,mp2003,'-',color='Gray',linewidth=0.5)
ax.plot_date(days,mp2004,'-',color='Gray',linewidth=0.5)
ax.plot_date(days,mp2005,'-',color='Gray',linewidth=0.5)
ax.plot_date(days,mp2006,'-',color='Gray',linewidth=0.5)
ax.plot_date(days,mp2007,'-',color='Gray',linewidth=0.5)
ax.plot_date(days,mp2008,'-',color='Gray',linewidth=0.5)
ax.plot_date(days,mp2009,'-',color='Gray',linewidth=0.5)
ax.plot_date(days,mp2010,'-',color='Gray',linewidth=0.5)
#xlabel('date')
ylabel('relative melt pond fraction')
legend(loc=8)
grid()
#title('melt pond fraction 2007 and 2011')
#ax.set_xlim(days2[0], days2[-1] )
ax.xaxis.set_major_formatter( DateFormatter('%d.%m') )
savefig('/pf/u/u241127/plots/zeitr_relmeans_bw.png')



###################flaechenskalierung
##NSIDC mean IC values (MONTHLY)
#a_ice05=array(loadtxt(homefolder+'Arctic_sea_ice_area_May.txt',skiprows=1, dtype='S')[-12:,-1],dtype=float)
#a_ice06=array(loadtxt(homefolder+'Arctic_sea_ice_area_Juni.txt',skiprows=1, dtype='S')[-12:,-1],dtype=float)
#a_ice07=array(loadtxt(homefolder+'Arctic_sea_ice_area_July.txt',skiprows=1, dtype='S')[-12:,-1],dtype=float)
#a_ice08=array(loadtxt(homefolder+'Arctic_sea_ice_area_Aug.txt',skiprows=1, dtype='S')[-12:,-1],dtype=float)
#a_ice09=array(loadtxt(homefolder+'Arctic_sea_ice_area_Sep.txt',skiprows=1, dtype='S')[-12:,-1],dtype=float)




#ic2000=a_ice05[0],a_ice06[0],a_ice07[0],a_ice08[0],a_ice09[0]
#ic2001=a_ice05[1],a_ice06[1],a_ice07[1],a_ice08[1],a_ice09[1]
#ic2002=a_ice05[2],a_ice06[2],a_ice07[2],a_ice08[3],a_ice09[2]
#ic2003=a_ice05[3],a_ice06[3],a_ice07[3],a_ice08[3],a_ice09[3]
#ic2004=a_ice05[4],a_ice06[4],a_ice07[4],a_ice08[4],a_ice09[4]
#ic2005=a_ice05[5],a_ice06[5],a_ice07[5],a_ice08[5],a_ice09[5]
#ic2006=a_ice05[6],a_ice06[6],a_ice07[6],a_ice08[6],a_ice09[6]
#ic2007=a_ice05[7],a_ice06[7],a_ice07[7],a_ice08[7],a_ice09[7]
#ic2008=a_ice05[8],a_ice06[8],a_ice07[8],a_ice08[8],a_ice09[8]
#ic2009=a_ice05[9],a_ice06[9],a_ice07[9],a_ice08[9],a_ice09[9]
#ic2010=a_ice05[10],a_ice06[10],a_ice07[10],a_ice08[10],a_ice09[10]
#ic2011=a_ice05[11],a_ice05[11],a_ice07[11],a_ice08[11],a_ice09[11]

#####interpolate to get weekly values
#ic2000=ndimage.zoom(ic2000,3.2)
#ic2001=ndimage.zoom(ic2001,3.2)
#ic2002=ndimage.zoom(ic2002,3.2)
#ic2003=ndimage.zoom(ic2003,3.2)
#ic2004=ndimage.zoom(ic2004,3.2)
#ic2005=ndimage.zoom(ic2005,3.2)
#ic2006=ndimage.zoom(ic2006,3.2)
#ic2007=ndimage.zoom(ic2007,3.2)
#ic2008=ndimage.zoom(ic2008,3.2)
#ic2009=ndimage.zoom(ic2009,3.2)
#ic2010=ndimage.zoom(ic2010,3.2)
#ic2011=ndimage.zoom(ic2011,3.2)

#mp area in mio km2
a_mp2000=ic2000*mp2000
a_mp2001=ic2001*mp2001
a_mp2002=ic2002*mp2002
a_mp2003=ic2003*mp2003
a_mp2004=ic2004*mp2004
a_mp2005=ic2005*mp2005
a_mp2006=ic2006*mp2006
a_mp2007=ic2007*mp2007
a_mp2008=ic2008*mp2008
a_mp2009=ic2009*mp2009
a_mp2010=ic2010*mp2010
a_mp2011=ic2011*mp2011

###########plot ice conc b&w for the cryosphere paper
fig = figure()
ax = fig.add_subplot(111)
#ax.plot_date(days,meanweek,'k-',label='mean',linewidth=4)
#ax.plot_date(days,stdplus,'k--',linewidth=1)
#ax.plot_date(days,stdminus,'k--',label='standard deviation',linewidth=1)
ax.plot_date(days,ic2007,'-',color='Gray',linewidth=0.5)
ax.plot_date(days,ic2011,'-',color='Gray',linewidth=0.5)
ax.plot_date(days,ic2000,'-',color='Gray',label='sea ice concentrations 2000-2011',linewidth=0.5)
ax.plot_date(days,ic2001,'-',color='Gray',linewidth=0.5)
ax.plot_date(days,ic2002,'-',color='Gray',linewidth=0.5)
ax.plot_date(days,ic2003,'-',color='Gray',linewidth=0.5)
ax.plot_date(days,ic2004,'-',color='Gray',linewidth=0.5)
ax.plot_date(days,ic2005,'-',color='Gray',linewidth=0.5)
ax.plot_date(days,ic2006,'-',color='Gray',linewidth=0.5)
ax.plot_date(days,ic2007,'-',color='Gray',linewidth=0.5)
ax.plot_date(days,ic2008,'-',color='Gray',linewidth=0.5)
ax.plot_date(days,ic2009,'-',color='Gray',linewidth=0.5)
ax.plot_date(days,ic2010,'-',color='Gray',linewidth=0.5)
#xlabel('date')
ylabel('sea ice concentration')
legend(loc=8)
grid()
#ax.set_xlim(days2[0], days2[-1] )
ax.xaxis.set_major_formatter( DateFormatter('%d.%m') )
savefig('/pf/u/u241127/plots/zeitr_iceconc_bw.png')


#STOPPPPP



###a_meanyear
a_meanyear=array([mean(a_mp2000),mean(a_mp2001),mean(a_mp2002),mean(a_mp2003),mean(a_mp2004),mean(a_mp2005),mean(a_mp2006),mean(a_mp2007),mean(a_mp2008),mean(a_mp2009),mean(a_mp2010),mean(a_mp2011)])

a_meanyear_std=array([scipy.stats.nanstd(a_mp2000),scipy.stats.nanstd(a_mp2001), scipy.stats.nanstd(a_mp2002), scipy.stats.nanstd(a_mp2003), scipy.stats.nanstd(a_mp2004), scipy.stats.nanstd(a_mp2005), scipy.stats.nanstd(a_mp2006), scipy.stats.nanstd(a_mp2007), scipy.stats.nanstd(a_mp2008), scipy.stats.nanstd(a_mp2009), scipy.stats.nanstd(a_mp2010), scipy.stats.nanstd(a_mp2011)])



a_maxyear=array([max(a_mp2000),max(a_mp2001),max(a_mp2002),max(a_mp2003),max(a_mp2004),max(a_mp2005),max(a_mp2006),max(a_mp2007),max(a_mp2008),max(a_mp2009),max(a_mp2010),max(a_mp2011)])


a_mp_mean=(a_mp2000+a_mp2001+a_mp2002+a_mp2003+a_mp2004+a_mp2005+a_mp2006+a_mp2007+a_mp2008+a_mp2009+
a_mp2010+a_mp2011)/12.
a_std=sqrt(1./(16.-1.)*((a_mp2000-a_mp_mean)**2+(a_mp2001-a_mp_mean)**2+(a_mp2002-a_mp_mean)**2+(a_mp2003-a_mp_mean)**2+(a_mp2004-a_mp_mean)**2+(a_mp2005-a_mp_mean)**2+(a_mp2006-a_mp_mean)**2+(a_mp2007-a_mp_mean)**2+(a_mp2008-a_mp_mean)**2+(a_mp2009-a_mp_mean)**2+(a_mp2010-a_mp_mean)**2+(a_mp2011-a_mp_mean)**2))
a_stdplus=a_mp_mean+a_std
a_stdminus=a_mp_mean-a_std

fig = figure()
ax = fig.add_subplot(111)
ax.plot_date(days,a_mp_mean,'k-',label='mean 2000-2011',linewidth=4)
ax.plot_date(days,a_stdplus,'k--',linewidth=1)
ax.plot_date(days,a_stdminus,'k--',label='standard deviation',linewidth=1)
ax.plot_date(days,a_mp2007,'r-',label='2007',linewidth=2)
ax.plot_date(days,a_mp2011,'m-',label='2011',linewidth=2)
ax.plot_date(days,a_mp2000,'-',color='gray',linewidth=0.5,label='melt ponds 2000-2011')
ax.plot_date(days,a_mp2001,'-',color='gray',linewidth=0.5)
ax.plot_date(days,a_mp2002,'-',color='gray',linewidth=0.5)
ax.plot_date(days,a_mp2003,'-',color='gray',linewidth=0.5)
ax.plot_date(days,a_mp2004,'-',color='gray',linewidth=0.5)
ax.plot_date(days,a_mp2005,'-',color='gray',linewidth=0.5)
ax.plot_date(days,a_mp2006,'-',color='gray',linewidth=0.5)
ax.plot_date(days,a_mp2007,'-',color='gray',linewidth=0.5)
ax.plot_date(days,a_mp2008,'-',color='gray',linewidth=0.5)
ax.plot_date(days,a_mp2009,'-',color='gray',linewidth=0.5)
ax.plot_date(days,a_mp2010,'-',color='gray',linewidth=0.5)
ax.plot_date(days,a_mp2011,'-',color='gray',linewidth=0.5)
#ax.plot_date(days,meanweek,'k-',label='mean',linewidth=4)
#ax.plot_date(days,stdplus,'k--',linewidth=1)
#ax.plot_date(days,stdminus,'k--',label='standard deviation',linewidth=1)
ax.xaxis.set_major_formatter( DateFormatter('%m.%d.'))
#ax.yaxis(fontsize=16)
#ax.xaxis(fontsize=14)
ylabel('Area in Mio km2')
#legend(loc=1)
grid()
title('weekly absolute melt pond area entire Arctic - zeitreihe_5.py')
savefig('/pf/u/u241127/plots/zeitr_mp_area_masked.png')







###datumsarray bauen
a=days.tolist()
a.append(257.0)
a=12*a


b=years.tolist()
b=17*b
b.sort()


daystr=[]
monthstr=[]
yearstr=[]
for i in range(17*12):
	y,m,d=JulDay2Date(int(b[i]),int(a[i]))
	daystr.append(d)
	monthstr.append(m)
	yearstr.append(y)

datestr=[]
for i in range(17*12):
	datestring=yearstr[i]+'-'+monthstr[i]+'-'+daystr[i]
	datestr.append(datestring)

#dummydatestr:
dd_anf=['2000-01-01','2000-01-09','2000-01-17','2000-01-25','2000-02-02','2000-02-10','2000-02-18','2000-02-26','2000-03-05','2000-03-13',
'2000-03-21','2000-03-29','2000-04-06','2000-04-14','2000-04-22','2000-04-30']  #16

ddd_anf=[]
for i in range(16):
	d=dt.datetime.strptime(dd_anf[i], "%Y-%m-%d")
	ddd_anf.append(d)
	
dd_end=['2011-09-22','2011-09-30','2011-10-08','2011-10-16','2011-10-24','2011-11-01','2011-11-09','2011-11-17','2011-11-25','2011-12-03',
'2011-12-11','2011-12-19','2011-12-27']  #13

ddd_end=[]
for i in range(13):
	d=dt.datetime.strptime(dd_end[i], "%Y-%m-%d")
	ddd_end.append(d)


dates=ddd_anf
for i in range(17*12):
	d=dt.datetime.strptime(datestr[i], "%Y-%m-%d")
	dates.append(d)



dates_all=dates+ddd_end
dates_all=array(dates_all)  #shape (642,)

nan16=[0.,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan]
nan13=[nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,0.]

###build relative mp fraction 12 yr array 
mp2000=mp2000.tolist()
mp2000=nan16+mp2000
mp2000.append(nan)
mp2001=mp2001.tolist()
mp2001.append(nan)
mp2002=mp2002.tolist()
mp2002.append(nan)
mp2003=mp2003.tolist()
mp2003.append(nan)
mp2004=mp2004.tolist()
mp2004.append(nan)
mp2005=mp2005.tolist()
mp2005.append(nan)
mp2006=mp2006.tolist()
mp2006.append(nan)
mp2007=mp2007.tolist()
mp2007.append(nan)
mp2008=mp2008.tolist()
mp2008.append(nan)
mp2009=mp2009.tolist()
mp2009.append(nan)
mp2010=mp2010.tolist()
mp2010.append(nan)
mp2011=mp2011.tolist()
mp2011.append(nan)
mp2011=mp2011+nan13
mp=mp2000+mp2001+mp2002+mp2003+mp2004+mp2005+mp2006+mp2007+mp2008+mp2009+mp2010+mp2011
mp=array(mp)

###build absulute mp area 12 yr array 
a_mp2000=a_mp2000.tolist()
a_mp2000=nan16+a_mp2000
a_mp2000.append(nan)
a_mp2001=a_mp2001.tolist()
a_mp2001.append(nan)
a_mp2002=a_mp2002.tolist()
a_mp2002.append(nan)
a_mp2003=a_mp2003.tolist()
a_mp2003.append(nan)
a_mp2004=a_mp2004.tolist()
a_mp2004.append(nan)
a_mp2005=a_mp2005.tolist()
a_mp2005.append(nan)
a_mp2006=a_mp2006.tolist()
a_mp2006.append(nan)
a_mp2007=a_mp2007.tolist()
a_mp2007.append(nan)
a_mp2008=a_mp2008.tolist()
a_mp2008.append(nan)
a_mp2009=a_mp2009.tolist()
a_mp2009.append(nan)
a_mp2010=a_mp2010.tolist()
a_mp2010.append(nan)
a_mp2011=a_mp2011.tolist()
a_mp2011.append(nan)
a_mp2011=a_mp2011+nan13
a_mp=a_mp2000+a_mp2001+a_mp2002+a_mp2003+a_mp2004+a_mp2005+a_mp2006+a_mp2007+a_mp2008+a_mp2009+a_mp2010+a_mp2011
a_mp=array(a_mp)







d_short=['2000-08-04','2001-08-04','2002-08-04','2003-07-03','2004-07-11','2005-06-25','2006-07-27','2007-06-09','2008-07-03','2009-08-04','2010-06-25','2011-06-25',]  #12

dd_short=[]
for i in range(12):
	d=dt.datetime.strptime(d_short[i], "%Y-%m-%d")
	dd_short.append(d)

dates_short=array(dd_short)
#fit



############Trendberechnung relativ mp fraction
t=arange(12)
#a=linregress(t,mpmax)

#Linear regressison -polyfit - polyfit can be used other orders polys
(ar,br)=polyfit(t,mpmax,1)
xr=polyval([ar,br],t)

#Linear regressison -polyfit - polyfit can be used other orders polys
(cr,dr)=polyfit(t,meanyear,1)
xxr=polyval([cr,dr],t)

###signifitanztest
m1=mean(meanyear)
std1=array([std(meanyear),std(meanyear),std(meanyear),std(meanyear),std(meanyear),std(meanyear),std(meanyear),std(meanyear),std(meanyear),std(meanyear),std(meanyear),std(meanyear)])

print 'relative mp fraction: '
print 'Mittelwert1: '+str(m1)
#Linear regression using stats.linregress
(a1,b1,r1,tt1,s1)=linregress(t,xxr)
print 'Steigung: '+str(a1)
print 'Stdabw1: '+str(s1)
print 'Stdabw2: '+str(std(meanyear))
print 'Korrelationscoeff1: '+str(r1)
print 'signifikanztest1: '+str(tt1)






######relative
#figure()
#ylabel('Melt pond fraction',fontsize=16)
#text(0.6,0.03,'mean: '+str(m))
#title('mean MODIS melt pond fraction - zeitreihe_4.py')
#ylim([0,0.4])
#plot(dates_all,mp,'k-',label='relative melt pond fraction')
#plot(dates_short,mpmax,'ro-',label='maximum relative melt pond fraction')
#plot(dates_short,xr,'r--',label='trend')
#plot(dates_short,meanyear,'bo-',label='mean relative melt pond fraction')
#plot(dates_short,xxr,'b--',label='trend')
#legend(loc=3)
#grid()
#savefig('/pf/u/u241127/plots/meltponds_timeline2000_2011_rel.png')


############Trendberechnung absolute mp area
t=arange(12)
#Linear regressison -polyfit - polyfit can be used other orders polys
(aar,bbr)=polyfit(t,a_maxyear,1)
xxxr=polyval([aar,bbr],t)
#Linear regressison -polyfit - polyfit can be used other orders polys
(ccr,ddr)=polyfit(t,a_meanyear,1)
xxxxr=polyval([ccr,ddr],t)

###signifitanztest
mm1=mean(a_meanyear)
std2=array([std(a_meanyear),std(a_meanyear),std(a_meanyear),std(a_meanyear),std(a_meanyear),std(a_meanyear),std(a_meanyear),std(a_meanyear),std(a_meanyear),std(a_meanyear),std(a_meanyear),std(a_meanyear)])
print ' '
print 'absolute mp area: '
print 'Mittelwert1: '+str(mm1)
#Linear regression using stats.linregress
(aa1,bb1,rr1,ttt1,ss1)=linregress(t,xxxxr)
print 'Steigung: '+str(aa1)
print 'Stdabw1: '+str(ss1)
print 'Stdabw2: '+str(std(a_meanyear))
print 'Korrelationscoeff1: '+str(rr1)
print 'signifikanztest1: '+str(ttt1)


#####absolute
#figure()
#ylabel('Melt pond area',fontsize=16)
#text(0.6,0.03,'mean: '+str(m))
#title('mean MODIS melt pond area - zeitreihe_4.py')
#ylim([0,2.6])
#plot(dates_all,a_mp,'k-',label='absolute melt pond area')
#plot(dates_short,a_maxyear,'ro-',label='maximum melt pond area')
#plot(dates_short,xxxr,'r--',label='trend')
#plot(dates_short,a_meanyear,'bo-',label='mean melt pond area')
#plot(dates_short,xxxxr,'b--',label='trend')
#legend(loc=3)
#grid()
#savefig('/pf/u/u241127/plots/meltponds_timeline2000_2011_abs.png')



########################twinx plot
#################double plot fot both days mean and max
#look at zeitreihe3.py

########################twinx plot
#################double plot fot both days mean

figure(figsize=(11,7))
ax1 = subplot(111)
ax1.set_ylabel('Melt pond fraction',color='DarkGreen',fontsize=16)
for tl in ax1.get_yticklabels():
    tl.set_color('DarkGreen')
    
    
ax1.set_ylim(0.00,0.3)

#text(0.6,0.03,'mean: '+str(m))
title('mean MODIS melt pond fraction and area - zeitreihe_5.py')
plot(dates_short,meanyear,'o-',color='Green',label='mean relative melt pond fraction')
plot(dates_short,xxr,'--',color='Green',label='trend')
legend(loc=3)
grid()
# 2nd y-axes 
ax2 = twinx()

ax2.set_ylabel('Melt pond area in Mio km2',color='DarkBlue',fontsize=16)
for tl in ax2.get_yticklabels():
    tl.set_color('DarkBlue')

#ylim([0,2.])
plot(dates_all,a_mp,'-',label='absolute melt pond area ',color='Gainsboro')
plot(dates_short,a_meanyear,'o-',label='mean melt pond area',color='Blue')
plot(dates_short,xxxxr,'--',label='trend',color='Blue')
legend(loc=4)
#grid()
savefig('/pf/u/u241127/plots/meltponds_timeline2000_2011_onlymeans_masked.png')

########################twinx plot in b&w
#################double plot fot both days mean


params = {'axes.labelsize': 13,
          'text.fontsize': 16,
          'legend.fontsize': 14,
          'xtick.labelsize': 15,
          'ytick.labelsize': 16,}
rcParams.update(params)

figure(figsize=(11,7))
ax1 = subplot(211)
ax1.set_ylabel('Melt pond fraction',fontsize=16)
ax1.set_ylim(0.21,0.28)
#text(0.6,0.03,'mean: '+str(m))
#title('mean MODIS melt pond fraction and area b&w - zeitreihe_5.py')
plot(dates_all,a_mp,'-',color='white')
#plot(dates_short,meanyear,'k^-',label='mean relative melt pond fraction',markersize=8)
errorbar(dates_short,meanyear,yerr=std1,fmt='k^-',label='mean relative melt pond fraction', markersize=8)
plot(dates_short,xxr,'k--',label='trend')
legend(numpoints=1,loc=3,frameon=False)
grid()

ax2=subplot(212, sharex=ax1)   #shared x-axis
ax2.set_ylim(1.2,1.9)
ax2.set_ylabel('Area [Mio km2]',fontsize=16)
plot(dates_all,a_mp,'-',color='white')
#plot(dates_short,a_meanyear,'ks-',label='mean melt pond area',markersize=8)
errorbar(dates_short,a_meanyear,yerr=std2,fmt='ks-',label='mean relative melt pond fraction', markersize=8)
plot(dates_short,xxxxr,'k--',label='trend')
legend(numpoints=1,loc=3,frameon=False)
grid()
savefig('/pf/u/u241127/plots/meltponds_timeline2000_2011_onlymeans_bw_masked.png')


#Trendcalculation:
deltaT1=xxr[-1]-xxr[1]
meanyear_first=meanyear[0]
T1=100-((meanyear_first-deltaT1)*100/meanyear_first)
print('Trend der rel melt pond fract zunahme ist '+str(T1)+' %')

#Trendcalculation:
deltaT=xxxxr[1]-xxxxr[-1]
a_meanyear_first=a_meanyear[0]
T=100-((a_meanyear_first-deltaT)*100/a_meanyear_first)
print('Trend der Flaechenabnahme ist '+str(T)+' %')
