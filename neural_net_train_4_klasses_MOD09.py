#training programm for neural net based on
#classifcation Tschudi 2008 for MODIS - haupt_mod_tschudi.py


import sys,pipes,struct,os,glob,fnmatch,pickle
from mpl_toolkits.basemap import NetCDFFile
from scipy.interpolate import *
from sattools import *
from scipy.ndimage import *



from ffnet import ffnet,mlgraph
import networkx

from pylab import *
from scipy.interpolate import *
import scipy.io as io
import scipy.optimize as opti

#####select date:
date='2008_153'

####---snow1:
workfolder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/GRIDDATAMOD09/MOD09A1/'+date+'_mosaic/'
tmpfolder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/tmp/'
############
##METALIST##
############/
metalist = glob.glob(workfolder+'*.nc')
print(metalist)

jahr=date[0:4]
julday=date[5:8]
workfile=workfolder+date+'_mosaic'
file1,file2,file3,file4=workfile+'1.nc',workfile+'2.nc',workfile+'3.nc',workfile+'4.nc'
#read Netcdfs and extract Bands
mb1,mb2,mb3,mb4=NetCDFFile(file1),NetCDFFile(file2),NetCDFFile(file3),NetCDFFile(file4)
print('read nc') 



#MX=array(mb1.variables['x'][:])
#MY=array(mb1.variables['y'][:])
#MB1=array(mb1.variables['z'][:])
#MB2=array(mb2.variables['z'][:,:])
#MB3=array(mb3.variables['z'][:,:])
#MB4=array(mb4.variables['z'][:,:])
#print('done')



#nur fuer TEST 2008_153 area:


Y0=9200
Y2=10100
X0=3200
X2=3800



#Y0=9200
#Y2=9400
#X0=3320
#X2=3400



#clip:
#nur fuer TEST 2008_6_06_05*/ area (small):
#Y0=2850
#Y2=2900
#X0=2700
#X2=2800
#nur fuer TEST 2008_6_06_05*/ area (large):
#Y0=2600
#Y2=3100
#X0=2600
#X2=3400

#nur fuer TEST 2008_6_07_0430*/ area (large):
#Y0=3000
#Y2=4000
#X0=2000
#X2=3000


MB1=array(mb1.variables['z'][Y0:Y2,X0:X2])
MB2=array(mb2.variables['z'][Y0:Y2,X0:X2])
MB3=array(mb3.variables['z'][Y0:Y2,X0:X2])
MB4=array(mb4.variables['z'][Y0:Y2,X0:X2])

figure(1)
imshow(MB1)
title(titlestr)
show()

#rgb darstellung
figure(2)
display_rgb(MB1,MB2,MB3)
title(titlestr)


MB1[MB1==0.]=NaN
MB2[MB2==0.]=NaN
MB3[MB3==0.]=NaN
MB4[MB4==0.]=NaN
eins=ones(MB1.shape)
eins[isnan(MB1)==True]=NaN


####tests

def extractXYZ(b1,b2,b3,y1,y2,x1,x2):
    X=b1[y1:y2,x1:x2].flatten()
    Y=b2[y1:y2,x1:x2].flatten()
    Z=b3[y1:y2,x1:x2].flatten()
    return X,Y,Z

"""
OW1,OW3,OW4=extractXYZ(MB1,MB3,MB4,570,620,300,360)
OW1=OW1.mean()
OW3=OW3.mean()
OW4=OW4.mean()
I1,I3,I4=extractXYZ(MB1,MB3,MB4,130,140,454,460)
I1=I1.mean()
I3=I3.mean()
I4=I4.mean()
MP1,MP3,MP4=extractXYZ(MB1,MB3,MB4,794,806,99,105)
MP1=MP1.mean()
MP3=MP3.mean()
MP4=MP4.mean()
print('OW1='+str(OW1)+',OW3='+str(OW3)+',OW4='+str(OW4))
print('I1='+str(I1)+',I3='+str(I3)+',I4='+str(I4))
print('MP1='+str(MP1)+',MP3='+str(MP3)+',MP4='+str(MP4))
"""

BAND1=MB1.flatten()
BAND2=MB2.flatten()
BAND3=MB3.flatten()
BAND4=MB4.flatten()
D=eins.flatten()

R=array([BAND1,BAND3,BAND4,D])
#r-werte (spectral reflectance werte) aus dem markus-bild abgelesen und gemittelt
#r=array([[pond1_b1,pond2_b1,ice_b1,water_b1],[pond1_b3,pond2_b3,ice_b3,water_b3],[pond1_b4,pond2_b4,ice_b4,water_b4],[1.,1.,1.,1.]])
#4klassen
#r=array([[0.74,0.46,0.8,0.18],[0.47,0.13,0.72,0.06],[0.31,0.06,0.67,0.03],[1.,1.,1.,1.]])

#3klassen - lief ganz gut
#r=array([[0.22,0.95,0.003],[0.16,0.95,0.02],[0.07,0.87,0.02],[1.,1.,1.]])

#values from tschudi2008 for MOD09 -----> komische Ergebnisse!!!!!!!
#4 klassen: pond, ice, snow, water
r=array([[0.22,0.76,0.95,0.003],[0.16,0.75,0.95,0.02],[0.07,0.56,0.87,0.02],[1.,1.,1.,1.]])

#fehleroptimierung

def null_bis_eins_smooth(x,alpha):
    return tanh(x*alpha)-tanh((x-1)*alpha)


alpha=150.
gewicht=0.1


def my_cost(x,Q,r):
    f=dot((dot(r,x)-Q).transpose(),(dot(r,x)-Q))
    neben=0.0  
    for i in range(3):
        neben=neben+1-null_bis_eins_smooth(x[i],alpha) #nebenbedingung (side condition)
    f=f+neben*gewicht
    return f


Z_opti=zeros((r.shape[1],R.shape[1]))
Z_opti[0][isnan(BAND1)==True]=NaN
Z_opti[1][isnan(BAND1)==True]=NaN
Z_opti[2][isnan(BAND1)==True]=NaN
Z_opti[3][isnan(BAND1)==True]=NaN
residuals=zeros(R.shape[1])
residuals[isnan(BAND1)==True]=NaN
f_orig=zeros(R.shape[1])
f_orig[isnan(BAND1)==True]=NaN
f_neu=zeros(R.shape[1])
f_neu[isnan(BAND1)==True]=NaN

########## Define network geometry#############
#conec = mlgraph((3,5,3)) # Define network geometry 3 input layers, 5 middle, 3 output neurons
conec = mlgraph((3,9,27,3))#--->result pickeld as net5.data -> 
net = ffnet(conec)

LS=5000 # Learning steps

# Optimization
x=array([0.25,0.25,0.25,0.25])


#zeitnahme start
import time
start = time.clock()
print('size of Matrix: '+str(R.shape[1])+' values ')

print 'Optimization'

def f(R):
    for i in range(R.shape[1]):
        print(str(i)+' of '+str(R.shape[1]))
        Q=R[:,i]
        if isnan(sum(Q))==True:
            x=array([0.25,0.25,0.25,0.25])
            z=array([ NaN,  NaN, NaN, NaN])
        else:
            x=array([0.25,0.25,0.25,0.25])
            z = opti.fmin_bfgs(my_cost, x, args = (Q, r),full_output=False,disp=False,retall=False,gtol=0.01)
        Z_opti[:,i]=z
        f_orig[i]=dot((dot(r,x)-Q).transpose(),(dot(r,x)-Q))
        f_neu[i]=dot((dot(r,z)-Q).transpose(),(dot(r,z)-Q))
        residuals[i]=sum(z)
    return Z_opti

Z_opti=f(R)

#zeitnahme ende
stend = time.clock()
print('dauer: '+str(stend-start))

#pickle Z_opti
#pickle.dump(Z_opti,open(tmpfolder+'Z_opti_neu.data','w')) # writes object to file

#plots
Z_opti_mp=Z_opti[0].reshape(MB1.shape)
Z_opti_ice=Z_opti[1].reshape(MB1.shape)
Z_opti_sn=Z_opti[2].reshape(MB1.shape)
Z_opti_ow=Z_opti[3].reshape(MB1.shape)

figure(3)
colorbar(orientation='horizontal')
subplot(221)
imshow(Z_opti_mp,origin='lower')
colorbar()
title('mp fraction')
subplot(222)
imshow(Z_opti_ice,origin='lower')
colorbar()
title('ice fraction')
subplot(223)
imshow(Z_opti_sn,origin='lower')
colorbar()
title('snow fraction')
subplot(224)
imshow(Z_opti_ow,origin='lower')
colorbar()
title('ow fraction')
show()




#creating values for neural network
Z_opti_NN=Z_opti.copy().flatten()
Z_opti_NN[isnan(Z_opti_NN)==True]=0
ind_Z_opti_NN=Z_opti_NN.nonzero()
Z_opti_NN=Z_opti_NN[Z_opti_NN.nonzero()]
Z_opti_NN=Z_opti_NN.reshape(3,Z_opti_NN.shape[0]/3)
print('Z_opti.shape: '+str(Z_opti.shape))
print('Z_opti_NN.shape: '+str(Z_opti_NN.shape))
R_NN =R[0:3,:].copy().flatten()
R_NN[isnan(R_NN)==True]=0
R_NN=R_NN[ind_Z_opti_NN]
R_NN=R_NN.reshape(3,R_NN.shape[0]/3)
print('R[0:3,:].shape: '+str(R[0:3,:].shape))
print('R_NN.shape: '+str(R_NN.shape))




start = time.clock()
print 'Training'
net.randomweights()
net.train_tnc(R_NN[:3,:].transpose(), Z_opti_NN.transpose(), maxfun = LS)

stend = time.clock()
print('dauer training: '+str(stend-start))
#pickle
pickle.dump(net,open(tmpfolder+'net5neu.data','w')) 


print  'Testing'
output, regression = net.test(R_NN[:3,:].transpose(),Z_opti_NN.transpose())

print 'Applying'
Z_neural=net(R[:3,:].transpose()).transpose()

Z_neural_plot=net(R_NN[:3,:].transpose()).transpose()

print 'Plotting'
for j in range(3):
    figure()
    plot(Z_opti_NN[j,:],Z_neural_plot[j,:],'.')


show()













mp_neu=Z_neural[0].reshape((MB1.shape))
ice_neu=Z_neural[1].reshape((MB1.shape))
water_neu=Z_neural[2].reshape((MB1.shape))

figure()
display_mono(mp_neu,cmap='jet')
title('meltpond1fraction neu')
colorbar(orientation='horizontal')
figure()
display_mono(ice_neu,cmap='jet')
title('icefraction neu')
colorbar()
figure()
display_mono(water_neu,cmap='jet')
title('waterfraction neu')
colorbar()


MPFraction=nansum(mp_neu.copy())
Waterfraction=nansum(water_neu.copy())
mp_neu1=mp_neu.copy()
water_neu1=water_neu.copy()
mp_neu1[isnan(mp_neu1)==False]=1.
water_neu1[isnan(water_neu1)==False]=1.
Pixels=nansum(mp_neu1)
mp=MPFraction/Pixels
water=Waterfraction/Pixels
totalarea1=0.25**2*Pixels
waterarea=totalarea1*water
meltpondarea=(totalarea1-waterarea)*mp

ratio1=meltpondarea/totalarea1
ratio2=waterarea/totalarea1
ratio_realMP=ratio1-ratio2
print('meltpondarea = '+str(meltpondarea)+' sqkm of '+str(totalarea1-waterarea)+' sqkm - MP Fraction = '+str(ratio1))

