# applay neural net from programm neural_net_train on full modis scene 
# based onclassifcation Tschudi 2008 for MODIS - haupt_mod_tschudi.py

#run neural_net_train_3_klasses_multiple_surfaces.py

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

import time

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
        #if x[2]>0.4:            #zur Korrektur der hohen Schmelztuempelwerte an der Wasserkante
            #print('x>0.4') 
            #=f+x[0]
    return f

def funktion(R):
    print('size of Matrix: '+str(R.shape[1])+' rows, '+str(R.shape[2])+' colums ')
    for l in range(R.shape[1]):
        print(str(l)+' of '+str(R.shape[1]))
        for n in range(R.shape[2]):
            Q=R[:,l,n]
            if isnan(sum(Q))==True:
                x=array([0.25,0.25,0.25])
                z=array([ NaN,  NaN,  NaN])
            else:
                x=array([0.25,0.25,0.25])
                z = opti.fmin_bfgs(my_cost, x, args = (Q,r),full_output=False,disp=False,retall=False,gtol=0.01)
                a_neu[:,l,n]=z
                f_orig[l,n]=dot((dot(r,x)-Q).transpose(),(dot(r,x)-Q))
                f_neu[l,n]=dot((dot(r,z)-Q).transpose(),(dot(r,z)-Q))
                residuals[l,n]=sum(z)
    return a_neu


####---snow1:
workfolder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/GRID_MOSAIC/'
tmpfolder='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/tmp/'
############
##METALIST##
############
#choose test_scene
metalist = glob.glob(workfolder+'2008*/')
print(metalist)
#nn_file='net5neu.data'
mosaic='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/GRID_MOSAIC/2008_169_mosaic/'
#only for no mp test
#mosaic='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/GRID_MOSAIC/2000_129_mosaic/'
print('now work on day '+mosaic)
jahr=mosaic[-16:-12]#2008
tag=mosaic[-11:-8]#'06'
workfile=mosaic+jahr+'_'+tag+'_mosaic'
file1,file2,file3,file4=workfile+'1.nc',workfile+'2.nc',workfile+'3.nc',workfile+'4.nc'
titlestr='MOD09_mosaic_'+str(tag)+'_'+str(jahr)
print(titlestr)
#read Netcdfs and extract Bands
mb1,mb2,mb3,mb4=NetCDFFile(file1),NetCDFFile(file2),NetCDFFile(file3),NetCDFFile(file4)
print('read nc') 

########################find subset coordinates possitions
"""#Franklin bay and surrounding
lon_y=-400. #center coordinates of the plot in polar stereo
lat_x=-2000.
f=500.
#canadian archipelaga fy ice plains
lon_y=-850. #center coordinates of the plot in polar stereo
lat_x=-1350.
f=100.
"""


  

#AEREAS={'canadian':[-850,-1350,100.],'fram':[-1200,700,200.],'multiyear':[-100,-1400,200.],'franklin':[-400,-2000,500.]} 
AEREAS={'multiyear':[-100,-1400,200.],'siberia_fy':[1300,500,200.],'canadian':[-850,-1350,100.],'fram':[-1200,700,200.]}
#'siberia_fy':[1300,500,200.]
#AEREAS={'franklin':[-400,-2000,500.]}
#only for no_mp test
#AEREAS={'no_mp':[-25,-2030,20.]}
#only ocean test
#AEREAS={'ocean':[-150,1750,50.]} 
for i,j in enumerate(AEREAS.values()):
    region=AEREAS.keys()[i]
    print region
    lon_y=j[0]
    lat_x=j[1]
    f=j[2]
    MX=array(mb1.variables['x'][:])
    MY=array(mb1.variables['y'][:])
    X0,X2,Y0,Y2=find_cut_area_coord(MY,MX,lon_y,lat_x,f)
    print('cut nc') 
    MX=array(mb1.variables['x'][X0:X2])
    MY=array(mb1.variables['y'][Y0:Y2])
    MB1=array(mb1.variables['z'][Y0:Y2,X0:X2])
    MB2=array(mb2.variables['z'][Y0:Y2,X0:X2])
    MB3=array(mb3.variables['z'][Y0:Y2,X0:X2])
    MB4=array(mb4.variables['z'][Y0:Y2,X0:X2])
    print('done')
    #plot of subset
    figure()
    subplot(212)
    imshow(MB1,origin='lower')
    subplot(211)
    display_rgb(MB1,MB2,MB3)
    title(titlestr+' '+region)
    show()
    MB1[MB1==0.]=NaN
    MB2[MB2==0.]=NaN
    MB3[MB3==0.]=NaN
    MB4[MB4==0.]=NaN
    #classification with slow classification####################
    eins=ones(MB1.shape)
    eins[isnan(MB1)==True]=NaN
    BAND1=MB1
    BAND2=MB2
    BAND3=MB3
    BAND4=MB4
    D=eins
    R=array([BAND1,BAND3,BAND4,D])
    #pickle.dump(R,open(tmpfolder+'R_'+titlestr+region+'.data','w')) 
    #3klassen
    r=array([[0.22,0.95,0.08],[0.16,0.95,0.08],[0.07,0.87,0.08],[1.,1.,1.]])
    #values from tschudi2008 for MOD09 
    a_neu=zeros((r.shape[1],R.shape[1],R.shape[2]))
    a_neu[0][isnan(MB1)==True]=NaN
    a_neu[1][isnan(MB1)==True]=NaN
    a_neu[2][isnan(MB1)==True]=NaN
    residuals=zeros((R.shape[1],R.shape[2]))
    residuals[isnan(MB1)==True]=NaN
    f_orig=zeros((R.shape[1],R.shape[2]))
    f_orig[isnan(MB1)==True]=NaN
    f_neu=zeros((R.shape[1],R.shape[2]))
    f_neu[isnan(MB1)==True]=NaN
    #zeitnahme start
    start = time.clock()
    a_neu=funktion(R)
    #zeitnahme ende
    stend = time.clock()
    print('dauer: '+str(stend-start))
    #####offest abziehen (ice test)
    #aus Test ergibt sich:
    #offset_mp=0.02424182
    #offset_ice=0.024644
    #offset_ow=0.0233015
    #mp_neu=a_neu[0]-offset_mp
    #ice_neu=a_neu[1]+offset_ice
    #water_neu=a_neu[2]-offset_ow
    #a_neu[0]=mp_neu
    #a_neu[1]=ice_neu
    #a_neu[2]=water_neu
    #####offest abziehen (water test)
    #aus Test ergibt sich
    #offset_mp=0.026807243631264285
    #offset_ice=0.01925828282936267
    #offset_ow=0.048327398747455286
    #mp_neu=a_neu[0]-offset_mp
    #ice_neu=a_neu[1]-offset_ice
    #water_neu=a_neu[2]+offset_ow
    #a_neu[0]=mp_neu
    #a_neu[1]=ice_neu
    #a_neu[2]=water_neu
    mp_neu=a_neu[0]
    ice_neu=a_neu[1]
    water_neu=a_neu[2]
    mp_real=mp_neu-water_neu
    mp_real[mp_real<0.]=0.
    mp_real=a_neu[0]
    #mp_neu[mp_neu>1.0]=1.
    #mp_neu[mp_neu<0.0]=0.00001
    #ice_neu[ice_neu>1.0]=1.
    #ice_neu[ice_neu<0.0]=0.00001
    #water_neu[water_neu>1.0]=1.
    #water_neu[water_neu<0.0]=0.00001
    mp_neu[mp_neu>1.0]=nan
    mp_neu[mp_neu<0.0]=nan
    ice_neu[ice_neu>1.0]=nan
    ice_neu[ice_neu<0.0]=nan
    water_neu[water_neu>1.0]=nan
    water_neu[water_neu<0.0]=nan
    pickle
    pickle.dump(a_neu,open(tmpfolder+'a_neu_fmin_'+titlestr+region+'.data','w')) 
    pickle.dump(R,open(tmpfolder+'R_'+titlestr+region+'.data','w')) 
    #####testplot
    figure(figsize=(15,10))
    subplot(131)
    imshow(mp_real,origin='lower')
    colorbar()
    title('mp real fraction '+region)
    subplot(132)
    imshow(ice_neu,origin='lower')
    colorbar()
    title('ice fraction '+region)
    subplot(133)
    imshow(water_neu,origin='lower')
    colorbar()
    title('ow fraction '+region)
    show()




#array anlegen um neurales Netz zu trainieren
a_fram=pickle.load(open(tmpfolder+'a_neu_fmin_'+titlestr+'fram.data','r'))
a_canadian=pickle.load(open(tmpfolder+'a_neu_fmin_'+titlestr+'canadian.data','r'))
a_multiyear=pickle.load(open(tmpfolder+'a_neu_fmin_'+titlestr+'multiyear.data','r'))
#a_franklin=pickle.load(open(tmpfolder+'a_neu_fmin_'+titlestr+'franklin.data','r'))
#a_ocean=pickle.load(open(tmpfolder+'a_neu_fmin_'+titlestr+'ocean.data','r'))
a_siberia_fy=pickle.load(open(tmpfolder+'a_neu_fmin_'+titlestr+'siberia_fy.data','r'))
#a_no_mp=pickle.load(open(tmpfolder+'a_neu_fmin_MOD09_mosaic_129_2000no_mp.data','r'))

R_fram=pickle.load(open(tmpfolder+'R_'+titlestr+'fram.data','r'))
R_canadian=pickle.load(open(tmpfolder+'R_'+titlestr+'canadian.data','r'))
R_multiyear=pickle.load(open(tmpfolder+'R_'+titlestr+'multiyear.data','r'))
#R_franklin=pickle.load(open(tmpfolder+'R_'+titlestr+'franklin.data','r'))
#R_ocean=pickle.load(open(tmpfolder+'R_'+titlestr+'ocean.data','r'))
R_siberia_fy=pickle.load(open(tmpfolder+'R_'+titlestr+'siberia_fy.data','r'))
#R_no_mp=pickle.load(open(tmpfolder+'R_MOD09_mosaic_129_2000no_mp.data','r'))

shape1=a_fram.shape[1]+a_canadian.shape[1]+a_multiyear.shape[1]+a_siberia_fy.shape[1]#+a_no_mp.shape[1]
shape2=a_fram.shape[2]+a_canadian.shape[2]+a_multiyear.shape[2]+a_siberia_fy.shape[2]#+a_no_mp.shape[2]

a_allTestSets=zeros((3,shape1,shape2))
a_allTestSets[:,:a_fram.shape[1],:a_fram.shape[2]]=a_fram
a_allTestSets[:,a_fram.shape[1]:a_fram.shape[1]+a_canadian.shape[1],a_fram.shape[2]:a_fram.shape[2]+a_canadian.shape[2]]=a_canadian
a_allTestSets[:,a_fram.shape[1]+a_canadian.shape[1]:a_fram.shape[1]+a_canadian.shape[1]+a_multiyear.shape[1],a_fram.shape[2]+a_canadian.shape[2]:a_fram.shape[2]+a_canadian.shape[2]+a_multiyear.shape[2]]=a_multiyear
#a_allTestSets[:,a_fram.shape[1]+a_canadian.shape[1]+a_multiyear.shape[1]:a_fram.shape[1]+a_canadian.shape[1]+a_multiyear.shape[1]+a_no_mp.shape[1],a_fram.shape[2]+a_canadian.shape[2]+a_multiyear.shape[2]:a_fram.shape[2]+a_canadian.shape[2]+a_multiyear.shape[2]+a_no_mp.shape[2]]=a_no_mp
a_allTestSets[:,a_fram.shape[1]+a_canadian.shape[1]+a_multiyear.shape[1]:,a_fram.shape[2]+a_canadian.shape[2]+a_multiyear.shape[2]:]=a_siberia_fy
a_allTestSets[a_allTestSets==0.]=NaN

R_allTestSets=zeros((4,shape1,shape2))
#R_allTestSets=zeros((3,shape1,shape2))
R_allTestSets[:,:a_fram.shape[1],:a_fram.shape[2]]=R_fram
R_allTestSets[:,a_fram.shape[1]:a_fram.shape[1]+a_canadian.shape[1],a_fram.shape[2]:a_fram.shape[2]+a_canadian.shape[2]]=R_canadian
R_allTestSets[:,a_fram.shape[1]+a_canadian.shape[1]:a_fram.shape[1]+a_canadian.shape[1]+a_multiyear.shape[1],a_fram.shape[2]+a_canadian.shape[2]:a_fram.shape[2]+a_canadian.shape[2]+a_multiyear.shape[2]]=R_multiyear

#R_allTestSets[:,a_fram.shape[1]+a_canadian.shape[1]+a_multiyear.shape[1]:a_fram.shape[1]+a_canadian.shape[1]+a_multiyear.shape[1]+a_no_mp.shape[1],a_fram.shape[2]+a_canadian.shape[2]+a_multiyear.shape[2]:a_fram.shape[2]+a_canadian.shape[2]+a_multiyear.shape[2]+a_no_mp.shape[2]]=R_no_mp

R_allTestSets[:,a_fram.shape[1]+a_canadian.shape[1]+a_multiyear.shape[1]:,a_fram.shape[2]+a_canadian.shape[2]+a_multiyear.shape[2]:]=R_siberia_fy
R_allTestSets[R_allTestSets==0.]=NaN
R=R_allTestSets.copy()
"""
###18 may
a_allTestSets=zeros((3,shape1,shape2))
a_allTestSets[:,:a_fram.shape[1],:a_fram.shape[2]]=a_fram
a_allTestSets[:,a_fram.shape[1]:a_fram.shape[1]+a_canadian.shape[1],a_fram.shape[2]:a_fram.shape[2]+a_canadian.shape[2]]=a_canadian
a_allTestSets[:,a_fram.shape[1]+a_canadian.shape[1]:a_fram.shape[1]+a_canadian.shape[1]+a_multiyear.shape[1],a_fram.shape[2]+a_canadian.shape[2]:a_fram.shape[2]+a_canadian.shape[2]+a_multiyear.shape[2]]=a_multiyear
a_allTestSets[:,a_fram.shape[1]+a_canadian.shape[1]+a_multiyear.shape[1]:,a_fram.shape[2]+a_canadian.shape[2]+a_multiyear.shape[2]:]=a_franklin
a_allTestSets[a_allTestSets==0.]=NaN

R_allTestSets=zeros((4,shape1,shape2))
#R_allTestSets=zeros((3,shape1,shape2))
R_allTestSets[:,:a_fram.shape[1],:a_fram.shape[2]]=R_fram
R_allTestSets[:,a_fram.shape[1]:a_fram.shape[1]+a_canadian.shape[1],a_fram.shape[2]:a_fram.shape[2]+a_canadian.shape[2]]=R_canadian
R_allTestSets[:,a_fram.shape[1]+a_canadian.shape[1]:a_fram.shape[1]+a_canadian.shape[1]+a_multiyear.shape[1],a_fram.shape[2]+a_canadian.shape[2]:a_fram.shape[2]+a_canadian.shape[2]+a_multiyear.shape[2]]=R_multiyear
R_allTestSets[:,a_fram.shape[1]+a_canadian.shape[1]+a_multiyear.shape[1]:,a_fram.shape[2]+a_canadian.shape[2]+a_multiyear.shape[2]:]=R_franklin
R_allTestSets[R_allTestSets==0.]=NaN
R=R_allTestSets.copy()

"""


#creating values for neural network
Z_opti_NN=a_allTestSets.copy().flatten()
Z_opti_NN[isnan(Z_opti_NN)==True]=0
ind_Z_opti_NN=Z_opti_NN.nonzero()
Z_opti_NN=Z_opti_NN[Z_opti_NN.nonzero()]
Z_opti_NN=Z_opti_NN.reshape(3,Z_opti_NN.shape[0]/3)
print('a_allTestSets.shape: '+str(a_allTestSets.shape))
print('Z_opti_NN.shape: '+str(Z_opti_NN.shape))
R_NN =R[0:3,:].copy().flatten()
R_NN[isnan(R_NN)==True]=0
R_NN=R_NN[ind_Z_opti_NN]
R_NN=R_NN.reshape(3,R_NN.shape[0]/3)
print('R[0:3,:].shape: '+str(R[0:3,:].shape))
print('R_NN.shape: '+str(R_NN.shape))


########## Define network geometry#############
#conec = mlgraph((3,5,3)) # Define network geometry 3 input layers, 5 middle, 3 output neurons
#conec = mlgraph((3,9,15,9,3))
#conec = mlgraph((3,5,7,5,3))
#conec = mlgraph((3,5,11,5,3))
#conec = mlgraph((3,21,21,3)) #--->ging ganz gut
#conec = mlgraph((3,12,12,3))
#conec = mlgraph((3,15,15,3))
#conec = mlgraph((3,27,27,3))#--->result pickeld as net4.data -> ging gut
conec = mlgraph((3,27,27,3))#--->result pickeld as net5.data -> 
#conec = mlgraph((3,20,3))
net = ffnet(conec)

LS=1000 # Learning steps
#LS=5000 # Learning steps

start = time.clock()
print 'Training'
net.randomweights()
net.train_tnc(R_NN[:3,:].transpose(), Z_opti_NN.transpose(), maxfun = LS)

stend = time.clock()
print('dauer training: '+str(stend-start))
#pickle
pickle.dump(net,open(tmpfolder+'neural_net_my_sib_can_fram','w')) 
#pickle.dump(net,open(tmpfolder+'neural_net_my_sib_can_fram_MP_valuesTSCHUDI','w')) 


print  'Testing'
output, regression = net.test(R_NN[:3,:].transpose(),Z_opti_NN.transpose())

print 'Applying'

R1=R[:3,:].flatten()
R=R1.reshape(3,R1.shape[0]/3)
Z_neural=net(R[:3,:].transpose()).transpose()
Z_neural_plot=net(R_NN[:3,:].transpose()).transpose()


print 'Plotting'
for j in range(3):
    figure()
    plot(Z_opti_NN[j,:],Z_neural_plot[j,:],'.')


show()



mp_neu=Z_neural[0].reshape((R_allTestSets.shape[1],R_allTestSets.shape[2]))
ice_neu=Z_neural[1].reshape((R_allTestSets.shape[1],R_allTestSets.shape[2]))
water_neu=Z_neural[2].reshape((R_allTestSets.shape[1],R_allTestSets.shape[2]))

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

"""


    #classification with neural network####################
    BAND1=MB1.flatten()
    BAND2=MB2.flatten()
    BAND3=MB3.flatten()
    BAND4=MB4.flatten()
    Rf=array([BAND1,BAND3,BAND4])

    net=pickle.load(open(tmpfolder+nn_file,'r'))



    #zeitnahme start
    start = time.clock()

    print 'Applying'
    Z_neural=net(Rf.transpose()).transpose()


    #zeitnahme ende
    stend = time.clock()
    print('dauer: '+str(stend-start))

    mp=Z_neural[0].reshape((MB1.shape))
    snow=Z_neural[1].reshape((MB1.shape))
    ow=Z_neural[2].reshape((MB1.shape))

    mp[mp>1.0]=NaN
    mp[mp<0.0]=NaN
    snow[snow>1.0]=NaN
    snow[snow<0.0]=NaN
    ow[ow>1.0]=NaN
    ow[ow<0.0]=NaN

    MPFraction=nansum(mp.copy())
    Waterfraction=nansum(ow.copy())
    mp_neu1=mp.copy()
    water_neu1=ow.copy()
    mp_neu1[isnan(mp_neu1)==False]=1.
    water_neu1[isnan(water_neu1)==False]=1.
    Pixels=nansum(mp_neu1)
    MP=MPFraction/Pixels
    water=Waterfraction/Pixels
    totalarea1=0.25**2*Pixels
    waterarea=totalarea1*water
    meltpondarea=(totalarea1-waterarea)*MP

    ratio1=meltpondarea/totalarea1
    ratio2=waterarea/totalarea1
    ratio_realMP=ratio1-ratio2
    print('meltpondarea = '+str(meltpondarea)+' sqkm of '+str(totalarea1-waterarea)+' sqkm - MP Fraction = '+str(ratio1))
    
    #plot of results neureal net method
    #figure()
    #subplot(311)
    #display_mono(mp,cmap='jet')
    #title(titlestr+' meltpond fraction NN: '+str(ratio1))
    #subplot(312)
    #display_mono(snow,cmap='jet')
    #title('ice fraction NN')
    #subplot(313)
    #display_mono(ow,cmap='jet')
    #title('water fraction NN')
    #colorbar(orientation='horizontal',shrink=0.5)
    #savefig('/pf/u/u241127/plots/'+titlestr+'_nn.png')

    print('PLOT')
    #plot all
    figure(figsize=(10,15))
    suptitle(titlestr+' - '+nn_file)
    subplot(321)
    display_mono(mp_neu,cmap='jet')
    title(' meltpond fraction fmin: '+str(ratio1_fmin))
    subplot(323)
    display_mono(ice_neu,cmap='jet')
    title('ice fraction fmin')
    subplot(325)
    display_mono(water_neu,cmap='jet')
    title('water fraction fmin')
    subplot(322)
    display_mono(mp,cmap='jet')
    title(' meltpond fraction NN: '+str(ratio1))
    subplot(324)
    display_mono(snow,cmap='jet')
    title('ice fraction NN')
    subplot(326)
    display_mono(ow,cmap='jet')
    title('water fraction NN')
    cax = axes([0.2, 0.03, 0.6, 0.01])
    colorbar(cax=cax, orientation='horizontal')
    #savefig('/pf/u/u241127/plots/'+titlestr+'_'+nn_file+'_results.png')
    #show()
    del MX,MY,MB1,MB2,MB3,MB4

    
"""

