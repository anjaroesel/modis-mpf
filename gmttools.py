from scipy import *
import scipy.ndimage 
from pyhdf.SD import *
import sys,pipes,struct,os,glob,fnmatch
from polar_projection import *
from caselist import *
import cPickle
from mpl_toolkits.basemap import NetCDFFile

home=os.getenv('HOME')

def cmd(cmdline):
    print cmdline
    os.system(cmdline)
    return

def tv(a):
    figure()
    imshow(a,interpolation='nearest',origin='lower')       
    colorbar()
    show()

def grid_par_reg(region,cs,G):
    regions=def_regions()
    print regions
    sgn,cds= regions[region][0],regions[region][1:5]
    x0,y0,x2,y2=cds[0],cds[1],cds[2],cds[3]
    xs,ys=(x0-x2)/cs,(y0-y2)/cs
    G['Rx']=' -R'+str(x2)+'/'+str(x0)+'/'+str(y2)+'/'+str(y0)+' '
    xwidth,ywidth=xs*cs,ys*cs
    lat2,lon2=mapxy(x0,y0,sgn)
    lat0,lon0=mapxy(x2,y2,sgn)
    map_x_size, map_y_size = 1.65, 2.45     #A4
    map_x_size, map_y_size = 1.5, 2.22     #A4 smaller
    map_ratio=map_x_size/map_y_size
    img_ratio=xwidth/ywidth
    width=xwidth
    paper=map_x_size
    #paper=map_y_size #Landscape
    latlon_corners=str(lon0)+'/'+str(lat0)+'/'+str(lon2)+'/'+str(lat2)
    scalell=long(width*10000.0/paper)
    if sgn==1:
        grid = '-Js-45/90/70/1:'
    if sgn==-1:
        grid = '-Js0/-90/-70/1:'
    latlon_par=' -R'+latlon_corners+'r '+grid+str(scalell)+' '
    latlon=' -R'+latlon_corners+'r '
    G['Rlatlon']=latlon
    G['Rll']=latlon_par
    scale=1./width*10.0*paper
    G['Jx']=' -Jx'+str(scale)+'c '
    G['I']=' -I'+str(cs)+' '
    G['S']=' -S'+str(cs*4)+' '  #Suchradius
    G['N']=' -N4 '
    #G['Bll']=' -Ba10g10/a2g2NesW ' # medium (use for LS single scene)
    #G['Bll']=' -Ba20g20/a10g10nesW ' # coarse (use for overviewplot)
    G['Bll']=' -Ba45g45/a10g10nesW ' # very coarse 
    #G['Bll']=' -Ba10mg10m/a5mg5mWESN ' # fine
    #G['Bx']=' -Ba100g50/a100g50nESw ' # (use for LS single scene)
    G['Bx']=' -Ba1000g500/a1000g500NEsw ' #coarse (use for overviewplot)
    G['Bx']=' -Ba1000/a2000NEsw ' #coarse no gridlines
    #G['coast']=' -Df -W1/0/255/0 '
    G['coast']=' -Dl -W0.5/0/0/0 ' #black 0/0/128  #-->blue
    return G

def grid_par_neu(sgn,x0,y0,x2,y2,cs,G): #(sgn,xmax,ymax,xmin,ymin,cs,G)
    xs,ys=(x0-x2)/cs,(y0-y2)/cs
    G['Rx']=' -R'+str(x2)+'/'+str(x0)+'/'+str(y2)+'/'+str(y0)+' '
    xwidth,ywidth=xs*cs,ys*cs
    lat2,lon2=mapxy(x0,y0,sgn)    
    lat0,lon0=mapxy(x2,y2,sgn)    
    map_x_size, map_y_size = 1.65, 2.45     #A4
    map_x_size, map_y_size = 1.5, 2.22     #A4 smaller
    map_ratio=map_x_size/map_y_size
    img_ratio=xwidth/ywidth
    width=xwidth
    paper=map_x_size
    #paper=map_y_size #Landscape
    latlon_corners=str(round(lon2,1))+'/'+str(round(lat0,1))+'/'+str(round(lon0,1))+'/'+str(round(lat2,1))             # geaendert
    scalell=long(width*10000.0/paper)
    if sgn==1:
        grid = '-Js-45/90/70/1:'
    if sgn==-1:
        grid = '-Js0/-90/-70/1:'
    latlon_par=' -R'+latlon_corners+'r '+grid+str(scalell)+' '
    latlon_pro=' -R'+str(round(lon2,1))+'/'+str(round(lon0,1))+'/'+str(round(lat0,1))+'/'+str(round(lat2,1))+' '+grid+str(scalell)+' '
    latlon=' -R'+latlon_corners+'r '
    G['Rlatlon']=latlon
    G['Rlatlonsam']=' -R'+str(round(lon2,1))+'/'+str(round(lon0,1))+'/'+str(round(lat0,1))+'/'+str(round(lat2,1))
    G['Rll']=latlon_par
    G['Rllpro']=latlon_pro
    scale=1/width*10.0*paper
    G['Jx']=' -Jx'+str(scale)+'c '
    G['I']=' -I'+str(cs)+'= '
    G['S']=' -S'+str(cs*4)+' '  #Suchradius
    G['N']=' -N4 '
    #G['Bll']=' -Ba10g10/a2g2NesW ' # medium
    #G['Bll']=' -Ba20g20/a10g10NesW ' # coarse
    G['Bll']=' -Ba40g40/a20g20NesW ' # very coarse
    G['Bll']=' -Ba10mg10m/a5mg5mWESN ' # fine
    G['Bx']=' -Ba100g50/a100g50nESw ' #medium
    #G['Bx']=' -Ba500g100/a500g100nESw ' #coarse
    #G['Bx']=' -Ba1000/a1000nESw ' #coarse no gridlines
    #G['coast']=' -Df -W2navy '#blue coastline
    G['coast']=' -Df -W2black '
    return G

#x0,y0,x2,y2=X_polar.max(),Y_polar.max(),X_polar.min(),Y_polar.min()
def grid_par_round(sgn,latmax,lonmax,latmin,lonmin,cs,G): #(sgn,latmax,lonmax,latmin,lonmin,cs,G)
    xmax,ymax=mapll(latmax,lonmax,sgn)    
    xmin,ymin=mapll(latmin,lonmin,sgn) 
    xs,ys=(xmax-xmin)/cs,(ymax-ymax)/cs
    G['Rx']=' -R'+str(xmin)+'/'+str(xmax)+'/'+str(ymin)+'/'+str(ymax)+' '
    xwidth,ywidth=xs*cs,ys*cs
    map_x_size, map_y_size = 1.65, 2.45     #A4
    map_x_size, map_y_size = 1.5, 2.22     #A4 smaller
    map_ratio=map_x_size/map_y_size
    img_ratio=xwidth/ywidth
    width=xwidth
    paper=map_x_size
    #paper=map_y_size #Landscape
    latlon_corners=str(round(lonmin,1))+'/'+str(round(lonmax,1))+'/'+str(round(latmin,1))+'/'+str(round(latmax,1))             # geaendert
    scalell=long(width*10000.0/paper)
    if sgn==1:
        grid = '-Js-45/90/70/1:'
    if sgn==-1:
        grid = '-Js0/-90/-70/1:'
    latlon_par=' -R'+latlon_corners+grid+str(scalell)+' '
    latlon=' -R'+latlon_corners+' '
    G['Rlatlon']=latlon
    G['Rll']=latlon_par
    scale=1/width*10.0*paper
    G['J']=' -Js-45/90/12c/70: '
    G['I']=' -I'+str(cs)+'= '
    G['S']=' -S'+str(cs*4)+' '  #Suchradius
    G['N']=' -N4 '
    G['B']=' -B30g30/15g15 ' # medium
    G['coast']=' -Df -W2black '
    return G

def makecpt(ctabfile,ctab,min,max,step,G,**kw):
    s=''
    if kw.has_key('I'):
        s=' -I '
    #cmd('makecpt -C'+ctab+s+' -T0/255/50 -D -Z >'+ctabfile)
    cmd('makecpt -C'+ctab+s+' -D -Z -T'+str(min)+'/'+str(max)+'/'+str(step)+' >'+ctabfile)
    G['C']=' -C'+ctabfile
    return G

def opts(G,opts):
    s=' '
    for o in opts:
        s=s+G[o]+' '
    return s

def nearneighbor(xab,yab,I,filename,G):
    cmd='nearneighbor  '+G['Rx']+' -V -F '+G['N']+G['I']+G['S']+' -bis -G'+filename
    pipe=pipes.Template()
    pipe.append(cmd,'-.')
    pipe.open(filename+'.pipe', 'w').write(array([xab.flatten(),yab.flatten(),I.flatten()]).swapaxes(0,1).tostring())
    return




def xyz2grd(xab,yab,I,filename,G):
    cmd='xyz2grd  '+G['Rx']+' -V -F  '+G['I']+' -bis -G'+filename
    pipe=pipes.Template()
    pipe.append(cmd,'-.')
    pipe.open(filename+'.pipe', 'w').write(array([xab.flatten(),yab.flatten(),I.flatten()]).swapaxes(0,1).tostring())
    return

def surface(xab,yab,I,filename,G):
    cmd='surface  '+G['Rx']+' -V -F '+G['N']+G['I']+G['S']+' -bis -G'+filename
    pipe=pipes.Template()
    pipe.append(cmd,'-.')
    pipe.open(filename+'.pipe', 'w').write(array([xab.flatten(),yab.flatten(),I.flatten()]).swapaxes(0,1).tostring())
    return

def readnc(filename):
    g=Nio.open_file(filename,'r')
    return g.variables['z'][:,:]

def write_table(filename,x,y,z):
    l=open(filename,'w')
    for i in range(z.shape[0]):
        l.write(str(x[i])+' '+str(y[i])+' '+str(z[i])+'\n') # e.g. Lat lon TB
    l.close()
    return

def x0y0x2y2(region):
    regions=def_regions()
    sgn,cds= regions[region][0],regions[region][1:5]
    return cds[0],cds[1],cds[2],cds[3]

def xsys(region,cs):
    x0,y0,x2,y2=x0y0x2y2(region)
    xs,ys=int((x0-x2)/cs),int((y0-y2)/cs)
    return xs,ys

def XYgrid(region,cs):
    x0,y0,x2,y2=x0y0x2y2(region)
    y=linspace(y2+cs/2,y0-cs/2,(y0-y2)/cs)
    x=linspace(x2+cs/2,x0-cs/2,(x0-x2)/cs)
    return meshgrid(x,y)

def XYgrid2(region,cs):
    x0,y0,x2,y2=x0y0x2y2(region)
    y=linspace(y2+cs/2,y0-cs/2,(y0-y2)/cs)
    x=linspace(x2+cs/2,x0-cs/2,(x0-x2)/cs)
    return x,y 

def get_landmask(geogr_file,polarstereogr_file): 
    path='/scratch/clisap/seaice/OWN_PRODUCTS/MELT_PONDS/tmp/'
    if not os.path.exists(path+polarstereogr_file):
        os.system('grdlandmask '+opts(G,['Rlatlon','I'])+' -Df -N1/NaN/1/NaN/NaN -G'+path+geogr_file+' -V') #-Nocean/land/lake/island/pond.
        #os.system('grdcut '+path+geogr_file+xy_region+'-G'+path+polarstereogr_file+' -V ')
        os.system('grdproject '+path+geogr_file+' '+opts(G,['Rll'])+' -A -C -G'+path+polarstereogr_file+' -V')
    data=NetCDFFile(path+polarstereogr_file)
   #print data,data.dimensions,data.variables.keys():
    x=array(data.variables['x'][:])
    y=array(data.variables['y'][:])   
    z=array(data.variables['z'][:,:])
    #z[z<=0.5]=0
    #z[z>0.5]=1
    #zz=array(z,dtype=int)
    return x,y,z



