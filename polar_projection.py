#module general/polar_projection.py
# LK 14/08/2008,4/2/2009

from numpy import *



# Polar stereographic (NSDIDC grid)
def mapxy(x, y, sgn):
    cdr  = 57.29577951
    slat = 70.0
    re   = 6378.273
    ec2   = .006693883
    ec    =  sqrt(ec2)
    if sgn == 1:
        delta = 45.0
    else:
        delta = 0.0
    sl = slat * pi/180.0
    rho = sqrt(x**2 + y**2)
    cm = cos(sl) / sqrt(1.0 - ec2 * (sin(sl)**2))
    t = tan((pi / 4.0) - (sl / 2.0)) / ((1.0 - ec * sin(sl)) / (1.0 + ec * sin(sl)))**(ec / 2.0)
    if  (absolute(slat-90.0) < 1.e-5):
        t = rho * sqrt((1. + ec)**(1. + ec) * (1. - ec)**(1. - ec)) / 2. / re
    else:
        t = rho * t / (re * cm)
    chi = (pi / 2.0) - 2.0 * arctan(t)
    alat = chi + ((ec2 / 2.0) + (5.0 * ec2**2.0 / 24.0) + (ec2**3.0 / 12.0)) * sin(2 * chi) + ((7.0 * ec2**2.0 / 48.0) + (29.0 * ec2**3 / 240.0)) *sin(4.0 * chi)+ (7.0 * ec2**3.0 / 120.0) * sin(6.0 * chi)
    alat = sgn * alat
    along = arctan2(sgn * x, -sgn * y)
    along = sgn * along

    along = along * 180. / pi
    alat  = alat * 180. / pi
    along = along - delta
    return [alat,along]

def mapll(lat,lon,sgn):
    cdr  = 57.29577951
    slat = 70.
    re   = 6378.273
    ec2   = .006693883
    ec    =  sqrt(ec2)
    if sgn == 1:
        delta = 45.0
    else:
        delta = 0.0
    latitude  = absolute(lat) * pi/180.
    longitude = (lon + delta) * pi/180.
    t =   tan(pi/4-latitude/2)/((1-ec*sin(latitude))/(1+ec*sin(latitude)))**(ec/2)
    if ((90-slat) < 1.e-5):
        rho = 2.*re*t/sqrt((1.+ec)**(1.+ec)*(1.-ec)**(1.-ec))
    else:
        sl  = slat*pi/180.
        tc  = tan(pi/4.-sl/2.)/((1.-ec*sin(sl))/(1.+ec*sin(sl)))**(ec/2.)
        mc  = cos(sl)/sqrt(1.0-ec2*(sin(sl)**2))
        rho = re*mc*t/tc

    y=-rho*sgn*cos(sgn*longitude)
    x= rho*sgn*sin(sgn*longitude)
    return [x,y]

def  mapll2(lat, lon):
    if lat>0:
        sgn=1
    else:
        sgn=-1
    return mapll(lat,lon,sgn)

def dist(lat1,lon1,lat2,lon2):
    x1,y1=mapll2(lat1,lon1)
    x2,y2=mapll2(lat2,lon2)
    #print x1,y1,x2,y2
    d=sqrt((x2-x1)**2+(y2-y1)**2)
    if sign(lat1)!=sign(lat2):
        d=10000.0
    return d
