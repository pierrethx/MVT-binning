import numpy as np
import matplotlib.pyplot as plt 
import tkinter
from tkinter.filedialog import askopenfilename
from astropy.io import fits
import bin_accretion,wvt_iteration

#175x130

widpix=128
heipix=128
xcent=int(widpix/2)
ycent=int(heipix/2)

xlist=np.arange(0,widpix,1)
ylist=np.arange(0,heipix,1)
xx,yy=np.meshgrid(xlist,ylist)

## generate variance, since all pixels are given by poisson dis, the variances=lambda

## Assume that lambda follows I=A [1+(x^2+y^2)/r^2)]**(0.5-3b)with x,y dist from center
A=50
r=32
b=0.67
I_b=100 #background signal
edge=80
noise2=1

var=A*((1+((xx-xcent)**2+(yy-ycent)**2)/r**2))**(0.5-3*b)
#var[((xx-xcent)**2+(yy-ycent)**2)>edge**2]=0

for j in range(len(xx)):
    for i in range(len(xx[0])):
        if ((xx[j][i]-xcent)**2+(yy[j][i]-ycent)**2)>edge**2:
            var[j][i]=0
var+= noise2*(np.random.rand(len(var),len(var[0])))
if edge<0.5*(np.sqrt((widpix**2+heipix**2))):
    sourcedir="/Users/pierre/Downloads/July30/size"+str(widpix)+"x"+str(heipix)+"/noise"+str(I_b)+"/edge"+str(edge)
else:
    sourcedir="/Users/pierre/Downloads/July30/size"+str(widpix)+"x"+str(heipix)+"/noise"+str(I_b)+"/noedge"
objname="testdata"+str(I_b)
print(sourcedir)
lamba=var+I_b

## then the signal is a poisson distribution with lambda

signal=np.random.poisson(lam=lamba)-I_b

hdu = fits.PrimaryHDU(signal)
hdul = fits.HDUList([hdu])
hdul.writeto(sourcedir+"/"+objname+"_signal.fits",overwrite=True)
hdu = fits.PrimaryHDU(var)
hdul = fits.HDUList([hdu])
hdul.writeto(sourcedir+"/"+objname+"_var.fits",overwrite=True)




