import numpy as np
import matplotlib.pyplot as plt 
import tkinter
from tkinter.filedialog import askopenfilename
from astropy.io import fits
import bin_accretion,wvt_iteration

#175x130

widpix=120
heipix=90

xlist=np.arange(0,widpix,1)
ylist=np.arange(0,heipix,1)
xx,yy=np.meshgrid(xlist,ylist)


## Assume that intensity follows I=A exp(-(x^2+y^2)/a^2) with x,y dist from center
centerx=int(widpix/2)
centery=int(heipix/2)
A=20
a=10
signal=A*np.exp(-((xx-centerx)**2+(yy-centery)**2)/(a**2))

## assume an approximately constant noise, with amplitude B, and range b
B=10
b=1
var=B+b*(np.random.rand(heipix,widpix)-0.5)
sourcedir="/Users/pierre/Downloads"
objname="testdata1"
target=3 #Target Signal To noise
wvt,init_generators, init_scalelengths=bin_accretion.cc_accretion(signal,var,target)
epsilon=10
wvt=wvt_iteration.iteration_func(target,signal,var,init_generators,init_scalelengths,epsilon,displaywvt=True)
hdu = fits.PrimaryHDU(wvt)
hdul = fits.HDUList([hdu])
hdul.writeto(sourcedir+"/iterated_"+objname+"_target_"+str(target)+".fits",overwrite=True)
fig,ax=plt.subplots()
image=ax.imshow(wvt,cmap="cubehelix")
fig.colorbar(image)
plt.show()



