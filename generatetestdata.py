import numpy as np
import matplotlib.pyplot as plt 
import tkinter
from tkinter.filedialog import askopenfilename
from astropy.io import fits
import bin_accretion,wvt_iteration

def gendata(widpix,heipix,A,r,b,I_b,edge,sourcedir,objname):
    xcent=int(widpix/2)
    ycent=int(heipix/2)

    xlist=np.arange(0,widpix,1)
    ylist=np.arange(0,heipix,1)
    xx,yy=np.meshgrid(xlist,ylist)

    model=A*((1+((xx-xcent)**2+(yy-ycent)**2)/r**2))**(0.5-3*b)
    for j in range(len(xx)):
        for i in range(len(xx[0])):
            if ((xx[j][i]-xcent)**2+(yy[j][i]-ycent)**2)>edge**2:
                model[j][i]=0
    model+=I_b
    model_w_noise=np.random.poisson(lam=model)+0.01*np.random.rand(len(model),len(model[0]))
    #+noise2*np.random.rand(len(xx),len(xx[0]))
    noise=model-model_w_noise
    var=model_w_noise
    signal=model_w_noise - I_b

    '''
    fig,ax = plt.subplots(1,1)
    ax.imshow(signal/np.sqrt(var),cmap='cubehelix')
    plt.show()
    '''
    

    ## A crispy little layer of noises to stop evil twins!!!! (When a 2px bin has StoN 0 because of equal but opposite StoNs)

    #if edge<0.5*(np.sqrt((widpix**2+heipix**2))):
    hdu = fits.PrimaryHDU(signal)
    hdul = fits.HDUList([hdu])
    hdul.writeto(sourcedir+"/"+objname+"_edge"+str(edge)+"_signal.fits",overwrite=True)
    hdu = fits.PrimaryHDU(var)
    hdul = fits.HDUList([hdu])
    hdul.writeto(sourcedir+"/"+objname+"_edge"+str(edge)+"_var.fits",overwrite=True)

## generate variance, since all pixels are given by poisson dis, the variances=lambda

## Assume that lambda follows I=A [1+(x^2+y^2)/r^2)]**(0.5-3b)with x,y dist from center
A=100
r=32
b=0.67
I_b=[20,40,60,80,100,120] #background signal
edge=[35,70,105]

widpix=128
heipix=128

for back in I_b:
    for ed in edge:
        sourcedir="/Users/pierre/Downloads/Aug4/"+str(widpix)+"x"+str(heipix)+"_peak"+str(A)+"/bg"+str(back)+"/"
        objname="testdata"

        gendata(widpix,heipix,A,r,b,back,ed,sourcedir,objname)

