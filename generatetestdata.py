import numpy as np
import matplotlib.pyplot as plt 
import tkinter
from tkinter.filedialog import askopenfilename
from astropy.io import fits
import bin_accretion,wvt_iteration
import os

def makeitfolder(subfolder,expt):
    slist=subfolder.split("/")
    for s in range(len(slist)):
        try:
            os.mkdir("/Users/pierre/Downloads/"+"/".join(slist[:s+1]))
        except:
            pass
    for ex in expt:
        try:
            os.mkdir("/Users/pierre/Downloads/"+subfolder+"/exp"+str(int(ex)))
        except:
            pass
        targs=["unbinned","target3","target5","target10"]
        for t in targs:
            try:
                os.mkdir("/Users/pierre/Downloads/"+subfolder+"/exp"+str(int(ex))  +"/"+t)
            except:
                pass
    print("already exists")

def generator(widpix,heipix,A,r,b,I_b,expt,edge):
    xcent=int(widpix/2)
    ycent=int(heipix/2)

    xlist=np.arange(0,widpix,1)
    ylist=np.arange(0,heipix,1)
    xx,yy=np.meshgrid(xlist,ylist)

    model=A*((1+((xx-xcent)**2+(yy-ycent)**2)/r**2))**(0.5-3*b)
    if np.isnan(edge):
        pass
    else:
        for j in range(len(xx)):
            for i in range(len(xx[0])):
                if ((xx[j][i]-xcent)**2+(yy[j][i]-ycent)**2)>edge**2:
                    model[j][i]=0
    model+=I_b
    model*=expt
    model_w_noise=(np.random.poisson(lam=model)+0.01*np.random.rand(len(model),len(model[0])))
    #+noise2*np.random.rand(len(xx),len(xx[0]))
    noise=model-model_w_noise
    var=model_w_noise/expt/expt
    #signal=model_w_noise - I_b
    #signal=model_w_noise/expt-0.9*I_b
    signal=model_w_noise/expt-I_b

    return signal,var

def ogenerator(widpix,heipix,center,A,r,b,I_b,expt,edge):
    ycent=center[0]
    xcent=center[1]

    xlist=np.arange(0,widpix,1)
    ylist=np.arange(0,heipix,1)
    xx,yy=np.meshgrid(xlist,ylist)

    model=A*((1+((xx-xcent)**2+(yy-ycent)**2)/r**2))**(0.5-3*b)
    if np.isnan(edge):
        pass
    else:
        for j in range(len(xx)):
            for i in range(len(xx[0])):
                if ((xx[j][i]-xcent)**2+(yy[j][i]-ycent)**2)>edge**2:
                    model[j][i]=0
    model+=I_b
    model*=expt
    model_w_noise=(np.random.poisson(lam=model)+0.01*np.random.rand(len(model),len(model[0])))
    #+noise2*np.random.rand(len(xx),len(xx[0]))
    noise=model-model_w_noise
    var=model_w_noise/expt/expt
    #signal=model_w_noise - I_b
    #signal=model_w_noise/expt-0.9*I_b
    signal=model_w_noise/expt-I_b

    return signal,var

def gendata(widpix,heipix,A,r,b,I_b,expt,edge,sourcedir,objname):
    

    signal,var=generator(widpix,heipix,A,r,b,I_b,expt,edge)
    '''
    fig,ax = plt.subplots(1,1)
    ax.imshow(signal/np.sqrt(var),cmap="cubehelix")
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

def genmoredata(num,widpix,heipix,A,r,b,I_b,expt,edge,subfolder):
    ## Assume that lambda follows I=A [1+(x^2+y^2)/r^2)]**(0.5-3b)with x,y dist from center
    if num==1:
        for back in I_b:
            for ex in expt:
                for ed in edge:
                    #sourcedir="/Users/pierre/Downloads/"+subfolder+"/"+str(widpix)+"x"+str(heipix)+"_peak"+str(A)+"/bg"+str(back)+"/unbinned"
                    objname="testdata"
                    sourcedir="/Users/pierre/Downloads/"+subfolder+"/exp"+str(int(ex))+"/unbinned"

                    gendata(widpix,heipix,A,r,b,back,int(ex),ed,sourcedir,objname)
    else:
        for m in range(num):
            for back in I_b:
                for ed in edge:
                    sourcedir="/Users/pierre/Downloads/"+subfolder+"/"+str(widpix)+"x"+str(heipix)+"_peak"+str(A)+"/bg"+str(back)+"/unbinned"
                    objname="testdata"+str(m+1)

                    gendata(widpix,heipix,A,r,b,back,expt,ed,sourcedir,objname)

if __name__ == "__main__":
    A=10.75
    r=20
    b=1.8
    I_b=[10] #background signal
    #edge=[25,50,75,100]
    edge=[50,75,100]
    #expt=[3e4,1e5,3e5,1e6]
    expt=[1e2,1e4,1e6]

    widpix=128
    heipix=128

    subfolder="Nov6/expansion"

    makeitfolder(subfolder,expt)
    num=1
    genmoredata(num,widpix,heipix,A,r,b,I_b,expt,edge,subfolder)