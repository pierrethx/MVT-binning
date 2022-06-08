import numpy as np
import matplotlib.pyplot as plt 
import tkinter
from tkinter.filedialog import askopenfilename
from astropy.io import fits
#import bin_accretion,wvt_iteration
import os

#destination="/Users/pierre/Downloads/"
destination="/Volumes/TOSHIBA EXT/MARTIN2/"

def makeitfolder(subfolder,text,expt):
    slist=subfolder.split("/")
    for s in range(len(slist)):
        try:
            os.mkdir(destination+"/".join(slist[:s+1]))
        except:
            pass
    for ex in expt:
        try:
            os.mkdir(destination+subfolder+"/"+text+ex)
        except:
            pass
        targs=["unbinned","target3","target5","target8","target10","target15"]
        for t in targs:
            try:
                os.mkdir(destination+subfolder+"/"+text+ex  +"/"+t)
            except:
                pass
                print("already exists")
    return [destination+subfolder+"/"+text+ex for ex in expt]

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
    print("A: "+str(A))
    print("r c: "+str(r))
    print("b: "+str(b))
    print("bg: "+str(I_b))
    print("min: "+str(np.min(model)))
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

def gendata2(widpix,heipix,A,r,b,I_b,expt,edge,sourcedir,objname):
    

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
    hdul.writeto(sourcedir+"/"+objname+"_rc"+str(r)+"_signal.fits",overwrite=True)
    hdu = fits.PrimaryHDU(var)
    hdul = fits.HDUList([hdu])
    hdul.writeto(sourcedir+"/"+objname+"_rc"+str(r)+"_var.fits",overwrite=True)

def gendatam(widpix,heipix,A,r,b,I_b,expt,edge,sourcedir,objname):
    
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

    '''
    fig,ax = plt.subplots(1,1)
    ax.imshow(signal/np.sqrt(var),cmap="cubehelix")
    plt.show()
    '''

    ## A crispy little layer of noises to stop evil twins!!!! (When a 2px bin has StoN 0 because of equal but opposite StoNs)

    #if edge<0.5*(np.sqrt((widpix**2+heipix**2))):
    hdu = fits.PrimaryHDU(model)
    hdul = fits.HDUList([hdu])
    hdul.writeto(sourcedir+"/"+objname+"_edge"+str(edge)+"_model.fits",overwrite=True)

## generate variance, since all pixels are given by poisson dis, the variances=lambda

def genmoredata(num,widpix,heipix,A,r,b,I_b,expt,edge,subfolder):
    ## Assume that lambda follows I=A [1+(x^2+y^2)/r^2)]**(0.5-3b)with x,y dist from center
    if num==1:
        for back in I_b:
            for ex in expt:
                for ed in edge:
                    #sourcedir="/Users/pierre/Downloads/"+subfolder+"/"+str(widpix)+"x"+str(heipix)+"_peak"+str(A)+"/bg"+str(back)+"/unbinned"
                    objname="testdata"
                    sourcedir=destination+subfolder+"/exp"+str(int(ex))+"/unbinned"

                    gendata(widpix,heipix,A,r,b,back,int(ex),ed,sourcedir,objname)
    else:
        for m in range(num):
            for back in I_b:
                for ed in edge:
                    sourcedir=destination+subfolder+"/"+str(widpix)+"x"+str(heipix)+"_peak"+str(A)+"/bg"+str(back)+"/unbinned"
                    objname="testdata"+str(m+1)

                    gendata(widpix,heipix,A,r,b,back,expt,ed,sourcedir,objname)

def genmoredata2(num,widpix,heipix,A,r,b,I_b,expt,edge,subfolder):
    ## Assume that lambda follows I=A [1+(x^2+y^2)/r^2)]**(0.5-3b)with x,y dist from center
    for a in range(len(A)):
            for ed in edge:
                #sourcedir="/Users/pierre/Downloads/"+subfolder+"/"+str(widpix)+"x"+str(heipix)+"_peak"+str(A)+"/bg"+str(back)+"/unbinned"
                objname="testdata"
                sourcedir=destination+subfolder+"/sig"+str(int(A[a]))+"/unbinned"

                gendata(widpix,heipix,A[a],r,b,I_b[a],int(expt[a]),ed,sourcedir,objname)
 
def genlotsdata(num,widpix,heipix,A,r,b,I_b,expt,edge,subfolder,text):
    for i in range( num):
        for ed in edge:
            makeitfolder(subfolder,text,[str(round(be,2)) for be in b])
            for be in b:
                sourcedir=destination+subfolder+"/"+text+str(round(be,2))+"/unbinned"
                objname="testdata"+str(i)
                gendata(widpix,heipix,A,r,be,I_b,int(expt),ed,sourcedir,objname)

def genlotsdata2(num,widpix,heipix,A,r,b,I_b,expt,edge,subfolder,text):
    for i in range( num):
        for rc in r:
            makeitfolder(subfolder,text,[str(round(be,2)) for be in b])
            for be in b:
                sourcedir=destination+subfolder+"/"+text+str(round(be,2))+"/unbinned"
                objname="testdata"+str(i)
                gendata2(widpix,heipix,A,rc,be,I_b,int(expt),edge,sourcedir,objname)

def genlotssigdata(num,widpix,heipix,A,r,b,I_b,expt,edge,subfolder,text):
    for i in range( num):
        for ed in edge:
            for j in range(len(A)):
                makeitfolder(subfolder,text,[str(Aj) for Aj in A])
                sourcedir=destination+subfolder+"/"+text+str(A[j])+"/unbinned"
                objname="testdata"+str(i)
                gendata(widpix,heipix,A[j],r,b,I_b[j],int(expt),ed,sourcedir,objname)

def genlotssigdatax(num,widpix,heipix,A,r,b,I_b,expt,edge,subfolder,text):
    for i in range( num):
        for ed in edge:
            for j in range(len(A)):
                makeitfolder(subfolder,text,[str(Aj) for Aj in A])
                sourcedir=destination+subfolder+"/"+text+str(A[j])+"/unbinned"
                objname="testdata"+str(i)
                gendata(widpix,heipix,A[j],r,b,I_b[j],int(expt[j]),ed,sourcedir,objname)

def genmodel(widpix,heipix,A,r,b,edge,subfolder,text):
    for ed in edge:
        for j in range(len(A)):
            sourcedir=destination+subfolder+"/"+text+str(A[j])+"/unbinned"
            gendatam(widpix,heipix,A[j],r,b,I_b[j],int(expt[j]),ed,sourcedir,"testdata")

if __name__ == "__main__":
    #A=[1,10,100,1000,10000]
    A=10
    #r=16
    r=[4,8,12,16,20]
    #b=0.67
    b=[0.67,1.17,1.67,2.17,2.67,3.17]
    #I_b=[0.3,3,30,300,3000] #background signal
    I_b=3
    #edge=[25,50,75,100]
    #edge=[50,75,100]
    #edge=[75]
    edge=75
    #expt=[3e4,1e5,3e5,1e6]
    #expt=[1e2,1e4,1e6]
    expt=10000
    #expt=[100000,10000,1000,100,10]

    widpix=160
    heipix=160

    subfolder="simdata/betatest"

    #makeitfolder(subfolder,"sig",A)
    #num=1
    #genmoredata2(num,widpix,heipix,A,r,b,I_b,expt,edge,subfolder)
    #genmoredata(20,widpix,heipix,A,r,b,I_b,expt,edge,subfolder)
    genlotsdata2(1,widpix,heipix,A,r,b,I_b,expt,edge,subfolder,"beta")
    #genlotssigdatax(10,widpix,heipix,A,r,b,I_b,expt,edge,subfolder,"sig")
    genmodel(widpix,heipix,A,r,b,edge,subfolder,"sig")