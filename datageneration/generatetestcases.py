import numpy as np
import matplotlib.pyplot as plt 
import tkinter,os
from tkinter.filedialog import askopenfilename
from astropy.io import fits
import bin_accretion,wvt_iteration,generatetestdata


## the test cases are:
## constant signal (+,0,-)
## linear gradient (++,+-)
## islands (+0,++-,+--,0+,-++,--+)

def makeitfolder(subfolder,expt):
    slist=subfolder.split("/")
    for s in range(len(slist)):
        try:
            os.mkdir("/Users/pierre/Downloads/"+"/".join(slist[:s+1]))
        except:
            pass
    try:
        os.mkdir("/Users/pierre/Downloads/"+subfolder+"/"+expt)
    except:
        pass
    targs=["unbinned","target3","target5","target10"]
    for t in targs:
        try:
            os.mkdir("/Users/pierre/Downloads/"+subfolder+"/"+expt+"/"+t)
        except:
            pass
    print("already exists")

def gendata(signal,var,subfolder,case):

    boof="testcases"

    makeitfolder(subfolder,boof)
    
    objname="testdata_"+case
    sourcedir="/Users/pierre/Downloads/"+subfolder+"/"+boof+"/unbinned"


    #if edge<0.5*(np.sqrt((widpix**2+heipix**2))):
    hdu = fits.PrimaryHDU(signal)
    hdul = fits.HDUList([hdu])
    hdul.writeto(sourcedir+"/"+objname+"_signal.fits",overwrite=True)
    hdu = fits.PrimaryHDU(var)
    hdul = fits.HDUList([hdu])
    hdul.writeto(sourcedir+"/"+objname+"_var.fits",overwrite=True)
    
## the test cases are:
## constant signal (+,0,-)
## linear gradient (++,+-)
## islands (+0,++-,+--,0+,-++,--+)

def lineargrad(left,right,widpix,heipix,expt):
    xlist=np.arange(0,widpix,1)
    ylist=np.arange(0,heipix,1)
    xx,yy=np.meshgrid(xlist,ylist)

    lazx=widpix-1
    lingrad=right*xx/lazx+left*(1-xx/lazx)
    return lingrad,np.abs(lingrad)/expt

def x3grad(left,right,widpix,heipix,expt):
    xlist=np.arange(0,widpix,1)
    ylist=np.arange(0,heipix,1)
    xx,yy=np.meshgrid(xlist,ylist)

    xl=np.cbrt(left)
    xr=np.cbrt(right)

    lazx=widpix-1
    lingrad=xr*xx/lazx+xl*(1-xx/lazx)
    lin3g=lingrad**3
    return lin3g,np.abs(lin3g)/expt

def hlineargrad(right,widpix,heipix,expt):
    xlist=np.arange(0,widpix,1)
    ylist=np.arange(0,heipix,1)
    xx,yy=np.meshgrid(xlist,ylist)

    lazx=widpix-1
    lingrad=right*xx/lazx-right*(1-xx/lazx)
    lingrad[xx<(widpix)/2]=0
    return lingrad,np.abs(lingrad)/expt

def hx3grad(right,widpix,heipix,expt):
    xlist=np.arange(0,widpix,1)
    ylist=np.arange(0,heipix,1)
    xx,yy=np.meshgrid(xlist,ylist)

    xr=np.cbrt(right)

    lazx=widpix-1
    lingrad=xr*xx/lazx-xr*(1-xx/lazx)
    lin3g=lingrad**3
    lin3g[xx<(widpix)/2]=0
    return lin3g,np.abs(lin3g)/expt

def hx2grad(right,widpix,heipix,expt):
    xlist=np.arange(0,widpix,1)
    ylist=np.arange(0,heipix,1)
    xx,yy=np.meshgrid(xlist,ylist)

    xr=np.sqrt(right)

    lazx=widpix-1
    lingrad=xr*xx/lazx-xr*(1-xx/lazx)
    lin2g=lingrad**2
    lin2g[xx<(widpix)/2]=0
    return lin2g,np.abs(lin2g)/expt

def island(island,ocean,radius,widpix,heipix,expt):
    xcent=int(widpix/2.0)
    ycent=int(heipix/2.0)
    print(xcent)

    xlist=np.arange(0,widpix,1.0)
    ylist=np.arange(0,heipix,1.0)
    xx,yy=np.meshgrid(xlist,ylist)

    ocean=0*xx+ocean
    island=0*xx+island
    
    rrsq=(xx-xcent)**2+(yy-ycent)**2
    for y in range(len(rrsq)):
        for x in range(len(rrsq[0])):
            if rrsq[y][x]<=radius**2:
                rrsq[y][x]=island[y][x]
            else:
                rrsq[y][x]=ocean[y][x]

    return rrsq,np.abs(rrsq)/expt

def wnoise(model,I_b,exptime):
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

if __name__ == "__main__":
    A=5
    r=48

    bg=5
    
    expt=50
    widpix=128
    heipix=128

    subfolder="Mar4"

    cases=["constp","const0","constn","lgradpp","lgradpn","islandp0","islandppn","islandpnn","island0p","islandnpp","islandnnp","lgradppn","lgradnnp"]
    '''
    gendata(*lineargrad(A,A,widpix,heipix,expt),subfolder,cases[0])
    gendata(*lineargrad(0,0,widpix,heipix,expt),subfolder,cases[1])
    gendata(*lineargrad(-A,-A,widpix,heipix,expt),subfolder,cases[2])

    gendata(*lineargrad(A/10,A,widpix,heipix,expt),subfolder,cases[3])
    gendata(*lineargrad(A,-A,widpix,heipix,expt),subfolder,cases[4])

    gendata(*island(A,0,r,widpix,heipix,expt),subfolder,cases[5])
    gendata(*island(A,-A*0.5,r,widpix,heipix,expt),subfolder,cases[6])
    gendata(*island(A,-2*A,r,widpix,heipix,expt),subfolder,cases[7])

    gendata(*island(0,A,r,widpix,heipix,expt),subfolder,cases[8])
    gendata(*island(-A,2*A,r,widpix,heipix,expt),subfolder,cases[9])
    gendata(*island(-A,A*0.5,r,widpix,heipix,expt),subfolder,cases[10])
    
    gendata(*lineargrad(-A*0.3,A,widpix,heipix,expt),subfolder,cases[11])
    gendata(*lineargrad(-A,A*0.3,widpix,heipix,expt),subfolder,cases[12])'''

    gendata(*hlineargrad(A,widpix,heipix,expt),subfolder,"c_hLgrad")
    gendata(*wnoise(hlineargrad(A,widpix,heipix,expt)[0],bg,1000*expt),subfolder,"d_hLgrad")

    #gendata(*hx3grad(A,widpix,heipix,expt),subfolder,"c_h3grad")
    #gendata(*wnoise(hx3grad(A,widpix,heipix,expt)[0],bg,1000*expt),subfolder,"d_h3grad")

    #gendata(*hx2grad(A,widpix,heipix,expt),subfolder,"c_h2grad")
    #gendata(*wnoise(hx2grad(A,widpix,heipix,expt)[0],bg,1000*expt),subfolder,"d_h2grad")