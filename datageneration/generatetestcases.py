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

if __name__ == "__main__":
    A=1
    r=48
    
    expt=10
    widpix=128
    heipix=128

    subfolder="Dec26"

    cases=["constp","const0","constn","lgradpp","lgradpn","islandp0","islandppn","islandpnn","island0p","islandnpp","islandnnp","lgradppn","lgradnnp"]
    
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
    gendata(*lineargrad(-A,A*0.3,widpix,heipix,expt),subfolder,cases[12])