import numpy as np
import matplotlib.pyplot as plt 
import tkinter
from tkinter.filedialog import askopenfilename
from astropy.io import fits
from astropy import wcs
import functions,bin_accretion,wvt_iteration,qradial
import scipy.spatial as sp
import time

import tkinter as tk
from tkinter import simpledialog

def getname(initial=""):
    ## for getting the user to name an object
    application_window = tk.Tk()
    application_window.withdraw()
    target = simpledialog.askstring("Input", "What is the objectname?",parent=application_window,initialvalue=initial)
    if target is None:
        raise NameError("No name, goodbye.")
    return target

def gettarget(initial=0):
    ## for getting the target StoN from user
    application_window = tk.Tk()
    application_window.withdraw()
    target = simpledialog.askinteger("Input", "What is the target Signal-to-Noise? \n For cvt instead of wvt make this number negative :]",parent=application_window,initialvalue=initial)
    if target is None:
        raise NameError("No target, goodbye.")
    return target

def trim(lis,targ):
    other=[]
    for i in range(len(lis)):
        if lis[i]<=targ:
            other.append(lis[i])
    return other

def mainfunc(signal,var,target,weighting=True,displayWVT=True,epsilon=10):
    ## this is the most important function
    ## first we do bin accretion. then iteration. the wvt and vwvt is not necessary here but im scared to remove it
    binlist,init_generators,init_scalelengths=bin_accretion.cc_accretion(signal,var,target)
    binlist,diflist=wvt_iteration.iteration_func(target,signal,var,init_generators,init_scalelengths,epsilon,displaywvt=displayWVT)
    
    if displayWVT:
        wvt,ston=functions.generate_wvt2(binlist,signal,var,displayWVT)
        vwvt=functions.generate_wvt(binlist,var)
        maketargetscatter(target,binlist,signal,var)
    #blockout(target,wvt,ston)
    return binlist,diflist

def saveiteratedfits(target,wcsx,wvt,vwvt,objname,sourcedir,subfolder,weighting=True):
    
    hdu = fits.PrimaryHDU(np.flipud(wvt),header=wcsx)
    hdul = fits.HDUList([hdu])
    #manipulate(hdul)
    if weighting:
        hdul.writeto(sourcedir+"/"+subfolder+"/"+objname+"_wit_sig.fits",overwrite=True,checksum=True)
    else:
        hdul.writeto(sourcedir+"/"+subfolder+"/"+objname+"_cit_sig.fits",overwrite=True,checksum=True)

    hdu2 = fits.PrimaryHDU(np.flipud(vwvt),header=wcsx)
    hdul2 = fits.HDUList([hdu2])
    #4manipulate(hdul2)
    if weighting:
        hdul2.writeto(sourcedir+"/"+subfolder+"/"+objname+"_wit_var.fits",overwrite=True,checksum=True)
    else:
        hdul2.writeto(sourcedir+"/"+subfolder+"/"+objname+"_cit_var.fits",overwrite=True,checksum=True)

def saveblockoutfits(target,ston,wcsx,wvt,vwvt,objname,sourcedir,subfolder,weighting=True):
    blockout(target,wvt,ston)

    ## for visualizing results
    gig,gax=plt.subplots()
    edge=qradial.endnumber(objname)
    gax.imshow(ston,cmap="cubehelix")
    center=(len(ston[0])/2-0.5,len(ston)/2-0.5)
    gax.plot([center[0],center[0]],[0,len(ston)],color="red")
    plt.show()
    plt.savefig(sourcedir+"/"+subfolder+"/block_"+objname+"_stonover.png")
    plt.close()

    hdu = fits.PrimaryHDU(np.flipud(wvt),header=wcsx)
    hdul = fits.HDUList([hdu])
    #manipulate(hdul)
    if weighting:
        hdul.writeto(sourcedir+"/"+subfolder+"/block_"+objname+"_wit_sig.fits",overwrite=True,checksum=True)
    else:
        hdul.writeto(sourcedir+"/"+subfolder+"/block_"+objname+"_cit_sig.fits",overwrite=True,checksum=True)

    hdu2 = fits.PrimaryHDU(np.flipud(vwvt),header=wcsx)
    hdul2 = fits.HDUList([hdu2])
    #manipulate(hdul2)
    if weighting:
        hdul2.writeto(sourcedir+"/"+subfolder+"/block_"+objname+"_wit_var.fits",overwrite=True,checksum=True)
    else:
        hdul2.writeto(sourcedir+"/"+subfolder+"/block_"+objname+"_cit_var.fits",overwrite=True,checksum=True)

def saveblockoutoldfits(target,ston,wcsx,wvt,vwvt,objname,sourcedir,subfolder,weighting=True):
    blockout_old(target,wvt,ston)
    hdu = fits.PrimaryHDU(np.flipud(wvt),header=wcsx)
    hdul = fits.HDUList([hdu])
    #manipulate(hdul)
    if weighting:
        hdul.writeto(sourcedir+"/"+subfolder+"/oblock_"+objname+"_wit_sig.fits",overwrite=True,checksum=True)
    else:
        hdul.writeto(sourcedir+"/"+subfolder+"/oblock_"+objname+"_cit_sig.fits",overwrite=True,checksum=True)

    hdu2 = fits.PrimaryHDU(np.flipud(vwvt),header=wcsx)
    hdul2 = fits.HDUList([hdu2])
    #manipulate(hdul2)
    if weighting:
        hdul2.writeto(sourcedir+"/"+subfolder+"/oblock_"+objname+"_wit_var.fits",overwrite=True,checksum=True)
    else:
        hdul2.writeto(sourcedir+"/"+subfolder+"/oblock_"+objname+"_cit_var.fits",overwrite=True,checksum=True)


def maketargetscatter(target,binlist,signal,var):
    fig,ax=plt.subplots()
    ax.set_xlabel("generator radius from center")
    ax.set_ylabel("bin signal to noise")
    centx=0
    centy=0
    totalm=0
    gcents=[]
    grads=[]
    gstons=[]
    for binn in binlist:
        gstons.append(functions.calculate_SN(binn,signal,var))
        gcent,gmass=functions.weighted_centroid(binn,signal)
        gcent=functions.geometric_center(binn)
        centx+=gcent[1]*gmass
        centy+=gcent[0]*gmass
        totalm+=gmass
        gcents.append(gcent)
    centx=centx/totalm
    centy=centy/totalm
    for cent in gcents:
        grads.append(np.sqrt((cent[0]-centy)**2+(cent[1]-centx)**2))
    ax.plot([0,np.nanmax(grads)],[target,target],linestyle="dashed",color="navy")
    ax.plot(grads,gstons,linewidth=0,marker="s")
    plt.show()

def blockout_old(target,wvt,ston):
    for y in range(len(ston)):
        for x in range(len(ston[y])):
            if ston[y][x]<=target*0.5:
                wvt[y][x]=0
                ston[y][x]=0

def blockout(target,wvt,ston):
    ybin=[]
    continuous=False
    zbg=True
    for y in range(len(ston)):
        for x in range(len(ston[y])):
            ybin.append(ston[y][x])
    binsz=100
    fig,ax=plt.subplots()
    negacc=0
    nn,bins,patches=ax.hist(trim(ybin,100),binsz,alpha=0.5)
    if not continuous:
        minb=np.argmax(nn)
        for c in range(len(nn)):
            if bins[c]>=target:
                break
            elif bins[c+1]<=0:
                negacc+=nn[c]
            else:
                if zbg and negacc>0:
                    negacc*=0.95
                    negacc-=nn[c]
                elif nn[c]<nn[minb]:
                    minb=c
        if minb+2>=len(nn) or (bins[minb+2]>=target and nn[minb+2]<=nn[minb]):
            while bins[minb]>0 and nn[minb-1]>nn[minb]:
                minb-=1
        ## this is a minimum
        cutoff=0.5*(bins[minb]+bins[minb+1])
        ax.plot([target,target],[0,2000],linestyle="dashed",color="green")
        ax.plot([cutoff,cutoff],[0,2000],linestyle="dashed",color="black")
        for y in range(len(ston)):
            for x in range(len(ston[y])):
                if ston[y][x]<cutoff:
                    wvt[y][x]=0
                    ston[y][x]=0
        #plt.show()
        #plt.close()
    return wvt,ston


def saveston(wcsx,ston,sourcedir,objname,subfolder="unbinned"):
    hdu3 = fits.PrimaryHDU(np.flipud(ston),header=wcsx)
    hdul3 = fits.HDUList([hdu3])
    #manipulate(hdul3)
    hdul3.writeto(sourcedir+"/"+subfolder+"/zston_"+objname+".fits",overwrite=True,checksum=True)

def saveassign(wcsx,assign,sourcedir,objname,subfolder="unbinned"):
    hdu3 = fits.PrimaryHDU(np.flipud(assign),header=wcsx)
    hdul3 = fits.HDUList([hdu3])
    #manipulate(hdul3)
    hdul3.writeto(sourcedir+"/"+subfolder+"/z_"+objname+"_assigned.fits",overwrite=True,checksum=True)

if __name__ == "__main__":
    wcsx,signal,var,sourcedir,objname=bin_accretion.initialize(enternew=True)

    signal2=np.copy(signal)
    var2=np.copy(var)
    #var2[signal2<=0]=1e10
    signal2[signal2<=0]=0
    
    

    objname=getname("_".join(objname.split("_")[:-1]))
    sourcedir="/".join(sourcedir.split("/")[:-1])
    #saveston(wcsx,signal,var,sourcedir,objname,subfolder="unbinned")
    #objname="J024815-081723"
    ##empty for directly into sourcedir
    target=gettarget()
    
    if target<0:
        weighting=False
        target=-target
    else:
        weighting=True
    binlist,diflist=mainfunc(signal2,var2,target,displayWVT=False,epsilon=-10)
    
    wvt,ston=functions.generate_wvt4(binlist,signal,var,np.full(len(binlist),1),False)
    wvt,ston=functions.generate_wvt3(binlist,signal,var,np.full(len(binlist),1),False)
    
    subfolder="target"+str(target)

    vwvt=functions.generate_wvt(binlist,var)
    saveiteratedfits(target,wcsx,wvt,vwvt,objname,sourcedir,subfolder=subfolder)
    saveblockoutfits(target,ston,wcsx,wvt,vwvt,objname,sourcedir,subfolder=subfolder)
    saveston(wcsx,ston,sourcedir,objname,subfolder=subfolder)
    assign=functions.assign(binlist,target,ston)
    saveassign(wcsx,assign,sourcedir,objname,subfolder=subfolder)