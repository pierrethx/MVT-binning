import numpy as np
import matplotlib.pyplot as plt 
import tkinter
from tkinter.filedialog import askopenfilename
from astropy.io import fits
import functions,bin_accretion,wvt_iteration
import scipy.spatial as sp
import time

import tkinter as tk
from tkinter import simpledialog

def getname(initial=""):
    application_window = tk.Tk()
    application_window.withdraw()
    target = simpledialog.askstring("Input", "What is the objectname?",parent=application_window,initialvalue=initial)
    if target is None:
        raise NameError("No name, goodbye.")
    return target

def gettarget(initial=0):
    application_window = tk.Tk()
    application_window.withdraw()
    target = simpledialog.askinteger("Input", "What is the target Signal-to-Noise? \n For cvt instead of wvt make this number negative :]",parent=application_window,initialvalue=initial)
    if target is None:
        raise NameError("No target, goodbye.")
    return target

def mainfunc(signal,var,target,weighting=True,displayWVT=True,epsilon=10):
    binlist,init_generators,init_scalelengths=bin_accretion.cc_accretion(signal,var,target)
    binlist,diflist=wvt_iteration.iteration_func(target,signal,var,init_generators,init_scalelengths,epsilon,displaywvt=displayWVT)
    wvt,ston=functions.generate_wvt2(binlist,signal,var,displayWVT)
    vwvt=functions.generate_wvt(binlist,var)
    if displayWVT:
        maketargetscatter(target,binlist,signal,var)
    #blockout(target,wvt,ston)
    return binlist,diflist

def saveiteratedfits(target,wcsx,wvt,vwvt,objname,sourcedir,subfolder,weighting=True):
    
    header=wcsx.to_header()
    hdu = fits.PrimaryHDU(np.flipud(wvt),header=header)
    hdul = fits.HDUList([hdu])
    if weighting:
        hdul.writeto(sourcedir+"/"+subfolder+"/"+objname+"_wit_sig.fits",overwrite=True)
    else:
        hdul.writeto(sourcedir+"/"+subfolder+"/"+objname+"_cit_sig.fits",overwrite=True)

    hdu2 = fits.PrimaryHDU(np.flipud(vwvt),header=header)
    hdul2 = fits.HDUList([hdu2])
    if weighting:
        hdul2.writeto(sourcedir+"/"+subfolder+"/"+objname+"_wit_var.fits",overwrite=True)
    else:
        hdul2.writeto(sourcedir+"/"+subfolder+"/"+objname+"_cit_var.fits",overwrite=True)

def saveblockoutfits(target,ston,wcsx,wvt,vwvt,objname,sourcedir,subfolder,weighting=True):
    blockout(target,wvt,ston)
    header=wcsx.to_header()
    hdu = fits.PrimaryHDU(np.flipud(wvt),header=header)
    hdul = fits.HDUList([hdu])
    if weighting:
        hdul.writeto(sourcedir+"/"+subfolder+"/block_"+objname+"_wit_sig.fits",overwrite=True)
    else:
        hdul.writeto(sourcedir+"/"+subfolder+"/block_"+objname+"_cit_sig.fits",overwrite=True)

    hdu2 = fits.PrimaryHDU(np.flipud(vwvt),header=header)
    hdul2 = fits.HDUList([hdu2])
    if weighting:
        hdul2.writeto(sourcedir+"/"+subfolder+"/block_"+objname+"_wit_var.fits",overwrite=True)
    else:
        hdul2.writeto(sourcedir+"/"+subfolder+"/block_"+objname+"_cit_var.fits",overwrite=True)

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

def blockout(target,wvt,ston):
    for y in range(len(ston)):
        for x in range(len(ston[y])):
            if ston[y][x]<=target*0.5:
                wvt[y][x]=0
                #ston[y][x]=0

def saveston(wcsx,ston,sourcedir,objname,subfolder="unbinned"):
    header=wcsx.to_header()
    hdu3 = fits.PrimaryHDU(np.flipud(ston),header=header)
    hdul3 = fits.HDUList([hdu3])
    hdul3.writeto(sourcedir+"/"+subfolder+"/zston_"+objname+".fits",overwrite=True)

def saveassign(wcsx,assign,sourcedir,objname,subfolder="unbinned"):
    header=wcsx.to_header()
    hdu3 = fits.PrimaryHDU(np.flipud(assign),header=header)
    hdul3 = fits.HDUList([hdu3])
    hdul3.writeto(sourcedir+"/"+subfolder+"/z_"+objname+"_assigned.fits",overwrite=True)

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
    
    wvt,ston=functions.generate_wvt4(binlist,signal,var,np.full(len(binlist),1),True)
    wvt,ston=functions.generate_wvt3(binlist,signal,var,np.full(len(binlist),1),True)
    
    #saveiteratedfits(target,binlist,wcsx,signal,var,objname,sourcedir,subfolder="target"+str(target))
    #saveblockoutfits(target,binlist,wcsx,signal,var,objname,sourcedir,subfolder="target"+str(target))
    #saveston(wcsx,signal,var,sourcedir,objname+"_wit",subfolder="target"+str(target))