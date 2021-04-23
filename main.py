import numpy as np
import matplotlib.pyplot as plt 
import tkinter
from tkinter.filedialog import askopenfilename
from astropy.io import fits
from astropy import wcs
import functions,bin_accretion,wvt_iteration,qradial
import scipy.spatial as sp
from scipy import ndimage
import time, os

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
    target = simpledialog.askfloat("Input", "What is the target Signal-to-Noise?",parent=application_window,initialvalue=initial)
    if target is None:
        raise NameError("No target, goodbye.")
    return target

def trim(lis,targ):
    other=[]
    for i in range(len(lis)):
        if lis[i]<=targ:
            other.append(lis[i])
    return other

def trim2(lis,minn,targ):
    other=[]
    for i in range(len(lis)):
        if lis[i]<=targ and lis[i]>=minn:
            other.append(lis[i])
    return other

def trim3(lis,minn,targ):
    other=[]
    for i in range(len(lis)):
        if lis[i][0]<=targ and lis[i][0]>=minn:
            other.append(lis[i])
    return np.array(other)

def endnumber(stri):
    beg=-1
    while stri[beg].isdigit():
        beg-=1
    if beg==-1:
        return 0
    else:
        print("length: "+(stri[beg+1:]))
        return int(stri[beg+1:])

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
    edge=endnumber(objname)
    st=gax.imshow(ston,cmap="cubehelix", vmin=0, vmax=target)
    center=(len(ston[0])/2-0.5,len(ston)/2-0.5)
    #gax.plot([center[0],center[0]],[0,len(ston)],color="red")
    plt.colorbar(st)
    #plt.show()
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
    ## convoluted convolution stuff first
    weights=[]
    xs=[]
    ys=[]
    weighting=ndimage.convolve(ston,np.ones((10,10)),mode='reflect') ##reflect to avoid edge effects
    ybin=[]
    option=1
    negaoption=False
    for y in range(len(ston)):
        for x in range(len(ston[y])):
            ybin.append(ston[y][x])
            weights.append(weighting[y][x])
            xs.append(x)
            ys.append(y)
    fig,ax=plt.subplots()
    if negaoption:
        ## moves up the minimum acceptable cutoff based off the neg values on assumption this comes from noise
        negacc=len(trim(ybin,0))
        if negacc==0:
            minacc=0
        else:
            minacc=np.sort(ybin)[int(1.95*negacc)]
    else:
        minacc=0
    if (minacc>=target):
        option=3
    else:
        trimmed=trim2(ybin,minacc,target)

    if option==0:
        ## progressive cutoff
        big=np.array([ybin,weights,ys,xs]).T
        trimmed=trim3(big,np.min(ybin),0)
        minsize=100
        binsize=int(len(trimmed)/100)
        nn,bins,patches=ax.hist(trimmed[:,0],binsize,alpha=0.5,edgecolor='black',linewidth=1)
        while(np.min(nn)<minsize and binsize>1):
            if(nn[0]<0.2*minsize and nn[1]<0.2*minsize):
                ## suggests we have outliers and should retract the bins
                print("SNIPP")
                trimmed=trim3(trimmed,bins[1],0)
            binsize-=1
            #print(str(binsize)+"  "+str(np.max(nn)))
            nn,bins,patches=ax.hist(trimmed[:,0],binsize,alpha=0.5)
            ax.clear()
        nn,bins,patches=ax.hist(trimmed[:,0],binsize,alpha=0.5,edgecolor='black',linewidth=1)
        ax.clear()
        ## bins has the binedges, nn has number per bin, big.T has the list of [ybin,weight]
        bins=-np.flip(bins)
        nn=np.flip(nn)
        bigsort=big[np.argsort(ybin)]
        temp=[]
        tem=1
        for q in range(len(bigsort)):
            if tem>=len(bins):
                break
            elif bigsort[q][0]<=bins[tem-1]:
                y=int(bigsort[q][2])
                x=int(bigsort[q][3])
                wvt[y][x]=0
                ston[y][x]=0
            elif bigsort[q][0]>bins[tem-1] and bigsort[q][0]<=bins[tem]:
                temp.append(bigsort[q])
            elif bigsort[q][0]>bins[tem]:
                temp2=np.array(temp)
                temp2=temp2[np.argsort(temp2[:,1])]
                i=0
                while(i<nn[tem-1] and i<len(temp2)):
                    y=int(temp2[i][2])
                    x=int(temp2[i][3])
                    wvt[y][x]=0
                    ston[y][x]=0
                    i+=1
                temp=[bigsort[q]]
                tem+=1
        
        ybin2=[]   
        for y in range(len(ston)):
            for x in range(len(ston[y])):
                if ston[y][x]!=0:
                    ybin2.append(ston[y][x])  
        
        nn,bins,patches=ax.hist(ybin,120,alpha=0.4,edgecolor='black',linewidth=1) 
        nn,bins,patches=ax.hist(ybin2,bins,alpha=0.5,edgecolor='black',linewidth=1) 
        plt.show()
    elif option==1:
        minsize=100
        ## equal length bin method:
        binsize=int(len(trimmed)/minsize)
        nn,bins,patches=ax.hist(trimmed,binsize,alpha=0.5)
        ax.clear()
        while(np.min(nn)<minsize and binsize>1):
            binsize-=1
            #print(str(binsize)+"  "+str(np.max(nn)))
            nn,bins,patches=ax.hist(trimmed,binsize,alpha=0.5)
            ax.clear()
        nn,bins,patches=ax.hist(trimmed,binsize,alpha=0.5,edgecolor='black',linewidth=1)
        minb=np.argmin(nn)
        cutoff=0.5*(bins[minb]+bins[minb+1])
    elif option==2:
        minsize0=100
        minsize=minsize0
        ## equal population bin method
        trimmed=np.sort(trimmed)
        repeat=True
        while(repeat):
            ends=[]
            i=0
            while(i<len(trimmed)):
                ends.append(trimmed[i])
                i+=minsize
            ends.append(trimmed[-1])
            ax.clear()
            nn,bins,patches=ax.hist(trimmed,ends,alpha=0.5,edgecolor='black',linewidth=1)
            repeat=np.min(nn)<minsize0 and minsize<len(trimmed)
            minsize=int(minsize+0.2*minsize0)
        intervals=bins[1:]-bins[:-1]
        minb=np.argmax(intervals)
        cutoff=0.5*(bins[minb]+bins[minb+1])
    else:
        ## hard cutoff
        cutoff=target
    
    if option!=0:
        nn,bins,patches=ax.hist(ybin,100,alpha=0.5,edgecolor='black',linewidth=1)
        ax.plot([target,target],[0,2000],linestyle="dashed",color="green")
        ax.plot([cutoff,cutoff],[0,2000],linestyle="dashed",color="black")
        ax.plot([minacc,minacc],[0,2000],linestyle="dashed",color="red")

        for y in range(len(ston)):
            for x in range(len(ston[y])):
                if ston[y][x]<cutoff:
                    wvt[y][x]=0
                    ston[y][x]=0
        #plt.show()
    plt.close()

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

def makesubfolder(sourcedir,target):
    subfolder="target"+str(round(target,5)).strip("0").strip(".")
    try:
        os.mkdir(sourcedir+"/"+subfolder)
    except:
        pass
        print("already exists")
    return subfolder


if __name__ == "__main__":
    wcsx,signal,var,sourcedir,objname=bin_accretion.initialize(enternew=True)

    signal2=np.copy(signal)
    var2=np.copy(var)
    #var2[signal2<=0]=1e10
    signal2[signal2<=0]=0
    
    

    objname=getname("_".join(objname.split("_")[:-1]))
    target=gettarget()
    #print(sourcedir)
    #print(os.path.exists(sourcedir))
    subfolder=makesubfolder(sourcedir,target)

    #saveston(wcsx,signal,var,sourcedir,objname,subfolder="unbinned")
    #objname="J024815-081723"
    ##empty for directly into sourcedir
    
    
    if target<0:
        weighting=False
        target=-target
    else:
        weighting=True
    binlist,diflist=mainfunc(signal2,var2,target,displayWVT=False,epsilon=-10)
    
    wvt,ston=functions.generate_wvt4(binlist,signal,var,np.full(len(binlist),1),False)
    wvt,ston=functions.generate_wvt3(binlist,signal,var,np.full(len(binlist),1),False)
    
    

    vwvt=functions.generate_wvt(binlist,var)
    #saveiteratedfits(target,wcsx,wvt,vwvt,objname,sourcedir,subfolder=subfolder)
    saveblockoutfits(target,ston,wcsx,wvt,vwvt,objname,sourcedir,subfolder=subfolder)
    saveston(wcsx,ston,sourcedir,objname,subfolder=subfolder)
    assign=functions.assign(binlist,target,ston)
    saveassign(wcsx,assign,sourcedir,objname,subfolder=subfolder)