import numpy as np
import matplotlib.pyplot as plt 
import tkinter as tk
from tkinter.filedialog import askopenfilename
from astropy.io import fits
from astropy import wcs
import functions,bin_accretion,wvt_iteration
import scipy.spatial as sp
from scipy import ndimage
import time, os
from collections import Counter


SMALL_SIZE=14
MEDIUM_SIZE=16
BIGGER_SIZE=20

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

def assign(binlist,target,ston):
    ## generates an assignment file from a binlist. This makes it easier to reconstruct a binning
    ## without saving the binlist as a separate file or anything.
    ## every bin gets a unique integer value
    ## here we have inadequate bins with negative integer value, and nothing as 0. this makes it easier to visualize how the 
    ## binning went when we have large amounts of invalid bins that we want to look closer at
    assign=np.zeros_like(ston)
    binlist2=binlist.copy()
    np.random.shuffle(binlist2)
    g=0
    for i in range(len(binlist2)):
        k0=binlist2[i][0]
        if ston[k0[0]][k0[1]]==0:
            for k in binlist2[i]:
                assign[k[0]][k[1]]=0
        else:
            g=g+1
            for k in binlist2[i]:
                assign[k[0]][k[1]]=g      
    return assign

def saveunblockedfits(wcsx,wvt,vwvt,objname,sourcedir,mode,subfolder="unbinned",check=1):
    ## for visualizing results

    hdu = fits.PrimaryHDU(wvt,header=wcsx)
    hdul = fits.HDUList([hdu])
    hdul.writeto(sourcedir+"/"+subfolder+"/"+objname+"_"+mode+"_sig.fits",overwrite=True,checksum=True)
   
    hdu2 = fits.PrimaryHDU(vwvt,header=wcsx)
    hdul2 = fits.HDUList([hdu2])

    hdul2.writeto(sourcedir+"/"+subfolder+"/"+objname+"_"+mode+"_var.fits",overwrite=True,checksum=True)
 
def savestonimg(target,ston,sourcedir,objname,subfolder,mode):
    gig,gax=plt.subplots()
    st=gax.imshow(ston,cmap="cubehelix",origin="lower", vmin=0, vmax=2*target)
    plt.colorbar(st)
    #plt.show()
    plt.savefig(sourcedir+"/"+subfolder+"/block_"+objname+"_"+mode+"_ston.png")
    plt.close()
    
def saveblockoutfits(target,ston,wcsx,wvt,vwvt,objname,sourcedir,mode,subfolder="unbinned",check=1):
    blockout(target,wvt,ston,sourcedir+"/"+subfolder+"/block_"+objname+"_"+mode+"_hist.png",check)

    ## for visualizing results
    if check>0:
        gig,gax=plt.subplots()
        st=gax.imshow(ston,cmap="cubehelix",origin="lower", vmin=0, vmax=2*target)
        plt.colorbar(st)
        #plt.show()
        plt.savefig(sourcedir+"/"+subfolder+"/block_"+objname+"_"+mode+"_ston.png")
        plt.close()

    hdu = fits.PrimaryHDU(wvt,header=wcsx)
    hdul = fits.HDUList([hdu])
    #manipulate(hdul)
    hdul.writeto(sourcedir+"/"+subfolder+"/block_"+objname+"_"+mode+"_sig.fits",overwrite=True,checksum=True)

    hdu2 = fits.PrimaryHDU(vwvt,header=wcsx)
    hdul2 = fits.HDUList([hdu2])
    #manipulate(hdul2)
    hdul2.writeto(sourcedir+"/"+subfolder+"/block_"+objname+"_"+mode+"_var.fits",overwrite=True,checksum=True)

def saveston(wcsx,ston,sourcedir,objname,mode,subfolder="unbinned"):
    hdu3 = fits.PrimaryHDU(ston,header=wcsx)
    hdul3 = fits.HDUList([hdu3])
    #manipulate(hdul3)
    hdul3.writeto(sourcedir+"/"+subfolder+"/zston_"+objname+"_"+mode+".fits",overwrite=True,checksum=True)

def saveassign(wcsx,assign,sourcedir,objname,mode,subfolder="unbinned"):
    hdu3 = fits.PrimaryHDU(assign,header=wcsx)
    hdul3 = fits.HDUList([hdu3])
    #manipulate(hdul3)
    hdul3.writeto(sourcedir+"/"+subfolder+"/z_"+objname+"_"+mode+"_assigned.fits",overwrite=True,checksum=True)

def convergencelist(epz,diflistlist,targlist,sourcedir,objname,mode,subfolder="_"):
    ## generates a chart to see how a function converges over many iterations. 
    # This is to check condition 4 of the function converging over
    if subfolder=="_":
        if len(targlist)==1:
            fpath=sourcedir+"/"+"target"+str(round(targlist[0],5)).strip("0").strip(".")+"_"+mode+"/y_"+objname+"_"+mode+"_convergence.png" 
        else:
            fpath=sourcedir+"/unbinned/y_"+objname+"_"+mode+"_convergence.png"
    else:
        fpath=sourcedir+"/"+subfolder+"/y_"+objname+"_"+mode+"_convergence.png"

    fig,ax=plt.subplots()
    ax.set_xlabel("iteration number")
    ax.set_ylabel("normalized difference")

    for i in range(len(diflistlist)):
        diflist=np.where(np.isnan(diflistlist[i]),np.zeros_like(diflistlist[i]),diflistlist[i])[1:]
        
        xes=range(1,len(diflist)+1)
        if(epz>0):
            ax.plot([1,len(diflist)],[epz,epz],linestyle="dashed")
            
        ax.plot(xes,diflist,marker="o",label=targlist[i])
    ax.set_yscale('log')
    fig.legend(title="Target SNR")
    plt.tight_layout()
    
    plt.savefig(fpath)
    plt.close()

def blockout(target,wvt,ston,fpath,check=1):
    ## convoluted convolution stuff first
    weights=[]
    xs=[]
    ys=[]
    weighting=ndimage.convolve(ston,np.ones((10,10)),mode='reflect') ##reflect to avoid edge effects
    ybin=[]
    option=2
    negaoption=True
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
        binnum=int(len(trimmed)/100)
        nn,bins,patches=ax.hist(trimmed[:,0],binnum,alpha=0.5,edgecolor='black',linewidth=1)
        while(np.min(nn)<minsize and binnum>1):
            if(nn[0]<0.2*minsize and nn[1]<0.2*minsize):
                ## suggests we have outliers and should retract the bins
                print("SNIPP")
                trimmed=trim3(trimmed,bins[1],0)
            binnum-=1
            #print(str(binnum)+"  "+str(np.max(nn)))
            nn,bins,patches=ax.hist(trimmed[:,0],binnum,alpha=0.5)
            ax.clear()
        nn,bins,patches=ax.hist(trimmed[:,0],binnum,alpha=0.5,edgecolor='black',linewidth=1)
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
        ## minimum histogram method:
        #print(len(trimmed))
        if len(trimmed)==0:
            cutoff=0
        else:
            nn,bins,binnum=histify(target,trimmed,minsize,minacc,ax)
            minb=np.argmin(nn)
            cutoff=0.5*(bins[minb]+bins[minb+1])
            blockoutchex(target,check,ax,minacc,cutoff,bins,binnum,ybin,fpath)
    elif option==2:
        ## maximum histogram method:
        minsize=100
        if len(trimmed)==0:
            cutoff=0
        else:
            nn,bins,binnum=histify(target,trimmed,minsize,minacc,ax)
            maxb=np.argmax(nn)
            current=maxb-1
            minb=maxb
            if(maxb==0):
                ## we are on a noise peak, roll forward
                while current<binnum and nn[current]<=nn[minb]:
                    minb=current
                    current+=1
            else:
                ## hopefully on a signal peak, roll backwards
                while current>=0 and nn[current]<=nn[minb]:
                    minb=current
                    current-=1
                ## and then find minimum
                while current>=0:
                    if nn[current]<=nn[minb]:
                        minb=current
                    current-=1
            cutoff=0.5*(bins[minb]+bins[minb+1])
            blockoutchex(target,check,ax,minacc,cutoff,bins,binnum,ybin,fpath)
    else:
        ## hard cutoff
        cutoff=target
    
    if option!=0:
        """
        nn,bins,patches=ax.hist(ybin,100,alpha=0.5,edgecolor='black',linewidth=1)
        ax.plot([target,target],[0,2000],linestyle="dashed",color="green")
        ax.plot([cutoff,cutoff],[0,2000],linestyle="dashed",color="black")
        ax.plot([minacc,minacc],[0,2000],linestyle="dashed",color="red")
        """
        for y in range(len(ston)):
            for x in range(len(ston[y])):
                if ston[y][x]<cutoff:
                    wvt[y][x]=0
                    ston[y][x]=0
        #plt.show()
    plt.close()

    return wvt,ston

def blockoutchex(target,check,ax,minacc,cutoff,bins,binnum,ybin,fpath):
    if check>0:
        ax.clear()
        step=(target-minacc)/binnum
        bins2=[i for i in bins]
        maxxx=target*2
        minnn=-target
        while(bins2[0]>minnn):
            bins2.insert(0,bins2[0]-step)
        while(bins2[-1]<maxxx):
            bins2.append(bins2[-1]+step)
        nn,bins,patches=ax.hist(trim2(ybin,minnn,maxxx),5*len(bins2),edgecolor='black',linewidth=1,histtype='step',fill=True,facecolor='black')
        nn,bins,patches=ax.hist(trim2(ybin,minnn,maxxx),bins2,edgecolor='red',linewidth=2,histtype='step')
        ax.vlines(cutoff,0,np.max(nn)*1.05,"black",linestyle='dashed',linewidth=3,label="minimum")
        #ax.vlines([minacc,target],0,np.max(nn),"blue",label="cut range")
        ax.axvspan(minnn,minacc,facecolor="blue",alpha=0.5)
        ax.axvspan(target,maxxx,facecolor="blue",alpha=0.5)
        ax.set_xlim(minnn,maxxx)
        ax.set_xlabel("SNR")
        ax.set_ylabel("frequency (pixels)")
        #plt.savefig(fpath[:-4]+"0.png")
        if check>1:
            plt.gcf().show()
            cup=input("New cutoff? ")
            try:
                cutoff=float(cup)
                ax.vlines(cutoff,0,np.max(nn)*1.05,"black",linestyle='dotted',linewidth=3,label="user-placed")
                print("New cutoff selected at "+cup)
            except:
                print("No new cutoff selected")
        #plt.tight_layout()
        plt.legend()
        plt.savefig(fpath,bbox_inches="tight")

def histify(target,trimmed,minsize,minacc,ax):
    binnum=int(len(trimmed)/minsize)
    while binnum==0:
        minsize=int(minsize/2)
        binnum=int(len(trimmed)/minsize)
        print(binnum)
    nn,bins,patches=ax.hist(trimmed,binnum,alpha=0.5,histtype='step')
    ax.clear()
    nindices=Counter(np.digitize(trimmed,np.linspace(minacc,target,binnum+1)))

    incompleteness=True
    while (nindices.most_common()[-1][1]<minsize or incompleteness) and binnum>1:
        binnum-=1
        print(binnum)
        nindices=Counter(np.digitize(trimmed,np.linspace(minacc,target,binnum+1)))
        incompleteness=False
        for n in range(1,binnum+1):
            if(nindices[n]==0):
                incompleteness=True
                break
    nn,bins,patches=ax.hist(trimmed,np.linspace(minacc,target,binnum+1),alpha=0.5,edgecolor='black',linewidth=1,histtype='step')
    return nn,bins,binnum

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

if __name__ == "__main__":

    ## negative epsilon means iterate that number of times
    ## positive epsilon gives us a tolerance our iterations must get below
    eps=-30

    ## minsize depends on smallest viable cell size for resolution, which depends on what instrument you are using
    ## minsize too large on WVT might break it. Try to avoid that 
    minsize=30

    ## this toggles between "WVT" and "CVT". See Diehl and Statler 2006 for WVT and Cappellari and Copin 2003 for CVT, 
    ## though is modified by us
    ## as well as "WVT2s" and "VT" which are minimum size special methods explained in senior thesis
    ## "WVT2s" bin accretes, saves interior bins, then iterates only outer, low SNR bins
    ## VT does WVT but does not use scalelengths, or CVT with geometric center rather than weighted centroid
    modetypes=[]

    modetypes.append("VT")
    modetypes.append("WVT2s")
    modetypes.append("WVT")
    modetypes.append("CVT")
    
    ## here for easy testing we allow you to select multiple options but 
    ## unless you were actually using multiple options you would just comment the rest out


    targhold="10,5,3"
    wcsxlist,siglist,varlist,sourcelist,objlist=bin_accretion.minitialize()
    if len(siglist)==0:
        print("No files loaded!")
    else:
        print("Files loaded!")
        targlist=bin_accretion.getmulttarget(targhold)
        for i in range(len(sourcelist)):
            for modetype in modetypes:
                diflistlist=[]
                wcsx=wcsxlist[i]
                signal=siglist[i]
                var=varlist[i]
                if "_sig.fits" in objlist[i].lower():
                    objname="_".join(objlist[i].split("_")[:-1])
                else:
                    objname=".".join(objlist[i].split(".")[:-1])
                sourcedir=sourcelist[i]
                print(sourcedir)
                for m in range(len(targlist)):

                    target=targlist[m]
                    
                    subfolder=bin_accretion.makesubfolder(sourcedir,target,modetype)

                    ## this is us applying the data mask. We block out negative values and then run the binning algorithm on it
                    signal2=np.copy(signal)
                    var2=np.copy(var)
                    var2[signal2<=0]=1e10
                    signal2[signal2<=0]=0

                    ## generally the incidence of masked out pixels is
                    test=np.zeros_like(signal)
                    test[signal<=0]=1
                    incidence=np.average(test) # this is the "percentage of masked pixels"
                    
                    ## this is the most important function
                    ## first we do bin accretion. then iteration. the wvt and vwvt is not necessary here but im scared to remove it
                    binlist,init_generators,init_scalelengths=bin_accretion.cc_accretion(signal2,var2,target,minsize=minsize,mode=modetype,display=False)
                    binlist,diflist=wvt_iteration.iteration_moderator(target,signal2,var2,init_generators,init_scalelengths,eps,minsize,incidence,mode=modetype,display=False)

                    ## now to generate the binned signal and variance
                    wvt,ston=functions.generate_wvt3(binlist,signal,var,np.full(len(binlist),1))
                    vwvt=functions.generate_wvt(binlist,var)

                    saveblockoutfits(target,ston,wcsx,wvt,vwvt,objname,sourcedir,modetype,subfolder,check=1)
                    #saveunblockedfits(wcsx,wvt,vwvt,objname,sourcedir,modetype,subfolder)
                    saveston(wcsx,ston,sourcedir,objname,modetype,subfolder)
                    #savestonimg(target,ston,sourcedir,objname,subfolder,modetype+"0")
                    assig=assign(binlist,target,ston)
                    saveassign(wcsx,assig,sourcedir,objname,modetype,subfolder)
                    diflistlist.append(diflist)
                if(len(diflistlist[0])>0):
                    convergencelist(eps,diflistlist,targlist,sourcelist[i],objname,modetype)
    print("Bye bye")