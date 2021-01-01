import scipy.spatial as sp

import numpy as np
import matplotlib.pyplot as plt
import tkinter
from tkinter.filedialog import askopenfilename
from astropy.io import fits
import bin_accretion
from astropy import wcs
from radial_profile import getcenter
def trim(lis,targ):
    other=[]
    for i in range(len(lis)):
        if lis[i]<=targ:
            other.append(lis[i])
    return other

def trim2(lis,lis2):
    other=[]
    for i in range(len(lis)):
        if lis2[i]>1:
            other.append(lis[i])
    return other

'''
# make sig noise -> ston
root=tkinter.Tk()
root.withdraw()
signalname= askopenfilename(message="Select unbinned ston")
root.update()
root.destroy()
'''

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

wcsx,slist,vlist,sdlist,olist=bin_accretion.minitialize()

for r in range(len(slist)):
    signal=slist[r]
    var=vlist[r]
    sourcedir=sdlist[r]
    objname=olist[r]

    ston=signal/np.sqrt(var)

    fig,ax=plt.subplots(2,2)
    xun=[]
    yun=[]
    '''
    with fits.open(signalname) as hdul:
        ston=np.flipud(hdul[0].data)
        wcsx=wcs.WCS(hdul[0].header)
    '''

    center=getcenter(np.ma.masked_where(np.isnan(ston),ston)) #as usual is (y,x)


    for y in range(len(ston)):
        for x in range(len(ston[0])):
            if not np.isnan(ston[y][x]):
                xun.append(np.sqrt((y-center[0])**2+(x-center[1])**2))
                yun.append(ston[y][x])
    ax[0][0].plot([0,100],[0,0],color="k")
    ax[0][0].set_ylim(-5,50)
    ax[0][0].set_xlim(0,100)
    ax[0][0].plot(xun,yun,linewidth=0,marker=".",label="unbinned")
    ax[1][0].set_xlabel("Radius from center")
    ax[1][0].set_ylabel("S/N from center")


    target=[3,5,10]
    xes=[]
    yes=[]
    aes=[]

    col=["tab:orange","tab:green","tab:red"]

    yes2=[]

    '''
    source="/".join(signalname.split("/")[:-2])
    objname=signalname.split("/")[-1]
    stonfits=".".join(objname.split(".")[:-1])+"_wit.fits"
    '''
    source="/".join(sourcedir.split("/")[:-1])
    stonfits="zston_"+"_".join(objname.split("_")[:-1])+".fits"
    assfits="z_"+"_".join(objname.split("_")[:-1])+"_assigned.fits"


    for t in range(len(target)):
        xes.append([])
        yes.append([])
        yes2.append([])
        aes.append([])
        signalname=source+"/target"+str(target[t])+"/"+stonfits
        print(signalname)
        with fits.open(signalname) as hdul:
            ston=np.flipud(hdul[0].data)
            wcsx=wcs.WCS(hdul[0].header)
        signalname=source+"/target"+str(target[t])+"/"+assfits
        print(signalname)
        with fits.open(signalname) as hdul:
            ass=np.flipud(hdul[0].data.astype(int))
            wcsx=wcs.WCS(hdul[0].header)

        bins=int(np.nanmax(ass))
        bcenters=[[0,0] for i in range(bins)]

        binlist=[0 for i in range(bins)]
        manual=[]
        mlist=[]
        mcenters=[]
        for y in range(len(ass)):
            for x in range(len(ass[y])):
                yes2[t].append(ston[y][x])
                if ass[y][x]>0:
                    binlist[ass[y][x]-1]+=1
                    bcenters[ass[y][x]-1][0]+=y
                    bcenters[ass[y][x]-1][1]+=x
                else:
                    manual.append((y,x))
            #print(str(y/len(ass))+" through stage 1")
        g=len(manual)
        for m in manual:
            g-=1
            print(g)
            found=False
            for l in mlist:
                for p in l:
                    if ((p[0]-m[0])**2<=1 or (p[1]-m[1])**2<=1) and ston[p[0]][p[1]]==ston[m[0]][m[1]]:
                        found=True
                        break
                if found:
                    l.append(m)
                    break
            if not found:
                mlist.append([m])
        for l in mlist:
            xx=0
            yy=0
            for p in l:
                xx+=p[1]
                yy+=p[0]
            binlist.append(len(l))
            bcenters.append([yy,xx])
        for l in range(len(binlist)):
            if binlist[l]>0:
                bcenters[l][0]/=binlist[l]
                bcenters[l][1]/=binlist[l]

        center=getcenter(ston) #as usual is (y,x)

        for point in range(len(bcenters)):
            x=int(bcenters[point][1])
            y=int(bcenters[point][0])
            if not np.isnan(ston[y][x]):
                xes[t].append(np.sqrt((y-center[0])**2+(x-center[1])**2))
                yes[t].append(ston[y][x])
        m=(t+1)%2
        n=int((t+1)/2)
        ax[n][m].plot([0,100],[0,0],color="k")
        ax[n][m].plot([0,100],[target[t],target[t]],linestyle="dashed",color="k")
        ax[n][m].set_ylim(-5,50)
        ax[n][m].set_xlim(0,100)
        ax[n][m].plot(xes[t],yes[t],linewidth=0,marker=".",color=col[t],label="target"+str(target[t]))

    fig.legend()
    #plt.savefig("/".join(sourcedir.split("/")[:-3])+"/images/a_"+"_".join(objname.split("_")[:-1])+"_".join(sourcedir.split("/")[-3:-1])+".png")
    plt.close()
    '''
    plt.show()
    '''
    fig2,ax2=plt.subplots(2,2)

    binsz=100



    a=ax2[0][0].hist(trim(yun,50),binsz,label="unbinned")
    ax2[0][0].plot([0,0],[0,np.max(a[0])],color="k")
    ax2[1][0].set_xlabel("S-to-N")
    ax2[1][0].set_ylabel("frequency")
    for t in (range(len(target))):
        m=(t+1)%2
        n=int((t+1)/2)
        
        a = ax2[n][m].hist(trim(yes2[t],50),binsz,color=col[t],label="target"+str(target[t]))
        ax2[n][m].plot([0,0],[0,np.max(a[0])],color="k")
        ax2[n][m].plot([target[t],target[t]],[0,np.max(a[0])],linestyle="dashed",color="k")
    fig2.legend()
    #plt.savefig("/".join(sourcedir.split("/")[:-3])+"/images/b_"+"_".join(objname.split("_")[:-1])+"_".join(sourcedir.split("/")[-3:-1])+".png")
    #plt.close()
    
    plt.show()
    
print("Bye Bye!")