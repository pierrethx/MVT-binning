import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys,os
sys.path.append(("/".join(os.path.realpath(__file__).split("/")[:-2])))
import bin_accretion, functions
from scipy.ndimage.measurements import center_of_mass

def trim(lis,targ):
    other=[]
    for i in range(len(lis)):
        if lis[i]<=targ:
            other.append(lis[i])
    return other

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

sublist=["target5"]
tlist=[5]
markers=["s"]

tag="CVT"

wcsxlist,siglist,varlist,sourcelist,objlist=bin_accretion.minitialize()
namelist=["_".join(obj.split("_")[:-1]) if "_sig.fits" in obj.lower() else ".".join(obj.split(".")[:-1]) for obj in objlist ]


for i in range(len(sourcelist)):
    fig,ax=plt.subplots()
    fig2,ax2=plt.subplots()
    unrad=[]
    unsn=[]

    ston=siglist[i]/np.sqrt(varlist[i])

    center=center_of_mass(np.ma.masked_where(np.isnan(ston),ston)) #as usual is (y,x)
    scale=0.5*np.sqrt(len(ston)**2+len(ston[0])**2)
    mask=np.copy(ston)
    for yy in range(len(mask)):
        for xx in range(len(mask[yy])):
            if np.sqrt((yy-center[0])**2+(xx-center[1])**2)>scale/3:
                mask[yy][xx]=0
    center=center_of_mass(mask)

    for y in range(len(ston)):
        for x in range(len(ston[0])):
            if not np.isnan(ston[y][x]):
                unrad.append(np.sqrt((y-center[0])**2+(x-center[1])**2))
                unsn.append(ston[y][x])

    unrad=np.array(unrad)
    unsn=np.array(unsn)
    ax.set_xlabel("Radius from center")
    ax.set_ylabel("SNR")
    ax.set_title(namelist[i])
    ax2.set_title(namelist[i])
    fig3,ax3=plt.subplots()
    ax3.set_title(namelist[i])
    minx,maxx=1000,0
    intt=-1
    for sub in sublist:
        intt+=1
        '''
        with fits.open(sourcelist[i]+"/"+sub+"/"+"block_"+namelist[i]+"_wit_sig.fits") as hdul:
            bsig=(np.flipud(hdul[0].data))
        with fits.open(sourcelist[i]+"/"+sub+"/"+"block_"+namelist[i]+"_wit_var.fits") as hdul:
            bvar=(np.flipud(hdul[0].data))
        '''
        with fits.open(sourcelist[i]+"/"+sub+"/"+"z_"+namelist[i]+"_"+tag+"_assigned.fits") as hdul:
            ass=(np.flipud(hdul[0].data))
        with fits.open(sourcelist[i]+"/"+sub+"/"+"zston_"+namelist[i]+"_"+tag+".fits") as hdul:
            bston=(np.flipud(hdul[0].data))

        
        brad=[]
        bsn=[]

        cents=[]
        bareas=[]
        bravs=[]
        binlist,bbl=functions.reverseassign(ass) 
        for bin in binlist:
            c=functions.geometric_center(bin)
            cents.append(c)
            bareas.append(len(bin))
            bravs.append(np.average([np.sqrt((c[0]-p[0])**2+(c[1]-p[1])**2) for p in bin]))
        bareas=np.array(bareas[1:])
        bravs=bravs[1:]
        bronds=[bravs[i]*np.sqrt(np.pi/bareas[i]) for i in range(len(bareas))]
        bronds=np.array(bronds)
        cents=cents[1:]
        
        for cen in cents:
            print(cen,cen[0],cen[1])
            brad.append(np.sqrt((cen[0]-center[0])**2+(cen[1]-center[1])**2))
            bsn.append(bston[int(cen[0]+0.5)][int(cen[1]+0.5)])
        #ax.plot(rads,np.full(len(rads),0),linestyle="dashed",color="red")
        bsn=np.array(bsn)
        brad=np.array(brad)
        ax.plot(brad[bareas>1],bsn[bareas>1],linewidth=0,marker=".",label=sub,zorder=10)
        #ax.plot(brad[brad>60],bsn[brad>60],linewidth=0,marker=".",label=sub,zorder=10)
        
        
        if (len(brad[bareas>1])==0):
            minx=0
            maxx=np.max(brad)
        else:
            if(np.min(brad[bareas>1])<minx):
                minx=np.min(brad[bareas>1])
            if(np.max(brad[bareas>1])>maxx):
                maxx=np.max(brad[bareas>1])
        ax2.plot(brad[bareas>1],bronds[bareas>1],linewidth=0,marker=".",label=sub)
        ax3.plot(brad[bareas>1],bareas[bareas>1],linewidth=0,marker=".",label=sub)
    ax2.plot(np.linspace(minx*0.9,maxx*1.1,10),0*np.linspace(minx*0.9,maxx*1.1,10)+0.667,linewidth=2,color="k",linestyle="dashed",label="ideal limit")
    #ax.set_xlim(minx*0.9,maxx*1.1)
    #ax.plot(unrad[unrad>minx*0.9],unsn[unrad>minx*0.9],linewidth=0,marker=".",label="unbinned",zorder=0)
    ax.plot(unrad[unrad>minx],unsn[unrad>minx],linewidth=0,marker=".",label="unbinned",zorder=0)
    for intt in range(len(sublist)):
        ax.plot(np.linspace(minx,maxx*1.1,10),0*np.linspace(minx,maxx*1.1,10)+tlist[intt],linewidth=1,color="k",linestyle="dashed",zorder=5)
        #ax.plot(np.linspace(minx*0.9,maxx*1.1,10),0*np.linspace(minx*0.9,maxx*1.1,10)+tlist[intt],linewidth=1,color="k",linestyle="dashed",zorder=5)
    ax.set_ylim(0,np.nanmax(unsn[unrad>minx*0.9]))
    ax2.set_xlabel("Radius from center")
    ax2.set_ylabel("Roundness Parameter")
    fig.legend()
    fig2.legend()
    
    #plt.show()
    fig.savefig(sourcelist[i]+"/"+namelist[i]+"_"+tag+"_SNR.png")
    fig2.savefig(sourcelist[i]+"/"+namelist[i]+"_"+tag+"_round.png")

    
    #ax3.plot([2,np.max(brad)],[30,30],"k:",linewidth=2,label="input minimum size")
    ax3.set_xlabel("Radius from center")
    ax3.set_ylabel("Bin Size (pixels)")
    fig3.legend()
    fig3.savefig(sourcelist[i]+"/"+namelist[i]+"_"+tag+"_areas.png")
    
