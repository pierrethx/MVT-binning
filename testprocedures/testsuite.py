import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys,os
sys.path.append(("/".join(os.path.realpath(__file__).split("/")[:-2])))
import bin_accretion, functions
from scipy.ndimage.measurements import center_of_mass
from scipy import stats

## There are some specific things about graph boundaries that are specific to the dataset 
## being processed here, but in general, this code displays how to
## track the properties mentioned in the paper. just change the graph/data boundaries


def trim(lis,targ):
    other=[]
    for i in range(len(lis)):
        if lis[i]<=targ:
            other.append(lis[i])
    return other

def trim2(lis1,lis2,mi,ma):
    other1=[]
    other2=[]
    for i in range(len(lis1)):
        if lis1[i]<=ma and lis1[i]>=mi:
            other1.append(lis1[i])
            other2.append(lis2[i])
    return np.array(other1),np.array(other2)

def makesubfolder(sourcedir,subfolder):
    for subfol in subfolder:
        try:
            os.mkdir(sourcedir+"/"+subfol)
        except:
            pass
            #print("already exists")
    return subfolder

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

wcsxlist,siglist,varlist,sourcelist,objlist=bin_accretion.minitialize()
namelist=["_".join(obj.split("_")[:-1]) if "_sig.fits" in obj.lower() else ".".join(obj.split(".")[:-1]) for obj in objlist ]

["CVT0","WVT0","CVT","WVT","VT","WVT2s"]

print("\n\n\n\n\n")
tagz=["CVT","VT","WVT","WVT2s"]

for i in range(len(sourcelist)):
    
    tlist=[10,5,3]
    markers=["s"]
    
    minsize=30
    minx,maxx=40,80
    '''
    minsize=1
    minx,maxx=10,80
    '''
    

    def theoretical_area(radius,target):
        I=1
        rc=10
        beta=0.67
        gain=1000
        bg=0.01
        signal=I*(1+(radius/rc)**2)**(0.5-3*beta)
        variance=(signal+bg)/gain
        sNR=signal/np.sqrt(variance)
        out=(target/sNR)**2
        out[out<minsize]=minsize
        return out

    styles=[]
    for tag in tagz:
        sublist=["target10_"+tag,"target5_"+tag,"target3_"+tag]
        makesubfolder(sourcelist[i],["SNR","round","areas"])

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
        
        fig3,ax3=plt.subplots()
        if tag=="CVT0":
            namez="CVT Mode"
        elif tag=="WVT0":
            namez="WVT Mode"
        else:
            namez= "Masked "+tag+" Mode"
        ax.set_title(namez)
        ax2.set_title(namez)
        ax3.set_title(namez)
        
        intt=-1
        maxnn=0
        bbronds=[]
        bbweights=[]

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
            
            for s in range(1,len(binlist)):
                cen=binlist[s][0]
                #print(cen,cen[0],cen[1])
                brad.append(np.sqrt((cen[0]-center[0])**2+(cen[1]-center[1])**2))
                stonk=bston[int(cen[0]+0.5)][int(cen[1]+0.5)]
                bsn.append(stonk)
            #ax.plot(rads,np.full(len(rads),0),linestyle="dashed",color="red")
            bsn=np.array(bsn)
            brad=np.array(brad)
            if (np.max(brad)>maxx):
                maxx=np.max(brad)
            brad_=brad[bareas>1]
            brad_s=np.sort(brad_)

            minrad=brad[brad>minx]
            minsn=bsn[brad>minx]
            minronds=bronds[brad>minx]
            ax.plot(minrad,minsn,linewidth=0,marker=".",label=sub,zorder=10-intt)
            
            #ax.plot(brad[brad>60],bsn[brad>60],linewidth=0,marker=".",label=sub,zorder=10)

            #ax2.hist(minrad,minronds,linewidth=0,marker=".",label=sub)
            bbronds.append(bronds[bareas>1])
            if intt==1:
                styles.append(bronds[bareas>1])
            bbweights.append(np.ones_like(bronds[bareas>1]) / len(bronds[bareas>1]))
            
            ax3.plot(brad_,bareas[bareas>1],linewidth=0,marker=".",label=sub)
            theor=theoretical_area(np.linspace(0,75,200),tlist[intt])
            ax3.plot(np.linspace(0,78,200),theor,linewidth=2,color="r")
        nn,binss,pathces=ax2.hist(bbronds,bins=np.linspace(0.6,0.8,13),histtype='bar',label=sub,weights=bbweights)
        if np.max(nn)>maxnn:
            maxnn=np.max(nn)
        ax2.plot([0.667,0.667],[0,85],linewidth=2,color="k",linestyle="dashed")
        ax2.set_xlim(0.6,0.8)
        ax2.set_ylim(0,0.55)

        ax.set_xlim(minx,85)
        #ax3.set_yscale("log")
        ax3.set_xlim(10,80)
        ax3.set_ylim(0,250)
        #ax3.set_ylim(-8,175)

        
        for k in range(len(bbronds)):
            for j in range(k):
                sss=stats.ks_2samp(bbronds[k],bbronds[j])
                if(sss[1]>0.01):
                    print("KS Test for "+sub+"_"+sourcelist[i])
                    print(tlist[k],tlist[j],sss)

        #ax.set_xlim(minx*0.9,maxx*1.1)
        #ax.plot(unrad[unrad>minx*0.9],unsn[unrad>minx*0.9],linewidth=0,marker=".",label="unbinned",zorder=0)
        
        munrad,munsn=trim2(unrad,unsn,minx,85)
        
        ax.plot(munrad,munsn,linewidth=0,marker=".",label="unbinned",zorder=0)
        for intt in range(len(sublist)):
            ax.plot(np.linspace(minx,85,10),0*np.linspace(minx,85,10)+tlist[intt],linewidth=1,color="k",linestyle="dashed",zorder=5)
            #ax.plot(np.linspace(minx*0.9,maxx*1.1,10),0*np.linspace(minx*0.9,maxx*1.1,10)+tlist[intt],linewidth=1,color="k",linestyle="dashed",zorder=5)
        ax2.set_ylabel("normalized frequency")
        ax2.set_xlabel("Roundness Parameter")
        #fig.legend(loc='upper right')
        #fig2.legend(loc='upper right')
        

        
        #plt.show()
        fig.savefig(sourcelist[i]+"/SNR/"+namelist[i]+"_"+tag+"_SNR.png")
        fig2.savefig(sourcelist[i]+"/round/"+namelist[i]+"_"+tag+"_round.png")

        
        #ax3.plot([2,np.max(brad)],[30,30],"k:",linewidth=2,label="input minimum size")
        ax3.plot(brad_s,minsize+0*brad_s,"k:",linewidth=3)
        
        ax3.set_xlabel("Radius from center")
        ax3.set_ylabel("Bin Size (pixels)")
        #fig3.legend(loc='upper left')
        fig3.savefig(sourcelist[i]+"/areas/"+namelist[i]+"_"+tag+"_areas.png")
        plt.close(fig)
        plt.close(fig2)
        plt.close(fig3)
        
    print("Different binnings")
    for h in range(len(styles)):
        for g in range(h):
            sss=stats.ks_2samp(styles[h],styles[g])
            if(sss[1]>0.01):
                print("KS Test for target 5 "+sourcelist[i])
                print(tagz[h],tagz[g],sss)