import numpy as np
import matplotlib.pyplot as plt 
import tkinter
from tkinter.filedialog import askopenfilename
from astropy.io import fits
import functions,bin_accretion,wvt_iteration,main
import scipy.spatial as sp
import time

def switch(version,target,binlist,signal,var):
    if version==4:
        binlist,geocarray,scalearray=functions.calculate_scales4(target,binlist,signal,var)
    elif version==2:
        binlist,geocarray,scalearray=functions.calculate_scales2(target,binlist,signal,var)
    elif version==3:
        binlist,geocarray,scalearray=functions.calculate_scales3(target,binlist,signal,var)
    else:
        binlist,geocarray,scalearray=functions.calculate_scales(target,binlist,signal,var)
    return binlist,geocarray,scalearray

def init(signdensity,ston):
    if signdensity:
        return ston*np.abs(ston)
    else:
        return ston*ston

def switch_cc_accretion(version,signdensity,signal,var,target):

    ## this is the signal-to-noise for each pixel
    ston=signal/np.sqrt(np.abs(var))

    
    assign=np.full_like(ston,-1)
            
    density = init(signdensity,ston)

    cellsleft=np.count_nonzero(assign == -1)
    print(cellsleft)

    supercentroid=np.unravel_index(ston.argmax(),ston.shape)
    supermass=density[supercentroid[0]][supercentroid[1]]
    
    viable=[]
    bin_accretion.validateappend(viable,(supercentroid[0],supercentroid[1]),assign)
    binlist=[]
    bcentroids=[]
    rcentroids=[]
    rebinlist=[]
    
    ## in case this borks np.count_nonzero(array == value) to return -1 in assign
    while len(viable)>0:
        centroid=functions.closest_point(supercentroid,viable,density)
        current=[]
        current.append(centroid)
        viable.remove(centroid)
        assign[centroid[0]][centroid[1]]=0
        viablecell=[]
        bin_accretion.validateappend(viablecell,(centroid[0]+1,centroid[1]),assign)
        bin_accretion.validateappend(viablecell,(centroid[0]-1,centroid[1]),assign)
        bin_accretion.validateappend(viablecell,(centroid[0],centroid[1]+1),assign)
        bin_accretion.validateappend(viablecell,(centroid[0],centroid[1]-1),assign)
        binmass=density[centroid[0]][centroid[1]]
        accrete=binmass<target**2 and len(viablecell)>0
        while accrete:
            nextpoint=functions.closest_point(centroid,viablecell,density)

            Rmax=0.3
            newbin=current+[nextpoint]
            nmass=binmass+density[nextpoint[0]][nextpoint[1]]
            ## replacing weighted centroids with geometric centers as suggested by Diehl 
            ## in Cappellari's implementation to address negative data
            ncentroid=functions.geometric_center(newbin)
            rmax=sp.distance.cdist([ncentroid],newbin).max()
            R=rmax*np.sqrt(np.pi/len(newbin))-1
            if R<=Rmax and np.abs(binmass-target**2)>np.abs(nmass-target**2):
                current.append(nextpoint)
                viablecell.remove(nextpoint)
                if viable.count(nextpoint)>0:
                    viable.remove(nextpoint)
                assign[nextpoint[0]][nextpoint[1]]=0
                binmass=nmass
                centroid=ncentroid
                bin_accretion.validateappend(viablecell,(nextpoint[0]+1,nextpoint[1]),assign)
                bin_accretion.validateappend(viablecell,(nextpoint[0]-1,nextpoint[1]),assign)
                bin_accretion.validateappend(viablecell,(nextpoint[0],nextpoint[1]+1),assign)
                bin_accretion.validateappend(viablecell,(nextpoint[0],nextpoint[1]-1),assign)
                accrete= len(viablecell)>0
            else: 
                accrete=False
        success=0.8
        if binmass/(target**2)<success:
            rebinlist.append(current)
            rcentroids.append(centroid)
        else:
            binlist.append(current)
            bcentroids.append(centroid)
            
        supercentroid=((supercentroid[0]*supermass+centroid[0]*binmass)/(supermass+binmass),(supercentroid[1]*supermass+centroid[1]*binmass)/(supermass+binmass))
        supermass=supermass+binmass
        viable.extend([v for v in viablecell if not v in viable])
        print(np.count_nonzero(assign == -1))

    print("Redistribution time")

    ## Then we assign each unsuccessfully binned pixel to a successfully bin
    #functions.generate_wvt3(binlist+rebinlist,signal,var,np.full(len(binlist)+len(rebinlist),1),displayWVT=True)
    functions.redistribute(binlist,rebinlist,bcentroids,density)
    #functions.generate_wvt3(binlist,signal,var,np.full(len(binlist),1),displayWVT=True)
    ## At this point, binlist should contain all of the original points
    ## Now I want to iterate through binlist to get the list of generators. This is really what this was for
    ## Though now is as good of a time as any to create the CVT
    
    
    return switch(version,target,binlist,signal,var)

def snext_iteration2(version,target,signal,var,geocarray,scalearray,weighting=True):
    ## generate all of the pixels to be slotted into a bin determined by the generators
    #allpix=[(y,x) for y in range(signal.shape[0]) for x in range(signal.shape[1]) ]
    
    #assign=[-1 for _ in range(len(allpix))]
    assign=np.full_like(signal,-1,dtype=int)
    viable=[]
    for g in range(len(geocarray)):
        point=(int(geocarray[g][0]),int(geocarray[g][1]))
        try:
            assign[point[0]][point[1]]=g
        except:
            print(point)
            print(geocarray[g])
            print(g)
            raise NameError("ouchi")
        viable.append([])
        wvt_iteration.append_validate((point[0]+1,point[1]),viable[g],assign)
        wvt_iteration.append_validate((point[0]-1,point[1]),viable[g],assign)
        wvt_iteration.append_validate((point[0],point[1]+1),viable[g],assign)
        wvt_iteration.append_validate((point[0],point[1]-1),viable[g],assign)
        #print(str(int(g*100/len(geocarray)))+" percent done with init pass")
    while wvt_iteration.checkneg(assign) or wvt_iteration.viabunempty(viable):
        for g in range(len(geocarray)):
            prune=True
            while prune and len(viable[g])>0:
                point=viable[g].pop(0)
                if assign[point[0]][point[1]]==g:
                    prune=True
                else:
                    prune=False
            if len(viable[g])>0:
                if assign[point[0]][point[1]]==-1:
                    assign[point[0]][point[1]]=g
                    wvt_iteration.append_validate((point[0]+1,point[1]),viable[g],assign)
                    wvt_iteration.append_validate((point[0]-1,point[1]),viable[g],assign)
                    wvt_iteration.append_validate((point[0],point[1]+1),viable[g],assign)
                    wvt_iteration.append_validate((point[0],point[1]-1),viable[g],assign)
                else:
                    if ((geocarray[g][0]-point[0])**2+(geocarray[g][1]-point[1])**2)/(scalearray[g]**2)<((geocarray[assign[point[0]][point[1]]][0]-point[0])**2+(geocarray[assign[point[0]][point[1]]][1]-point[1])**2)/(scalearray[assign[point[0]][point[1]]]**2):
                        assign[point[0]][point[1]]=g
                        wvt_iteration.append_validate((point[0]+1,point[1]),viable[g],assign)
                        wvt_iteration.append_validate((point[0]-1,point[1]),viable[g],assign)
                        wvt_iteration.append_validate((point[0],point[1]+1),viable[g],assign)
                        wvt_iteration.append_validate((point[0],point[1]-1),viable[g],assign)
                    else:
                        pass
    binlist=[ [] for _ in range(len(geocarray)) ]
    for j in range(len(assign)):
        for i in range(len(assign[0])):
            binlist[assign[j][i]].append((j,i))
    
    return switch(version,target,binlist,signal,var)

def siteration_func(version,target,signal,var,geocarray,scalearray,epsilon,weighting=True,displaywvt=False):
    wvt=np.zeros_like(signal)
    target=5
    start=time.time()
    ## have to manually kill terminal is does not converge
    difference=2*epsilon

    repeat=True
    numit=0

    while repeat:
        print("another iteration")
        wvt2=np.copy(wvt)
        binlist,geocarray,scalearray=snext_iteration2(version,target,signal,var,geocarray,scalearray)
        #wvt=functions.generate_wvt(binlist,signal,displayWVT=displaywvt)
        wvt,ston=functions.generate_wvt3(binlist,signal,var,scalearray)

        if epsilon<0:
            numit+=1
            if numit+epsilon>=0:
                repeat=False
        else:
            difference=np.sqrt(np.sum((wvt-wvt2)**2))
            print("dif",difference)
            repeat=difference>epsilon

    print("elapsed time "+str(time.time()-start))
    #functions.generate_wvt3(binlist,signal,var,scalearray,displayWVT=True)
    return binlist

def smaketargetscatter(ax,target,binlist,signal,var,version,signdensity):
    
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
    
    ax.plot(grads,gstons,linewidth=0,marker=".",label="ver"+str(version)+"initsign:"+str(signdensity))
    

if __name__ == "__main__":
    wcsx,signal,var,sourcedir,objname=bin_accretion.initialize(enternew=True)
    objname=main.getname("_".join(objname.split("_")[:-1]))
    sourcedir="/".join(sourcedir.split("/")[:-1])
    #objname="J024815-081723"
    ##empty for directly into sourcedir
    target=main.gettarget()
    
    if target<0:
        target=-target

    fig,ax2=plt.subplots()
    
    ax2.plot([0,100],[target,target],linestyle="dashed",color="navy")
    binlistlist={}

    for v in [3]:
        for b in [True,False]:
            binlist,init_generators,init_scalelengths=switch_cc_accretion(v,b,signal,var,target)
            binlist=siteration_func(v,target,signal,var,init_generators,init_scalelengths,-10)
            smaketargetscatter(ax2,target,binlist,signal,var,v,b)
            binlistlist[str(v)+str(b)]=binlist
    for key in binlistlist.values():
        print("This is: "+key)
        functions.generate_wvt3(binlistlist[key],signal,var,np.full(len(binlistlist),1),displayWVT=True)
    
    ax2.set_xlabel("generator radius from center")
    ax2.set_ylabel("bin signal to noise")
    plt.legend()
    plt.show()