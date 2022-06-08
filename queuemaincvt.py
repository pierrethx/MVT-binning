import bin_accretion,main,functions
import numpy as np
import scipy.spatial as sp

def checkneg(assign):
    num=np.count_nonzero(assign==-1)
    #print(num)
    if num>0:
        return True
    else:
        return False

def viabunempty(viable):
    for binn in viable:
        if len(binn)>0:
            return True
    return False

def append_validate(candidate,target,check):
    try:
        check[candidate[0]][candidate[1]]
        if candidate[0]<0 or candidate[0]>=len(check) or candidate[1]<0 or candidate[1]>=len(check[0]):
            raise NameError("brrr overflow bro")
        target.append(candidate)
    except:
        pass

def cc_accretion(signal,var,target):

    ## this is the signal-to-noise for each pixel
    ston=signal/np.sqrt(np.abs(var))
    ## to prevent situations where 0/0=nan by making any place where sig=0, ston=0.
    ston=np.where(np.equal(signal,np.zeros_like(signal)),signal,ston)

    var[signal<=0]=0
    
    ## we set up an assignment array to easier check bin loyalties for pixels
    assign=np.full_like(ston,-1)
            
    density=ston*np.abs(ston)
    ## density is defined in CC03 as SN squared, sign unimportant since CC03 uses pos data,
    ##  above def would preserve sign, but it doesnt seem like this makes binnings better
    #density=ston*ston
    cellsleft=np.count_nonzero(assign == -1)
    print(cellsleft)

    ## we start from pixel with highest signal to noise, assumed to be center of nebula
    ## we also need to keep track of the center of ALL binned pixels. at step 1 this is the only binned pix
    supercentroid=np.unravel_index(ston.argmax(),ston.shape)
    supermass=density[supercentroid[0]][supercentroid[1]]

    minsize=30
    #minsize=1
    
    viable=[]
    bin_accretion.validateappend(viable,(supercentroid[0],supercentroid[1]),assign)
    binlist=[]
    bcentroids=[]
    rcentroids=[]
    rebinlist=[]
    
    ## in case this borks np.count_nonzero(array == value) to return -1 in assign
    ## the process is as follows. we have a list of binned pixels, and a list of adjacent, "viable" pixels
    ## we choose whichever viable pixel is closest to the supercentroid to start the next bin
    ## as we add pixels to new bin, we check to see if accretion continues, then add adj pixels to viablecell
    ## we choose new pixels to add to the bin from viablecell
    ## once bin accretion stops (violates conditions), viablecell is merged into viable and a new pixel is selected to start next bin
    while len(viable)>0:
        ## picking new point to start next bin
        centroid=functions.closest_point(supercentroid,viable,density)
        ## begin with empty cell. Add pixel to new cell and remove from viable
        current=[]
        current.append(centroid)
        viable.remove(centroid)
        assign[centroid[0]][centroid[1]]=0
        ## introduce list of candidates for bin and populate with adjacent cell, checking to see if it is already assigned
        viablecell=[]
        bin_accretion.validateappend(viablecell,(centroid[0]+1,centroid[1]),assign)
        bin_accretion.validateappend(viablecell,(centroid[0]-1,centroid[1]),assign)
        bin_accretion.validateappend(viablecell,(centroid[0],centroid[1]+1),assign)
        bin_accretion.validateappend(viablecell,(centroid[0],centroid[1]-1),assign)
        binmass=density[centroid[0]][centroid[1]]
        ## do we continue to accrete? Have we reached minimum StoN?
        accrete=(len(current)<minsize or binmass<target**2 )and len(viablecell)>0
        while accrete:
            nextpoint=functions.closest_point(centroid,viablecell,density)
            ## Rmax is maxumum roundness. Defined in CC03
            Rmax=0.3
            newbin=current+[nextpoint]
            nmass=binmass+density[nextpoint[0]][nextpoint[1]]
            ## replacing weighted centroids with geometric centers as suggested by Diehl 
            ## in Cappellari's implementation to address negative data
            ncentroid=functions.geometric_center(newbin)
            rmax=sp.distance.cdist([ncentroid],newbin).max()
            R=rmax*np.sqrt(np.pi/len(newbin))-1
            #if R<=Rmax and np.abs(binmass-target**2)>np.abs(nmass-target**2): here we are saying that the new mass brings the bin closer to the target
            #if binmass >= nmass then the new pixel is either 0 or negative (some prescription). We should just add
            if len(newbin)<minsize or binmass>=nmass or (R<=Rmax and not np.isnan(nmass) and np.abs(binmass-target**2)>=np.abs(nmass-target**2)):
                ## pixel has been chosen to be accreted. we remove from viablecell (and viable) and add to bin
                current.append(nextpoint)
                viablecell.remove(nextpoint)
                if viable.count(nextpoint)>0:
                    viable.remove(nextpoint)
                assign[nextpoint[0]][nextpoint[1]]=0
                ## we update the binmasses, which will be used to calculate the scales and stuff
                binmass=nmass
                centroid=ncentroid
                ## then add new pixels to viablecell for continued accretion
                bin_accretion.validateappend(viablecell,(nextpoint[0]+1,nextpoint[1]),assign)
                bin_accretion.validateappend(viablecell,(nextpoint[0]-1,nextpoint[1]),assign)
                bin_accretion.validateappend(viablecell,(nextpoint[0],nextpoint[1]+1),assign)
                bin_accretion.validateappend(viablecell,(nextpoint[0],nextpoint[1]-1),assign)
                accrete= len(viablecell)>0
            else: 
                ## or, we have chosen to end accretion and exit the loop for this bin
                accrete=False
                
        ## if bin reaches minimum SN of 0.8target, we accept, otherwise it will be redispersed.
        success=0.8
        if binmass/(target**2)<success or len(current)<minsize-1 or np.isnan(binmass):
            rebinlist.append(current)
            rcentroids.append(centroid)
        else:
            binlist.append(current)
            bcentroids.append(centroid)
            
        ## construct centroid of all binned pixels to find start of next bin. and combine viablecell with viable
        if (supermass+binmass==0):
            pass
        else:
            supercentroid=((supercentroid[0]*supermass+centroid[0]*binmass)/(supermass+binmass),(supercentroid[1]*supermass+centroid[1]*binmass)/(supermass+binmass))
        supermass=supermass+binmass
        viable.extend([v for v in viablecell if not v in viable])
        print(np.count_nonzero(assign == -1))

    print("Redistribution time")

    #wvt,ston=functions.generate_wvt3(binlist+rebinlist,signal,var,np.full(len(binlist)+len(rebinlist),1),True)
    #wvt,ston=functions.generate_pass(binlist+rebinlist,signal,var,True)
    ## Then we assign each unsuccessfully binned pixel to a successfully bin
    functions.redistribute(binlist,rebinlist,bcentroids,density)
    ## At this point, binlist should contain all of the original points
    ## Now I want to iterate through binlist to get the list of generators. This is really what this was for
    ## Though now is as good of a time as any to create the CVT
    signal2,var2=makecopy(signal,var)
   #binlist,geocarray=functions.calculate_cvt(target,binlist,signal2,var2)
    binlist,geocarray,scalearray=functions.calculate_scales(target,binlist,signal,var)
    #wvt,ston=functions.generate_wvt3(binlist,signal,var,np.full(len(binlist,1)),True)
    #wvt,ston=functions.generate_pass(binlist,signal,var,True)
    
    return binlist,geocarray

def next_iteration2(target,signal,var,geocarray,weighting=True):
    ## generate all of the pixels to be slotted into a bin determined by the generators
    #allpix=[(y,x) for y in range(signal.shape[0]) for x in range(signal.shape[1]) ]
    
    #assign=[-1 for _ in range(len(allpix))]
    assign=np.full_like(signal,-1,dtype=int)
    viable=[]
    scalearray=np.full(len(geocarray),1)
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
        append_validate((point[0]+1,point[1]),viable[g],assign)
        append_validate((point[0]-1,point[1]),viable[g],assign)
        append_validate((point[0],point[1]+1),viable[g],assign)
        append_validate((point[0],point[1]-1),viable[g],assign)
        #print(str(int(g*100/len(geocarray)))+" percent done with init pass")
    while checkneg(assign) or viabunempty(viable):
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
                    append_validate((point[0]+1,point[1]),viable[g],assign)
                    append_validate((point[0]-1,point[1]),viable[g],assign)
                    append_validate((point[0],point[1]+1),viable[g],assign)
                    append_validate((point[0],point[1]-1),viable[g],assign)
                else:
                    if ((geocarray[g][0]-point[0])**2+(geocarray[g][1]-point[1])**2)/(scalearray[g]**2)<((geocarray[assign[point[0]][point[1]]][0]-point[0])**2+(geocarray[assign[point[0]][point[1]]][1]-point[1])**2)/(scalearray[assign[point[0]][point[1]]]**2):
                        assign[point[0]][point[1]]=g
                        append_validate((point[0]+1,point[1]),viable[g],assign)
                        append_validate((point[0]-1,point[1]),viable[g],assign)
                        append_validate((point[0],point[1]+1),viable[g],assign)
                        append_validate((point[0],point[1]-1),viable[g],assign)
                    else:
                        pass
    binlist=[ [] for _ in range(len(geocarray)) ]
    for j in range(len(assign)):
        for i in range(len(assign[0])):
            binlist[assign[j][i]].append((j,i))

    #binlist,geocarray=functions.calculate_cvt(target,binlist,signal,var)
    binlist,geocarray,scalearray=functions.calculate_scales(target,binlist,signal,var)
    return binlist,geocarray

def iteration_func(target,signal,var,geocarray,epsilon,displaywvt=False):
    wvt=np.zeros_like(signal)
    
    ## have to manually kill terminal is does not converge
    difference=2*epsilon
    diflist=[]
    repeat=True
    numit=0
    maxit=50
    iterator=0

    while repeat:
        print("another iteration")
        iterator+=1
        wvt2=np.copy(wvt)
        binlist,geocarray=next_iteration2(target,signal,var,geocarray)
        wvt,ston=functions.generate_wvt3(binlist,signal,var,np.full(len(binlist),1),displayWVT=displaywvt)
        difference=np.sqrt(np.sum((wvt-wvt2)**2)/np.sum(var))
        print("dif",difference)

        diflist.append(difference)
        if epsilon<=0:
            numit+=1
            if numit+epsilon>=0:
                repeat=False
        else:
            repeat=difference>epsilon
        if maxit<=iterator:
            repeat=False

    return binlist,diflist

def mainfunc(signal,var,target,displayWVT=False,epsilon=10):
    ## this is the most important function
    ## first we do bin accretion. then iteration. the wvt and vwvt is not necessary here but im scared to remove it
    binlist,init_generators=cc_accretion(signal,var,target)
    #wvt,ston=functions.generate_wvt3(binlist,signal,var,np.full(len(binlist),1),True)
    binlist,diflist=iteration_func(target,signal,var,init_generators,epsilon,displaywvt=displayWVT)
    #wvt,ston=functions.generate_pass(binlist,signal,var,True)

    return binlist,diflist

def makecopy(signal,var):
    signal2=np.copy(signal)
    var2=np.copy(var)
    var2[signal2<=0]=1e10
    signal2[signal2<=0]=0
    return signal2,var2

if __name__ == "__main__":
    targhold=0
    targlist=[]
    wcsxlist,siglist,varlist,sourcelist,objlist=bin_accretion.minitialize()
    contqueue=True
    while contqueue:
        try:
            target=main.gettarget(targhold)
            targlist.append(target)
        except:
            contqueue=False
    print("Files loaded!")
    for i in range(len(sourcelist)):
        for m in range(len(targlist)):
            
            wcsx=wcsxlist[i]
            signal=siglist[i]
            var=varlist[i]
            objname="_".join(objlist[i].split("_")[:-1])
            target=targlist[m]
            sourcedir=sourcelist[i]
            print(sourcedir)
            subfolder=main.makesubfolder(sourcedir,target)
            #main.saveston(wscxlist[i],siglist[i],varlist[i],sourcelist[i],objlist[i],subfolder="unbinned")

            ## this is us applying the data mask. We block out negative values and then run the binning algorithm on it
            signal2,var2=makecopy(signal,var)
            eps=-10
            binlist,diflist=mainfunc(signal2,var2,target,displayWVT=False,epsilon=-10)
            
            #main.saveblockoutfits(targlist[m],binlist,wscxlist[i],siglist[i],varlist[i],objlist[i],sourcelist[i],subfolder=subfolder)
            #wvt,ston=functions.generate_wvt2(binlist,siglist[i],varlist[i])
            ## then we apply the bins to the actual data and save all our files.
            wvt,ston=functions.generate_wvt3(binlist,signal,var,np.full(len(binlist),1))
            vwvt=functions.generate_wvt(binlist,var)
            main.saveiteratedfits(target,wcsx,wvt,vwvt,objname+"_c",sourcedir,subfolder=subfolder,weighting=False)
            functions.convergence(eps,diflist,sourcedir,objname+"_c",subfolder=subfolder)
            main.saveblockoutfits(target,ston,wcsx,wvt,vwvt,objname+"_c",sourcedir,subfolder=subfolder,weighting=False,check=1)
            #main.saveblockoutoldfits(target,ston,wcsx,wvt,vwvt,objname,sourcedir,subfolder=subfolder)
            main.saveston(wcsx,ston,sourcedir,objname+"_c",subfolder=subfolder)
            assign=functions.assign(binlist,target,ston)
            main.saveassign(wcsx,assign,sourcedir,objname+"_c",subfolder=subfolder)
    print("Bye bye")