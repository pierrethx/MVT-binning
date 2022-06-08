import numpy as np
import matplotlib.pyplot as plt 
import tkinter
from tkinter.filedialog import askopenfilename
from astropy.io import fits
import functions 
import scipy.spatial as sp
import time
from astropy import wcs

#from astropy.utils.data import download_file
def initialize(enternew=True):

    ## First let's select a signal file
    ## first you can select signal, then variance. 
    ## for convenience you can select both files in the first window but they need to have "sig" or "var" in the name
    if enternew:
        validatesignal=True
        validatevar=True
        while validatesignal:
            root=tkinter.Tk()
            root.withdraw()
            placeholder= askopenfilename(message="Select signal",multiple=True)
            root.update()
            root.destroy()
            if type(placeholder) is tuple:
                if len(placeholder)>2:
                    print("Too many files selected.")
                if len(placeholder)==1:
                    if ".fits" in placeholder[0]:
                        signalname=placeholder[0]
                        print("Good file! Moving on")
                        while validatevar:
                            root=tkinter.Tk()
                            root.withdraw()
                            place= askopenfilename(message="Select variance")
                            root.update()
                            root.destroy()
                            print("place: ",place)
                            if ".fits" in place:
                                varname=place
                                print("Good file! Moving on")
                                validatevar=False
                                validatesignal=False
                            elif place.strip()=="":
                                print("var input canceled, back to sig input")
                                validatevar=False
                            else:
                                print("invalid var file type")
                    else:
                        print("Invalid sig file type")
                else:
                    if ".fits" in placeholder[0] and ".fits" in placeholder[1]:
                        place0=placeholder[0].split("/")[-1]
                        place1=placeholder[1].split("/")[-1]
                        if "var" in place0.lower() or "sig" in place1.lower():
                            varname=placeholder[0]
                            signalname=placeholder[1]
                            validatesignal=False
                            print("Good files! Moving on.")
                        else:
                            signalname=placeholder[0]
                            varname=placeholder[1]
                            validatesignal=False
                            print("Good files! Moving on")
                    else:
                        print("Invalid file type")
            else:
                raise NameError("Bye bye, see you later.")  
    else:
        signalname="/Users/pierre/Downloads/image.J024815-081723_icubes.wc.c5008_29.fits"
        varname="/Users/pierre/Downloads/image.J024815-081723_icubes.wc.c5008_29_VAR.fits"

    objname=signalname.split("/")[-1]
    foldername=signalname.split("/")[-2]
    if "unbinned" in foldername.lower():
        sourcedir="/".join(signalname.split("/")[:-2])
        print("unbinned folder name, going back one directory level")
    else:
        sourcedir="/".join(signalname.split("/")[:-1])
    
    ## to get the directory of the file so that we can save new files there

    ##filename is a string locating the selected file
    with fits.open(signalname,checksum=True) as hdul:
        signal=np.flipud(hdul[0].data)
        wcsx=hdul[0].header
    with fits.open(varname,checksum=True) as hdul:
        var=np.flipud(hdul[0].data)
        
    return wcsx,signal,var,sourcedir,objname

## this function allows one to select multiple sets of files at once. it is better
def minitialize():
    wcsxlist=[]
    signallist=[]
    varlist=[]
    sourcedirlist=[]
    objnamelist=[]
    
    validatesignal=True
    while validatesignal:
        root=tkinter.Tk()
        root.withdraw()
        placeholder= askopenfilename(message="Select signal",multiple=True)
        root.update()
        root.destroy()
        fails=[]
        if type(placeholder) is tuple:
            if len(placeholder)>1:
                num=0
                placeholder=list(placeholder)
                while len(placeholder)>0:
                    try:
                        if ".fits" in placeholder[0]:
                            i=1
                            print("_".join(placeholder[0].split("_")[:-1]))
                            while not "_".join(placeholder[0].split("_")[:-1]) in placeholder[i]:
                                print(placeholder[i])
                                i+=1

                            if "var" in placeholder[0].lower():
                                signalname=placeholder.pop(i)
                                varname=placeholder.pop(0)
                            else:
                                varname=placeholder.pop(i)
                                signalname=placeholder.pop(0)

                            sourcedir="/".join(signalname.split("/")[:-1])
                            objname=signalname.split("/")[-1]
                            foldername=signalname.split("/")[-2]
                            with fits.open(signalname) as hdul:
                                signal=np.flipud(hdul[0].data)
                                wcsx=hdul[0].header
                            with fits.open(varname) as hdul:
                                var=np.flipud(hdul[0].data)
                            if "unbinned" in foldername.lower():
                                sourcedir="/".join(signalname.split("/")[:-2])
                                print("unbinned folder name, going back one directory level")
                            else:
                                sourcedir="/".join(signalname.split("/")[:-1])
                            wcsxlist.append(wcsx)
                            signallist.append(signal)
                            varlist.append(var)
                            sourcedirlist.append(sourcedir)
                            objnamelist.append(objname)
                            num+=1
                        else:
                            raise NameError("")
                    except:
                        fails.append(placeholder.pop(0))
                if len(fails)>0:
                    print("rejected: "+str(fails))
                print("Good "+str(num)+" pairs of files! Moving on")

            elif len(placeholder)==1:
                validatevar=True
                if ".fits" in placeholder[0]:
                    signalname=placeholder[0]
                    print("Good signal file! Moving on")
                    while validatevar:
                        root=tkinter.Tk()
                        root.withdraw()
                        place= askopenfilename(message="Select variance")
                        root.update()
                        root.destroy()
                        print("place: ",place)
                        if ".fits" in place:
                            varname=place
                            print("Good file! Moving on")
                            validatevar=False
                            validatesignal=False
                        elif place.strip()=="":
                            print("var input canceled, back to sig input")
                            validatevar=False
                        else:
                            print("invalid var file type")
                    sourcedir="/".join(signalname.split("/")[:-1])
                    objname=signalname.split("/")[-1]
                    foldername=signalname.split("/")[-2]
                    with fits.open(signalname) as hdul:
                        signal=np.flipud(hdul[0].data)
                        wcsx=hdul[0].header
                    with fits.open(varname) as hdul:
                        var=np.flipud(hdul[0].data)
                    if "unbinned" in foldername.lower():
                        sourcedir="/".join(signalname.split("/")[:-2])
                        print("unbinned folder name, going back one directory level")
                    else:
                        sourcedir="/".join(signalname.split("/")[:-1])
                    wcsxlist.append(wcsx)
                    signallist.append(signal)
                    varlist.append(var)
                    sourcedirlist.append(sourcedir)
                    objnamelist.append(objname)
                else:
                    print("Invalid sig file type")
            else:
                if ".fits" in placeholder[0] and ".fits" in placeholder[1]:
                    if "var" in placeholder[0].lower() or "sig" in placeholder[1].lower():
                        varname=placeholder[0]
                        signalname=placeholder[1]
                        validatesignal=False
                        print("Good files! Moving on")
                    else:
                        signalname=placeholder[0]
                        varname=placeholder[1]
                        validatesignal=False
                        print("Good files! Moving on")
                    objname=signalname.split("/")[-1]
                    foldername=signalname.split("/")[-2]
                    
                    with fits.open(signalname,checksum=True) as hdul:
                        signal=np.flipud(hdul[0].data)
                        wcsx=hdul[0].header
                    with fits.open(varname,checksum=True) as hdul:
                        var=np.flipud(hdul[0].data)
                    if "unbinned" in foldername.lower():
                        sourcedir="/".join(signalname.split("/")[:-2])
                        print("unbinned folder name, going back one directory level")
                    else:
                        sourcedir="/".join(signalname.split("/")[:-1])
                    wcsxlist.append(wcsx)
                    signallist.append(signal)
                    varlist.append(var)
                    sourcedirlist.append(sourcedir)
                    objnamelist.append(objname)

                else:
                    print("Invalid file type")
        else:
            validatesignal=False

    return wcsxlist,signallist,varlist,sourcedirlist,objnamelist

## this is used in bin accretion and wvt construction, to check to see that tuples are actually
## members of the image (Since they both construct bins using a wavefront approach, edges need to be checked)
def validateappend(target,candidate,check):
    ## candidate is as (y,x)
    try:
        if candidate[0]<0 or candidate[0]>=len(check) or candidate[1]<0 or candidate[1]>=len(check[0]):
            raise NameError("bluhbluh")
        if check[candidate[0]][candidate[1]]==-1:
            try:
                target.index(candidate)
            except:
                target.append(candidate)
    except:
        pass

## this is the primary bin accretion function
def cc_accretion(signal,var,target):

    ## this is the signal-to-noise for each pixel
    ston=signal/np.sqrt(np.abs(var))
    ## to prevent situations where 0/0=nan by making any place where sig=0, ston=0.
    ston=np.where(np.equal(signal,np.zeros_like(signal)),signal,ston)
    var[signal<=0]=0
    
    ## we set up an assignment array to easier check bin loyalties for pixels
    assign=np.full_like(ston,-1)
            
    #density=ston*np.abs(ston)
    ## density is defined in CC03 as SN squared, sign unimportant since CC03 uses pos data,
    ##  above def would preserve sign, but it doesnt seem like this makes binnings better
    density=ston*ston
    cellsleft=np.count_nonzero(assign == -1)
    print(cellsleft)

    ## we start from pixel with highest signal to noise, assumed to be center of nebula
    ## we also need to keep track of the center of ALL binned pixels. at step 1 this is the only binned pix
    supercentroid=np.unravel_index(ston.argmax(),ston.shape)
    supermass=density[supercentroid[0]][supercentroid[1]]
    
    viable=[]
    validateappend(viable,(supercentroid[0],supercentroid[1]),assign)
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
        validateappend(viablecell,(centroid[0]+1,centroid[1]),assign)
        validateappend(viablecell,(centroid[0]-1,centroid[1]),assign)
        validateappend(viablecell,(centroid[0],centroid[1]+1),assign)
        validateappend(viablecell,(centroid[0],centroid[1]-1),assign)
        binmass=density[centroid[0]][centroid[1]]
        ## do we continue to accrete? Have we reached minimum StoN?
        accrete=binmass<target**2 and len(viablecell)>0
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
            if R<=Rmax and not np.isnan(nmass) and np.abs(binmass-target**2)>=np.abs(nmass-target**2):
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
                validateappend(viablecell,(nextpoint[0]+1,nextpoint[1]),assign)
                validateappend(viablecell,(nextpoint[0]-1,nextpoint[1]),assign)
                validateappend(viablecell,(nextpoint[0],nextpoint[1]+1),assign)
                validateappend(viablecell,(nextpoint[0],nextpoint[1]-1),assign)
                accrete= len(viablecell)>0
            else: 
                ## or, we have chosen to end accretion and exit the loop for this bin
                accrete=False
                
        ## if bin reaches minimum SN of 0.8target, we accept, otherwise it will be redispersed.
        success=0.8
        if binmass/(target**2)<success or np.isnan(binmass):
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

    ## Then we assign each unsuccessfully binned pixel to a successfully bin
    functions.redistribute(binlist,rebinlist,bcentroids,density)
    ## At this point, binlist should contain all of the original points
    ## Now I want to iterate through binlist to get the list of generators. This is really what this was for
    ## Though now is as good of a time as any to create the CVT
    binlist,geocarray,scalearray=functions.calculate_scales(target,binlist,signal,var)
    
    #wvt,ston=functions.generate_wvt3(binlist,signal,var,scalearray,True)
    
    return binlist,geocarray,scalearray
    

if __name__ == "__main__":
    wcsx,signal,var,sourcedir,objname=initialize()
    target=300
    mid=time.time()
    binlist,geocarray,scalearray=cc_accretion(signal,var,target)
    print("elapsed time spread method"+str(time.time()-mid))
    #wvt,ston=functions.generate_wvt2(binlist,signal,var,displayWVT=True)
    wvt,ston=functions.generate_wvt4(binlist,signal,var,[1]*len(binlist),displayWVT=True)
    fig,ax=plt.subplots()
    g=ax.imshow(ston,cmap="cubehelix", vmin=np.nanmin(ston))
    plt.colorbar(g)
    plt.show()
    