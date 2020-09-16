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
                    else:
                        print("Invalid file type")
            else:
                raise NameError("Bye bye, see you later.")  
    else:
        signalname="/Users/pierre/Downloads/image.J024815-081723_icubes.wc.c5008_29.fits"
        varname="/Users/pierre/Downloads/image.J024815-081723_icubes.wc.c5008_29_VAR.fits"

    sourcedir="/".join(signalname.split("/")[:-1])
    objname=signalname.split("/")[-1]
    ## to get the directory of the file so that we can save new files there

    ##filename is a string locating the selected file
    with fits.open(signalname) as hdul:
        signal=np.flipud(hdul[0].data)
        wcsx=wcs.WCS(hdul[0].header)
    with fits.open(varname) as hdul:
        var=np.flipud(hdul[0].data)
    
    return wcsx,signal,var,sourcedir,objname

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

def cc_accretion(signal,var,target):

    ## this is the signal-to-noise for each pixel
    ston=signal/np.sqrt(np.abs(var))

    
    assign=np.full_like(ston,-1)
            
    #density=ston*np.abs(ston)
    density=ston*ston
    cellsleft=np.count_nonzero(assign == -1)
    print(cellsleft)

    supercentroid=np.unravel_index(ston.argmax(),ston.shape)
    supermass=density[supercentroid[0]][supercentroid[1]]
    
    viable=[]
    validateappend(viable,(supercentroid[0],supercentroid[1]),assign)
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
        validateappend(viablecell,(centroid[0]+1,centroid[1]),assign)
        validateappend(viablecell,(centroid[0]-1,centroid[1]),assign)
        validateappend(viablecell,(centroid[0],centroid[1]+1),assign)
        validateappend(viablecell,(centroid[0],centroid[1]-1),assign)
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
                validateappend(viablecell,(nextpoint[0]+1,nextpoint[1]),assign)
                validateappend(viablecell,(nextpoint[0]-1,nextpoint[1]),assign)
                validateappend(viablecell,(nextpoint[0],nextpoint[1]+1),assign)
                validateappend(viablecell,(nextpoint[0],nextpoint[1]-1),assign)
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
    functions.redistribute(binlist,rebinlist,bcentroids,density)
    ## At this point, binlist should contain all of the original points
    ## Now I want to iterate through binlist to get the list of generators. This is really what this was for
    ## Though now is as good of a time as any to create the CVT
    binlist,geocarray,scalearray=functions.calculate_scales(target,binlist,signal,var)
    
    return binlist,geocarray,scalearray
    

if __name__ == "__main__":
    wcsx,signal,var,sourcedir,objname=initialize(enternew=True)
    target=5
    mid=time.time()
    binlist,geocarray=cc_accretion(signal,var,target)
    print("elapsed time spread method"+str(time.time()-mid))
    wvt,ston=functions.generate_wvt2(binlist,signal,var,displayWVT=True)
    fig,ax=plt.subplots()
    g=ax.imshow(ston,cmap="cubehelix", vmin=np.nanmin(ston),vmax=30)
    plt.colorbar(g)
    plt.show()
    