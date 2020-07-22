import numpy as np
import matplotlib.pyplot as plt 
import tkinter
from tkinter.filedialog import askopenfilename
from astropy.io import fits
import functions 
import scipy.spatial as sp
import time

#from astropy.utils.data import download_file
def initialize(enternew=True):

    ## First let's select a signal file
    ## first you can select signal, then variance. 
    ## for convenience you can select both files in the first window but they need to have "sig" or "var" in the name
    if enternew:
        validatesignal=True
        validatevar=True
        while validatesignal:
            tkinter.Tk().withdraw()
            placeholder= askopenfilename(message="Select signal",multiple=True)
            if type(placeholder) is tuple:
                if len(placeholder)>2:
                    print("Too many files selected.")
                if len(placeholder)==1:
                    if ".fits" in placeholder[0]:
                        signalname=placeholder[0]
                        while validatevar:
                            tkinter.Tk().withdraw()
                            place= askopenfilename(message="Select variance")
                            if ".fits" in place:
                                varname=place
                                validatevar=False
                                validatesignal=False
                    else:
                        print("Invalid file type")
                else:
                    if ".fits" in placeholder[0] and ".fits" in placeholder[1]:
                        if "var" in placeholder[0].lower() or "sig" in placeholder[1].lower():
                            varname=placeholder[0]
                            signalname=placeholder[1]
                            validatesignal=False
                        else:
                            signalname=placeholder[0]
                            varname=placeholder[1]
                            validatesignal=False
                    else:
                        print("Invalid file type")
            else:
                raise NameError("Bye bye, see you later.")  
    else:
        signalname="/Users/pierre/Downloads/image.J024815-081723_icubes.wc.c5008_29.fits"
        varname="/Users/pierre/Downloads/image.J024815-081723_icubes.wc.c5008_29_VAR.fits"


    sourcedir="/".join(signalname.split("/")[:-1])
    ## to get the directory of the file so that we can save new files there

    ##filename is a string locating the selected file
    with fits.open(signalname) as hdul:
        signal=np.flipud(hdul[0].data)
    with fits.open(varname) as hdul:
        var=np.flipud(hdul[0].data)
    return sourcedir,signal,var

def cc_accretion0(signal,var,target):

    ## this is the signal-to-noise for each pixel
    ston=signal/np.sqrt(np.abs(var))

    ## we replace each unusable value with 0. These unusable values would occur when 
    ## the variance is 0, which should not occur for real data
    ## also i feel like there is a faster way to do this step
    unbin=[]
    for y in range(len(ston)):
        for x in range(len(ston[y])):
            if np.isnan(ston[y][x]):
                print(var[y][x])
            unbin.append((y,x))
            
    ## Coordinate tuples are inexplicably (y,x). This is how the builtin functions have it
    ## unbinned is the list of all unbinned pixels, which at this stage is all pixels

    '''Bin Accretion algorithm part'''

    ## Funnily enough, we don't use the signal to noise, we actually use the density
    ## For the bin-accretion, density = (S/N)^2 + combines additively. This is adequate.
    density=ston**2
    print(len(unbin))

    ## So we need to start from the pixel with the highest S/N which is at 
    firstcoord=np.unravel_index(ston.argmax(),ston.shape)
    ## We initialize our lists
    ## binlist is the list of bins, which are augmentable lists that hold coord tuples
    ## the first bin is initially just the first point
    binlist=[[firstcoord]]
    rebinlist=[]
    ## and to be proactive, a list to put the one's we'll rebin due to too little mass

    ## with that, we want to remove this point from unbinned
    unbin.remove(firstcoord)
    ## binfo holds the centroid of each (completed) bin so that we don't have to calc it each time
    ## here it is (0,0) and 0 since it will be updated once the first bin is complete
    ## we also have rebinfo for the bins that do not meet success condition
    binfo=[(0,0)]
    rebinfo=[(0,0)]
    ## And the same of the bin masses
    binm=[0]
    rebinm=[0]

    accrete=True
    while accrete:
        accrete=functions.binning(target,binlist,rebinlist,binfo,rebinfo,binm,rebinm,unbin,density)

    ## This is the bin accretion function. It self-iterates. All input objects are mutated
    while len(unbin)>0:
        
        print(len(unbin))

        ## Calculate centroid of all thusfar binned pixels to select new bin start
        supermass=sum(binm)
        print(binm)
        print(supermass)
        raise NameError("Stop!!!")
        supercentroid=(sum([binm[i]*binfo[i][0]/supermass for i in range(len(binfo))]),sum([binm[i]*binfo[i][1]/supermass for i in range(len(binfo))]))
        newpoint=functions.closest_point(supercentroid,unbin,density)

        print("scentroid is "+str(supercentroid))
        print("newpoint is "+str(newpoint),end=" ")
        

        ## then start a new bin. And its corresponding binfo entry
        binlist.append([newpoint])
        unbin.remove(newpoint)
        binfo.append((0,0))
        binm.append(0)
        functions.binning(target,binlist,rebinlist,binfo,rebinfo,binm,rebinm,unbin,density)
        accrete=True
        while accrete:
            accrete=functions.binning(target,binlist,rebinlist,binfo,rebinfo,binm,rebinm,unbin,density)
    print("Redistribution time")

    ## Then we assign each unsuccessfully binned pixel to a successfully bin
    functions.redistribute(binlist,rebinlist,binfo,density)
    ## At this point, binlist should contain all of the original points
    ## Now I want to iterate through binlist to get the list of centroids. This is really what this was for
    ## Though now is as good of a time as any to create the CVT
    wvt=np.zeros_like(signal)
    wvt,geocarray=functions.calculate_cvt(target,binlist,signal,var,wvt)
    return wvt,geocarray

def validateappend(target,candidate,check):
    ## candidate is as (y,x)
    try:
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
            
    density=ston**2
    cellsleft=np.count_nonzero(assign == -1)
    print(cellsleft)

    supercentroid=np.unravel_index(ston.argmax(),ston.shape)
    supermass=density[supercentroid[0]][supercentroid[1]]
    
    viable=[]
    validateappend(viable,(supercentroid[0],supercentroid[1]),assign)
    binlist=[]
    bcentroids=[]
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
            ncentroid, nmass=functions.weighted_centroid(newbin,density)
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
    wvt=np.zeros_like(signal)
    wvt,geocarray=functions.calculate_cvt(target,binlist,signal,var,wvt)
    return wvt,geocarray
    

if __name__ == "__main__":
    sourcedir,signal,var=initialize(enternew=True)
    target=5
    start=time.time()
    
    wvt,geocarray=cc_accretion0(signal,var,target)
    mid=time.time()
    wvt2,geocarray2=cc_accretion(signal,var,target)
    print("elapsed time first method"+str(mid-start))
    print("elapsed time spread method"+str(time.time()-mid))
    fig,ax=plt.subplots(1,1)
    image=ax.imshow(wvt-wvt2,cmap="cubehelix")
    fig.colorbar(image)
    plt.show()
    '''
    hdu = fits.PrimaryHDU(np.flipud(wvt))
    hdul = fits.HDUList([hdu])
    hdul.writeto(sourcedir+"July20/new_style/uniterated_cvt1.fits",overwrite=True)
    #np.save(sourcedir+"/gcenters_of_the_image1",geocarray)
    hdu = fits.PrimaryHDU(np.flipud(wvt2))
    hdul = fits.HDUList([hdu])
    hdul.writeto(sourcedir+"July20/original_recipe/uniterated_wvt2.fits",overwrite=True)
    #np.save(sourcedir+"/gcenters_of_the_image2",geocarray)
    
    fig,ax=plt.subplots()
    image=ax.imshow(wvt,cmap="cubehelix")
    fig.colorbar(image)
    plt.show()'''