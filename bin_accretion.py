import numpy as np
import matplotlib.pyplot as plt 
import tkinter
from tkinter.filedialog import askopenfilename
from astropy.io import fits
import functions 
import scipy.spatial as sp

#from astropy.utils.data import download_file
def initialize():
    enternew=False

    ##First let's select a signal file
    if enternew:
        ## this is to obtain a file
        tkinter.Tk().withdraw()
        signalname = askopenfilename()
    else:
        signalname="/Users/pierre/Downloads/image.J024815-081723_icubes.wc.c5008_29.fits"

    if enternew:
        ## this is to obtain a file
        tkinter.Tk().withdraw()
        varname = askopenfilename()
    else:
        varname="/Users/pierre/Downloads/image.J024815-081723_icubes.wc.c5008_29_VAR.fits"


    sourcedir="/".join(signalname.split("/")[:-1])
    ## to get the directory of the file so that we can save new files there

    ##filename is a string locating the selected file
    with fits.open(signalname) as hdul:
        signal=np.abs(hdul[0].data)
    with fits.open(varname) as hdul:
        var=np.abs(hdul[0].data)
    return sourcedir,signal,var

def cc_accretion(signal,var,target):

    ## this is the signal-to-noise for each pixel
    ston=signal/np.sqrt(var)

    '''
    ## To show the signal to noise map
    fig,ax=plt.subplots()
    image=ax.imshow(ston,cmap="cubehelix")
    fig.colorbar(image)
    plt.show()

    quit()
    '''

    ## we replace each unusable value with 0. These unusable values would occur when 
    ## the variance is 0, which should not occur for real data
    ## also i feel like there is a faster way to do this step
    unbin=[]
    for y in range(len(ston)):
        for x in range(len(ston[y])):
            if np.isnan(ston[y][x]):
                ston[y][x]=1e6
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
        supercentroid=(sum([binm[i]*binfo[i][0]/supermass for i in range(len(binfo))]),sum([binm[i]*binfo[i][1]/supermass for i in range(len(binfo))]))
        newpoint=functions.closest_point(supercentroid,unbin,density)
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
    wvt,geocarray,scalearray=functions.calculate_cvt(target,binlist,signal,var,wvt)
    return wvt,geocarray,scalearray
    

if __name__ == "__main__":
    sourcedir,signal,var=initialize()
    target=5
    wvt,geocarray,scalearray=cc_accretion(signal,var,target)
    hdu = fits.PrimaryHDU(wvt)
    hdul = fits.HDUList([hdu])
    hdul.writeto(sourcedir+"/iterated_wvt.fits",overwrite=True)
    np.save(sourcedir+"/gcenters_of_the_image",geocarray)
    np.save(sourcedir+"/scalelengths_first_iteration",scalearray)
    
    fig,ax=plt.subplots()
    image=ax.imshow(wvt,cmap="cubehelix")
    fig.colorbar(image)
    plt.show()
    

