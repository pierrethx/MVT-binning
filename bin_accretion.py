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

def cc_accretion(sourcedir,signal,var,target):

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
                ston[y][x]=0
            unbin.append((y,x))
            
    ## Coordinate tuples are inexplicably (y,x). This is how the builtin functions have it
    ## unbinned is the list of all unbinned pixels, which at this stage is all pixels

    '''Bin Accretion algorithm part'''

    ## Funnily enough, we don't use the signal to noise, we actually use the density
    ## For the bin-accretion, density = (S/N)^2 + combines additively. This is adequate.
    density=ston**2

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

    ## This is the bin accretion function. It self-iterates. All input objects are mutated
    while len(unbin)>0:
        accrete=True
        while accrete:
            accrete=functions.binning(target,binlist,rebinlist,binfo,rebinfo,binm,rebinm,unbin,density)
        
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

    '''
    This version creates a map, but I don't wan't this version
    ## And lastly we'll generate an empty set to help make our first step VT:
    pseudovt=np.full(density.shape,0)

    ## This is the bin accretion function. It self-iterates. All input objects are mutated
    while len(unbin)>0:
        accrete=True
        while accrete:
            accrete=functions.binningpseudo(target,binlist,binfo,unbin,density,pseudovt)
        print(len(unbin))

        ## Calculate centroid of all thusfar binned pixels to select new bin start
        supermass=sum(list(map(lambda x:x[2],binfo)))
        supercentroid=(sum(list(map(lambda x:x[0]*x[2],binfo)))/supermass,sum(list(map(lambda x:x[1]*x[2],binfo)))/supermass)
        newpoint=functions.closest_point(supercentroid,unbin,density)
        ## then start a new bin. And its corresponding binfo entry
        binlist.append([newpoint])
        unbin.remove(newpoint)
        binfo.append((0,0,0))
        pseudovt[newpoint[0]][newpoint[1]]=len(binfo)-1
        functions.binningpseudo(target,binlist,binfo,unbin,density,pseudovt)
    binfo[-1]=(newpoint[0],newpoint[1],density[newpoint[0]][newpoint[1]])


    emassvt=np.zeros_like(pseudovt)
    for y in range(len(pseudovt)):
        for x in range(len(pseudovt[y])):
            emassvt[y][x]=binfo[int(pseudovt[y][x])][2]
    ## What we really want is the centroid list. But it is instructive to generate the VT:
    fig,ax=plt.subplots()
    image=ax.imshow(emassvt,cmap="cubehelix")
    #cross=ax.plot(cent[1],cent[0],marker="+",color="k",markersize=8)
    fig.colorbar(image)
    plt.show()
    '''
    ## Then we assign each unsuccessfully binned pixel to a successfully bin
    functions.redistribute(binlist,rebinlist,binfo,density)
    ## At this point, binlist should contain all of the original points
    ## Now I want to iterate through binlist to get the list of centroids. This is really what this was for
    ## Though now is as good of a time as any to create the CVT
    centroids=[]
    scalelengths=[]
    cvt01=np.zeros_like(density)
    numberpoints=0
    for bindex in range(len(binlist)):
        StoN=functions.calculate_SN(binlist[bindex],signal,var)
        cent=functions.weighted_centroid(binlist[bindex],density)[0]
        centroids.append(cent)
        ## Define q to be some constant. Acc to Diehl&Statler, should not have effect
        q=np.pi ## for circular bins, which is generally what we are trying to achieve
        delta=np.sqrt(len(binlist[bindex])*target/(q*StoN))
        scalelengths.append(delta)
        for poindex in range(len(binlist[bindex])):
            point=binlist[bindex][poindex]
            cvt01[point[0]][point[1]]=StoN
            numberpoints=numberpoints+1
    print('number of points: ',numberpoints)
    centarray=np.array(centroids)
    scalearray=np.array(scalelengths)
    fig,ax=plt.subplots()
    image=ax.imshow(cvt01,cmap="cubehelix")
    #ax.plot([centroids[i][1] for i in range(len(centroids))],[centroids[i][0] for i in range(len(centroids))],marker=".",color="black",linewidth=0,markersize=0.5)
    fig.colorbar(image)
    plt.show()
    return centarray,scalearray

if __name__ == "__main__":
    sourcedir,signal,var=initialize()
    target=5
    centarray,scalearray=cc_accretion(sourcedir,signal,var,target)
    np.save(sourcedir+"/centroids_of_the_image",centarray)
    np.save(sourcedir+"/scalelengths_first_iteration",scalearray)

