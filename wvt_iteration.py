import numpy as np
import matplotlib.pyplot as plt 
import tkinter
from tkinter.filedialog import askopenfilename
from astropy.io import fits
import functions,bin_accretion
import scipy.spatial as sp
import time

def initialize():
    enternew=False

    ##First let's select a signal file
    if enternew:
        ## this is to obtain a file
        tkinter.Tk().withdraw()
        geocenters = askopenfilename()
    else:
        geocenters="/Users/pierre/Downloads/gcenters_of_the_image.npy"

    if enternew:
        ## this is to obtain a file
        tkinter.Tk().withdraw()
        scalelengths = askopenfilename()
    else:
        scalelengths="/Users/pierre/Downloads/scalelengths_first_iteration.npy"

    ## to get the directory of the file so that we can save new files there

    ##filename is a string locating the selected file
    geocarray=np.load(geocenters)
    scalearray=np.load(scalelengths)
    return geocarray,scalearray

def next_iteration(target,signal,var,geocarray,scalearray,wvt):
    ## generate all of the pixels to be slotted into a bin determined by the generators
    binlist=[ [] for _ in range(len(geocarray)) ]
    for y in range(len(signal)):
        print("line"+str(y))        
        for x in range(len(signal[y])):
            pixel=(y,x)
            
            ## assign pixel to the generator with the minimum scaled length
            binlist[functions.scaled_closest(pixel,geocarray,scalearray)].append(pixel)
    print("image scanned")
    wvt,geocarray2,scalearray2=functions.calculate_scales(target,binlist,signal,var,wvt)
    for r in range(len(binlist)):
        if binlist[r]==[]:
            geocarray2[r]=geocarray[r]
            scalearray2[r]=scalearray[r]
    return wvt,geocarray2,scalearray2

def next_iteration2(target,signal,var,geocarray,scalearray,wvt):
    ## generate all of the pixels to be slotted into a bin determined by the generators
    allpix=[(y,x) for y in range(signal.shape[0]) for x in range(signal.shape[1]) ]
    
    assign=[-1 for _ in range(len(allpix))]
    viable=[]
    for g in range(len(geocarray)):
        point=(int(geocarray[g][0]),int(geocarray[g][1]))
        assign[allpix.index(point)]=g
        viable.append([])
        append_validate((point[0]+1,point[1]),viable[g],allpix)
        append_validate((point[0]-1,point[1]),viable[g],allpix)
        append_validate((point[0],point[1]+1),viable[g],allpix)
        append_validate((point[0],point[1]-1),viable[g],allpix)
        #print(str(int(g*100/len(geocarray)))+" percent done with init pass")
    print(len(allpix))
    while checkneg(assign):
        for g in range(len(geocarray)):
            prune=True
            while prune and len(viable[g])>0:
                point=viable[g].pop(0)
                indx=allpix.index(point)
                if assign[indx]==g:
                    prune=True
                else:
                    prune=False
            if len(viable[g])>0:
                if assign[indx]==-1:
                    assign[indx]=g
                    append_validate((point[0]+1,point[1]),viable[g],allpix)
                    append_validate((point[0]-1,point[1]),viable[g],allpix)
                    append_validate((point[0],point[1]+1),viable[g],allpix)
                    append_validate((point[0],point[1]-1),viable[g],allpix)
                else:
                    if ((geocarray[g][0]-point[0])**2+(geocarray[g][1]-point[1])**2)/(scalearray[g]**2)<((geocarray[assign[indx]][0]-point[0])**2+(geocarray[assign[indx]][1]-point[1])**2)/(scalearray[assign[indx]]**2):
                        assign[indx]=g
                        append_validate((point[0]+1,point[1]),viable[g],allpix)
                        append_validate((point[0]-1,point[1]),viable[g],allpix)
                        append_validate((point[0],point[1]+1),viable[g],allpix)
                        append_validate((point[0],point[1]-1),viable[g],allpix)
    binlist=[ [] for _ in range(len(geocarray)) ]
    for i in range(len(allpix)):
        binlist[assign[i]].append(allpix[i])
    wvt,geocarray2,scalearray2=functions.calculate_scales(target,binlist,signal,var,wvt)
    for r in range(len(binlist)):
        if len(binlist[r])==0:
            print("empty index"+str(r))
            geocarray2[r]=geocarray[r]
            scalearray2[r]=scalearray[r]
    return wvt,geocarray2,scalearray2

def checkneg(assign):
    '''try:
        assign.index(-1)
        print(min(assign))
        return True
    except ValueError:
        print("failed :(")
        return False'''
    num=assign.count(-1)
    print(num)
    if num>0:
        return True
    else:
        return False

def append_validate(candidate,target,check):
    try:
        check.index(candidate)
        target.append(candidate)
    except:
        pass

def iteration_func(target,signal,var,geocarray,scalearray,epsilon):
    wvt=np.zeros_like(signal)
    target=5
    start=time.time()
    ## have to manually kill terminal is does not converge
    epsilon=10
    difference=2*epsilon

    

    while difference>epsilon:
        print("another iteration")
        wvt2=np.copy(wvt)
        wvt,geocarray,scalearray=next_iteration2(target,signal,var,geocarray,scalearray,wvt)
        
        difference=np.sqrt(np.sum((wvt-wvt2)**2))
        print("dif",difference)

    print("elapsed time "+str(time.time()-start))
    return wvt
    
    

if __name__ == "__main__":
    sourcedir,signal,var=bin_accretion.initialize()
    geocarray,scalearray=initialize()
    target=5
    epsilon=10
    wvt=iteration_func(target,signal,var,geocarray,scalearray,epsilon)

    fig,ax=plt.subplots()
    image=ax.imshow(wvt,cmap="cubehelix")
    fig.colorbar(image)
    plt.show()
    np.save(sourcedir+"/iterated_gcenters2",geocarray)
    np.save(sourcedir+"/iterated_scalelengths2",scalearray)
    hdu = fits.PrimaryHDU(wvt)
    hdul = fits.HDUList([hdu])
    hdul.writeto(sourcedir+"/iterated_wvt2.fits",overwrite=True)