import numpy as np
import matplotlib.pyplot as plt 
import tkinter
from tkinter.filedialog import askopenfilename
from matplotlib.widgets import Cursor
from astropy.io import fits
import random

import scipy.spatial as sp 

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

def reverseassign(map):
    m=int(np.nanmax(map))+1
    
    binlist=[[] for r in range(m)]
    bbl=[]
    for y in range(len(map)):
        for x in range(len(map[0])):
            binlist[int(map[y,x])].append((y,x))
            bbl.append((y,x))  
    return binlist,bbl

def closest_point(point,pointset,weightmap):
    ## as name suggests, we are calculating the closest point. Returns index
    try:
        ind=sp.distance.cdist([point],pointset).argmin()
    except:
        print(point)
        print(pointset)
        print('cdist failed, using manually calculating distance')
        pointsetlist=np.asarray(list(map(lambda x:[x[0],x[1]],pointset)))
        pointvect=np.array([point[0],point[1]])
        differences=pointsetlist-pointvect
        ind=np.einsum('ij,ij->i',differences,differences).argmin()
        
    ## if multiple points on list fulfill this criteria multiple will be returned as ind
    ## thus we need to reduce to one entry.
    ## for now, it chooses the entry that has the highest intrinsic S/N
    if isinstance(ind,list):
        temp=list(map(lambda x:weightmap[x[0]][x[1]],ind))
        return pointset[max(temp)]
    else:
        return pointset[ind]

def closest_index(point,pointset,weightmap):
    ## Closest point but for the index
    try:
        ind=sp.distance.cdist([point],pointset).argmin()
    except:
        print('cdist failed, using manually calculating distance')
        pointsetlist=np.asarray(list(map(lambda x:[x[0],x[1]],pointset)))
        pointvect=np.array([point[0],point[1]])
        differences=pointsetlist-pointvect
        ind=np.einsum('ij,ij->i',differences,differences).argmin()
        
    ## if multiple points on list fulfill this criteria multiple will be returned as ind
    ## thus we need to reduce to one entry. This particular function is used to split unsuc. bins
    ## so for this one, have it choose the point with the least intrinsic S/N
    ## It should rarely come to this but ideally if it does, it would make bins more even
    if isinstance(ind,list):
        temp=list(map(lambda x:weightmap[x[0]][x[1]],ind))
        return min(temp)
    else:
        return ind

def scaled_closest(point,pointset,scaleset):
    ## Closest point but for the index
    pointsetlist=np.asarray(list(map(lambda x:[x[0],x[1]],pointset)))
    pointvect=np.array([point[0],point[1]])
    differences=pointsetlist-pointvect
    sqrdistances=np.einsum('ij,ij->i',differences,differences)
    return np.divide(sqrdistances,scaleset**2).argmin()

def weighted_centroid(pointset,weightmap):
    ## We are calculating the weighted centroid for a set of points. Returns index and set mass.
    if len(pointset)==1:
        return pointset[0],weightmap[pointset[0][0]][pointset[0][1]]
    else:
        mass=sum(list(map(lambda x:weightmap[x[0]][x[1]],pointset)))
        indlist=list(map(lambda y:tuple(map(lambda x:weightmap[y[0]][y[1]]*x/mass,y)),pointset))
        ind=tuple(map(sum,zip(*indlist)))
        return ind, mass

def geometric_center(pointset):
    ## We are calculating the geometric center for a set of points. Returns index.
    if len(pointset)==1:
        return pointset[0]
    else:
        ind=tuple(map(lambda x:sum(x)/len(x),zip(*pointset)))
        return ind

def calculate_signal(binn,sigmap):
    ## sums the signal in a bin from a list of the tuples and the signalmap
    numerator=0
    for tupple in binn:
        numerator=numerator+sigmap[tupple[0]][tupple[1]]
    return numerator

def calculate_SN(binn,sigmap,varmap):
    ## as defined in CC03 and DS06
    numerator=0
    denominator=0
    for tupple in binn:
        numerator=numerator+sigmap[tupple[0]][tupple[1]]
        denominator=denominator+varmap[tupple[0]][tupple[1]]
    SN=numerator/np.sqrt(denominator)
    return SN

def calculate_cvt(binlist,signal,var):
    ## like calculate scales, but if we were constructing a CVT, scales all=1 so we just need the geometric centers

    geomcentres=[]
    bin2=[]

    unweightedSN=signal/np.abs(var)**0.5
    unweightedSN=np.where(np.isnan(unweightedSN),0,unweightedSN)
    unweightedmass=unweightedSN**2
    #weightedSN=signal**2/np.abs(var)**1.5
    #weightedSN=np.where(np.isnan(weightedSN),0,weightedSN)

    for bindex in range(len(binlist)):
        if len(binlist[bindex])==0:
            print("empty bin, this is being passed over.")
        else:
            #geoc=geometric_center(binlist[bindex])
            geoc,mass=weighted_centroid(binlist[bindex],unweightedmass)
            if(np.isnan(geoc[0]) or np.isnan(geoc[1])) or np.isnan(mass) or mass==0:
                print("issue but we keep on trucking")
                geoc=geometric_center(binlist[bindex])
            #geoc=geometric_center(binlist[bindex])
            geomcentres.append(geoc)
            bin2.append(binlist[bindex])
            
    geocarray=np.array(geomcentres)
    return bin2,geocarray

def calculate_scales(target,binlist,signal,var):
    ## calculates the scales as formulated in DS06.
    geomcentres=[]
    scalelengths=[]
    binlist2=[]
    for bindex in range(len(binlist)):
        if len(binlist[bindex])==0:
            pass
        else:
            StoN=calculate_SN(binlist[bindex],signal,var)
            geoc=geometric_center(binlist[bindex])
            
            ## Define q to be some constant. Acc to Diehl&Statler, should not have effect
            q=np.pi ## for circular bins, which is generally what we are trying to achieve
            delta=np.sqrt(len(binlist[bindex])*target/(q*StoN))
            scalelengths.append(delta)
            geomcentres.append(geoc)
            binlist2.append(binlist[bindex])

    geocarray=np.array(geomcentres)
    scalearray=np.array(scalelengths)
 
    return binlist2,geocarray,scalearray

def generate_wvt(binlist,signal,displayWVT=False):
    ## generates a WVT from a datamap and displays it
    wvt=np.zeros_like(signal,dtype=float)
    for bindex in range(len(binlist)):
        sig=calculate_signal(binlist[bindex],signal)
        for point in binlist[bindex]:
                wvt[point[0]][point[1]]=sig/len(binlist[bindex])
    if displayWVT:
        fig,ax=plt.subplots()
        def onclick(event):
            x1,y1=event.xdata,event.ydata
            for binn in binlist:
                if (int(y1+.5),int(x1+.5)) in binn:
                    print("other pixels in this bin: ")
                    for tup in binn:
                        print("x: "+str(tup[1])+", y:"+str(tup[0]))
                    break
        cursor=Cursor(ax, horizOn=False,vertOn=False,color='red',linewidth=2.0)
        fig.canvas.mpl_connect('button_press_event',onclick)
        image=ax.imshow(wvt,cmap="cubehelix",origin="lower")
        fig.colorbar(image)
        plt.show()
    return wvt

def generate_wvt3(binlist,signal,var,scalearray,maxx=-1,displayWVT=False):
    ## generates a WVT and StoN map from signal and var maps and displays the StoN
    ## and allows you to interact with graph and see the scalelengths
    wvt=np.zeros_like(signal,dtype=float)
    ston=np.zeros_like(signal,dtype=float)
    for bindex in range(len(binlist)):
        sig=calculate_signal(binlist[bindex],signal)
        StoN=calculate_SN(binlist[bindex],signal,var)
        for point in binlist[bindex]:
            wvt[point[0]][point[1]]=sig/len(binlist[bindex])
            if sig==0:
                ston[point[0]][point[1]]=0
            else:
                ston[point[0]][point[1]]=StoN
    if displayWVT:
        fig,ax=plt.subplots()
        def onclick(event):
            x1,y1=event.xdata,event.ydata
            for bint in range(len(binlist)):
                if (int(y1+.5),int(x1+.5)) in binlist[bint]:
                    binn=binlist[bint]
                    
                    print("other pixels in this bin: ")
                    for tup in binn:
                        print("x: "+str(tup[1])+", y:"+str(tup[0]))
                    print("StoN is for this bin: "+str(ston[binn[0][0]][binn[0][1]]))
                    print("scale for this bin: "+str(scalearray[bint]))
                    print("index for this: "+str(bint))
                    print("density for this: "+str(ston[binn[0][0]][binn[0][1]]/len(binn)))
                    break
        cursor=Cursor(ax, horizOn=False,vertOn=False,color='red',linewidth=2.0)
        fig.canvas.mpl_connect('button_press_event',onclick)
        if maxx<0:
            image=ax.imshow(ston,cmap="cubehelix",origin="lower")
        else:
            image=ax.imshow(ston,cmap="cubehelix",vmax=maxx,origin="lower")
        fig.colorbar(image)
        plt.show()
    return wvt, ston 

def justshow(signal):
    fig,ax=plt.subplots()
    image=ax.imshow(signal,cmap="cubehelix",origin="lower")
    fig.colorbar(image)
    plt.show()

if __name__ == "__main__":
    signal=np.full((4,4),1)
    var=np.full((4,4),0.1)
    binlist=[[(3,0),(3,1)],[(2,2),(2,3),(3,2),(3,3)],[(2,0),(2,1),(1,0),(1,1),(0,0),(0,1)],[(1,2)],[(0,2),(0,3),(1,3)]]
    generate_wvt3(binlist,signal,var,[1,1,1,1,1],maxx=-1,displayWVT=False)