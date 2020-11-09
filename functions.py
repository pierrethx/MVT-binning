import numpy as np
import matplotlib.pyplot as plt 
import tkinter
from tkinter.filedialog import askopenfilename
from matplotlib.widgets import Cursor
from astropy.io import fits

import scipy.spatial as sp 

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


def checktoadd(point,binn,weightmap,target,binmass):
    ## There are three conditions we have to check for. All must be met.
    ## First is the 'topological' condition. Must be adjacent.
    if not (( (point[0]+1,point[1]) in binn)or( (point[0]-1,point[1]) in binn)or( (point[0],point[1]+1) in binn)or( (point[0],point[1]-1) in binn)):
        return False
        ## It's unelegant but basically I'm checking that the pixel shares a side with the bin
    
    ## Next condition is 'morphological', that it does not exceed a certain roundness factor set at
    Rmax=0.3
    ## This is suggested by Cappellari and Copin 2003.
    ## Where we need the furthest pixel from the new center
    newbinn=binn+[point]
    ncentroid, nmass=weighted_centroid(newbinn,weightmap)
    rmax=sp.distance.cdist([ncentroid],newbinn).max()
    ## And the area, which is simply the number of pixels in the bin. So
    R=rmax*np.sqrt(np.pi/len(newbinn))-1
    if R>Rmax:
        return False
    ## Lastly the uniformity condition. If the next point brings the S/N further from target we reject:
    if np.abs(binmass-target**2)<np.abs(nmass-target**2):
        return False
    ## otherwise it has passed all three tests. Huzzah!
    return True

def redistribute(binlist,rebinlist,binfo,weightmap):
    for bindex in range(len(rebinlist)):
        for poindex in range(len(rebinlist[bindex])):
            ## Finding the index of the bin with the closest center
            centroind=closest_index(rebinlist[bindex][poindex],binfo,weightmap)
            ## We add the point to that bin
            ## Don't bother with binfo/binm, last used here:
            binlist[centroind].append(rebinlist[bindex][poindex])

def calculate_SN(binn,sigmap,varmap):
    numerator=0
    denominator=0
    for tupple in binn:
        numerator=numerator+sigmap[tupple[0]][tupple[1]]
        denominator=denominator+varmap[tupple[0]][tupple[1]]
    SN=numerator/np.sqrt(denominator)
    return SN

def calculate_signal(binn,sigmap):
    numerator=0
    for tupple in binn:
        numerator=numerator+sigmap[tupple[0]][tupple[1]]
    return numerator

def calculate_scales(target,binlist,signal,var):
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
    '''
            if not np.isnan(delta):
                scalelengths.append(delta)
                geomcentres.append(geoc)
                binlist2.append(binlist[bindex])
            else:
                print("DELETED")
                unbinlist.append(binlist[bindex])
    redistribute(binlist2,unbinlist,geomcentres,signal)
    '''
        

    geocarray=np.array(geomcentres)
    scalearray=np.array(scalelengths)
    
    
    return binlist2,geocarray,scalearray

def calculate_cvt(target,binlist,signal,var):

    geomcentres=[]
    for bindex in range(len(binlist)):
        if len(binlist[bindex])==0:
            print("issue aaaahhhhhh")
            geomcentres.append((0,0))
        else:
            geoc=geometric_center(binlist[bindex])
            geomcentres.append(geoc)
            
    geocarray=np.array(geomcentres)
    return geocarray

def generate_wvt(binlist,signal,displayWVT=False):
    wvt=np.zeros_like(signal)
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
        image=ax.imshow(wvt,cmap="cubehelix")
        fig.colorbar(image)
        plt.show()
    return wvt

def generate_wvt2(binlist,signal,var,displayWVT=False):
    wvt=np.zeros_like(signal)
    ston=np.zeros_like(signal)
    for bindex in range(len(binlist)):
        sig=calculate_signal(binlist[bindex],signal)
        StoN=calculate_SN(binlist[bindex],signal,var)
        for point in binlist[bindex]:
            wvt[point[0]][point[1]]=sig/len(binlist[bindex])
            ston[point[0]][point[1]]=StoN
    if displayWVT:
        fig,ax=plt.subplots()
        def onclick(event):
            x1,y1=event.xdata,event.ydata
            for binn in binlist:
                if (int(y1+.5),int(x1+.5)) in binn:
                    
                    print("other pixels in this bin: ")
                    for tup in binn:
                        print("x: "+str(tup[1])+", y:"+str(tup[0]))
                    print("StoN is for this bin: "+str(ston[binn[0][0]][binn[0][1]]))
                    break
        cursor=Cursor(ax, horizOn=False,vertOn=False,color='red',linewidth=2.0)
        fig.canvas.mpl_connect('button_press_event',onclick)
        image=ax.imshow(wvt,cmap="cubehelix")
        fig.colorbar(image)
        plt.show()
    return wvt, ston

def generate_wvt3(binlist,signal,var,scalearray,displayWVT=False):
    wvt=np.zeros_like(signal)
    ston=np.zeros_like(signal)
    for bindex in range(len(binlist)):
        sig=calculate_signal(binlist[bindex],signal)
        StoN=calculate_SN(binlist[bindex],signal,var)
        for point in binlist[bindex]:
            wvt[point[0]][point[1]]=sig/len(binlist[bindex])
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
                    break
        cursor=Cursor(ax, horizOn=False,vertOn=False,color='red',linewidth=2.0)
        fig.canvas.mpl_connect('button_press_event',onclick)
        image=ax.imshow(ston,cmap="cubehelix")
        fig.colorbar(image)
        plt.show()
    return wvt, ston 

def generate_wvt4(binlist,signal,var,scalearray,displayWVT=False):
    wvt=np.zeros_like(signal)
    ston=np.zeros_like(signal)
    for bindex in range(len(binlist)):
        sig=calculate_signal(binlist[bindex],signal)
        StoN=calculate_SN(binlist[bindex],signal,var)
        for point in binlist[bindex]:
            wvt[point[0]][point[1]]=sig/len(binlist[bindex])
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
                    break
        cursor=Cursor(ax, horizOn=False,vertOn=False,color='red',linewidth=2.0)
        fig.canvas.mpl_connect('button_press_event',onclick)
        image=ax.imshow(wvt,cmap="cubehelix")
        fig.colorbar(image)
        plt.show()
    return wvt, ston         

def numdif(x,y):
    diff=[]
    for i in range(len(y)):
        if i==0 or i==len(y)-1:
            diff.append(0)
        else:
            diff.append((y[i+1]-y[i-1])/(x[i+1]-x[i-1]))
    return np.array(diff)

def smoother(x,n):
    ## no smoothing is n=0
    y=[]
    for j in range(len(x)):
        upper=j+n+1
        lower=j-n
        if lower<0:
            lower=0
        if upper>len(x):
            upper=len(x)
        y.append(np.average(x[lower:upper]))
    return y
                

def lerp(x):
    arr=[]
    arr.append(x[0])
    for i in range(1,len(x)):
        arr.append((x[i]+x[i-1])/2)
        arr.append(x[i])
    return np.array(arr)

def assign(binlist,target,ston,signal):
    assign=np.zeros_like(signal)
    binlist2=binlist.copy()
    np.random.shuffle(binlist2)
    g=1
    for i in range(len(binlist2)):
        for k in binlist2[i]:
            if ston[k[0]][k[1]]<0.5*target:
                assign[k[0]][k[1]]=0
            else:
                assign[k[0]][k[1]]=g
        g=g+1
    return assign