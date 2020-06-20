import numpy as np
import matplotlib.pyplot as plt 
import tkinter
from tkinter.filedialog import askopenfilename
from astropy.io import fits
import scipy.spatial as sp 

def closest_point(point,pointset,weightmap):
    ## as name suggests, we are calculating the closest point. Returns index
    try:
        ind=sp.distance.cdist([point],pointset).argmin()
    except:
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

def weighted_centroid(pointset,weightmap):
    ## We are calculating the weighted centroid for a set of points. Returns index and set mass.
    mass=sum(list(map(lambda x:weightmap[x[0]][x[1]],pointset)))
    indlist=list(map(lambda y:tuple(map(lambda x:weightmap[y[0]][y[1]]*x/mass,y)),pointset))
    ind=tuple(map(sum,zip(*indlist)))
    return ind, mass


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

def binningpseudo(target,binlist,binfo,unbin,stonmap,tessellation):
    ## immediate stop if there are no more points
    if len(unbin)==0:
        return False
    
    ## Look at the next possible point to be accreted as the closest point to bin centroid
    centroid, mass=weighted_centroid(binlist[-1],stonmap)
    nextpoint=closest_point(centroid,unbin,stonmap)
    ## Apply target condition and accretion condition. If the first is unmet and the second holds,
    ## We add the point to the bin and remove it from unbinned.
    ## The first part is just for expediency. We know that it will fail the uniformity condition.
    if mass<target**2 and checktoadd(nextpoint,binlist[-1],stonmap,target,mass):
        binlist[-1].append(nextpoint)
        unbin.remove(nextpoint)
        tessellation[nextpoint[0]][nextpoint[1]]=len(binfo)-1
        
        return True
    ## If target condition is met or accretion condition unmet, then we terminate bin
    ## Then we move to the next bin
    else:
        ## update binfo with the finished bin centroids and masses
        binfo[-1]=(centroid[0],centroid[1],mass)
        return False

def binning(target,binlist,rebinlist,binfo,rebinfo,binm,rebinm,unbin,stonmap):

    ## We define the success threshold to be some fraction. CC03 suggests:
    ## This will be relevant after accretion has stopped
    success=0.8
    
    ## Look at the next possible point to be accreted as the closest point to bin centroid
    centroid, mass=weighted_centroid(binlist[-1],stonmap)

    ## immediate stop if there are no more points
    if len(unbin)==0:
        ## Assign binfo
        binfo[-1]=(centroid[0],centroid[1])
        binm[-1]=mass
        ## And designate as unsuccessful if condition unmet
        if mass/(target**2)<success:
            rebinfo.append(binfo.pop(-1))
            rebinlist.append(binlist.pop(-1))
            rebinm.append(binm.pop(-1))
        return False
    ## else calculate the next point
    nextpoint=closest_point(centroid,unbin,stonmap)


    ## Apply target condition and accretion condition. If the first is unmet and the second holds,
    ## We add the point to the bin and remove it from unbinned.
    ## The first part is just for expediency. We know that it will fail the uniformity condition.
    if mass<target**2 and checktoadd(nextpoint,binlist[-1],stonmap,target,mass):
        binlist[-1].append(nextpoint)
        unbin.remove(nextpoint)
        
        return True
    ## If target condition is met or accretion condition unmet, then we terminate bin
    ## Then we move to the next bin
    else:
        ## update binfo with the finished bin centroids and masses
        binfo[-1]=(centroid[0],centroid[1])
        binm[-1]=mass
        ## Now is also a time where we can designate binning success relative to the target
        if mass/(target**2)<success:
            rebinfo.append(binfo.pop(-1))
            rebinlist.append(binlist.pop(-1))
            rebinm.append(binm.pop(-1))
        return False

def redistribute(binlist,rebinlist,binfo,weightmap):
    for bindex in range(len(rebinlist)):
        for poindex in range(len(rebinlist[bindex])):
            ## Finding the index of the bin with the closest centroid
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
    return numerator/np.sqrt(denominator)

            