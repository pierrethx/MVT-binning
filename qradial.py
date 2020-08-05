import numpy as np
import matplotlib.pyplot as plt 
import tkinter
from tkinter.filedialog import askopenfilename
from matplotlib.widgets import Cursor,Slider
from astropy.io import fits
import functions,bin_accretion,wvt_iteration,main
from scipy.ndimage.measurements import center_of_mass
import time

def alignwcs(wcs,angles):
    tags=[]
    for ang in (angles*180/np.pi):
        if ang<45:
            tags.append(str(round(ang,1))+"º North of East")
        elif ang<90:
            tags.append(str(round(90-ang,1))+"º East of North")
        elif ang<135:
            tags.append(str(round(ang-90,1))+"º West of North")
        elif ang<180:
            tags.append(str(round(180-ang,1))+"º North of West")
        elif ang<225:
            tags.append(str(round(ang-180,1))+"º South of West")
        elif ang<270:
            tags.append(str(round(270-ang,1))+"º West of South")
        elif ang<315:
            tags.append(str(round(ang-270,1))+"º East of South")
        else:
            tags.append(str(round(360-ang,1))+"º South of East")
    taglen=len(tags)
    if taglen%2==0:
        tags[0]="True East"
        tags[int(taglen/2)]="True West"
        if taglen%4==0:
            tags[int(taglen/4)]="True North"
            tags[int(3*taglen/4)]="True South"
            if taglen%8==0:
                tags[int(taglen/8)]="True Northeast"
                tags[int(3*taglen/8)]="True Northwest"
                tags[int(5*taglen/8)]="True Southwest"
                tags[int(7*taglen/8)]="True Southeast"
    rdirections=np.array([np.cos(argus),np.sin(argus)]).T
    enddirs=wcs.wcs_pix2world(rdirections,0)
    wdirs=np.array([ enddirs[i]-wcs.wcs_pix2world([[0,0]],0)[0] for i in range(len(enddirs))]).T
    wangs=np.angle(wdirs[0]+1j*wdirs[1])
    directions=np.array([np.sin(wangs),np.cos(wangs)]).T
    return directions,tags


sourcelist=[]
wscxlist=[]
siglist=[]
varlist=[]
objlist=[]
contqueue=True
while contqueue:
    try:
        wcsx,signal,var,sourcedir,objname=bin_accretion.initialize(enternew=True)
        sourcelist.append(sourcedir.split("/")[-1])
        objlist.append(objname)
        wscxlist.append(wcsx)
        siglist.append(signal)
        varlist.append(var)
    except:
        contqueue=False
print("Files loaded!")
for i in range(len(sourcelist)):
    signal=siglist[i]
    var=varlist[i]
    wcsx=wscxlist[i]

    center=center_of_mass(signal) #as usual is (y,x)
    numdirections=56
    numsteps=50
    argus=np.linspace(0,2*np.pi,numdirections+1)[:-1]
    directions,tags=alignwcs(wcsx,argus)
    scale=0.5*np.sqrt(len(signal)**2+len(signal[0])**2)
    steps=np.linspace(0,1,numsteps)*scale

    histvalues=np.full((numdirections,numsteps),np.nan)
    varvalues=np.full((numdirections,numsteps),np.nan)
    #print(len(histvalues[0]))
    for dire in range(numdirections):
        for st in range(numsteps):
            rawpoint=directions[dire]*steps[st]
            yy=int(rawpoint[0]+0.5+center[0])
            xx=int(rawpoint[1]+0.5+center[1])
            if yy>=len(signal) or yy<0 or xx>=len(signal[0]) or xx<0:
                pass
            else:
                histvalues[dire][st]=signal[yy][xx]
                varvalues[dire][st]=np.abs(var[yy][xx])
        #np.ma.masked_where(vibearray==0,vibearray)

    x=np.ma.masked_where(np.isnan(histvalues[0]),steps)
    y=np.ma.masked_where(np.isnan(histvalues[0]),histvalues[0])
    yv=np.ma.masked_where(np.isnan(histvalues[0]),np.sqrt(varvalues[0]))

    def smartav(limst):
        neg=np.count_nonzero(np.isnan(limst))
        return (np.nansum(limst))/(len(limst)-neg)

    averages=np.array(list(map(lambda x:smartav(x),histvalues.T)))
    xa=np.ma.masked_where(np.isnan(averages),steps)
    ya=np.ma.masked_where(np.isnan(averages),averages)
    combinedstd=np.array(list(map(lambda x:smartav(x),varvalues.T)))
    yav=np.ma.masked_where(np.isnan(averages),np.sqrt(combinedstd))



    edge=False
    for val in range(len(ya)):
        edgex=xa[val]
        edgey=(ya-yav)[val]
        if edgey<=0:
            edge=True
            break
    print(sourcelist[i]+"/"+objlist[i],end=": ")
    if edge:
        print(edgex)
    else:
        print("no edge")

print("Bye Bye!")