import numpy as np
import matplotlib.pyplot as plt 
import tkinter
from tkinter.filedialog import askopenfilename
from matplotlib.widgets import Cursor,Slider
from astropy.io import fits
from astropy import wcs
import functions,bin_accretion,wvt_iteration
from scipy.ndimage.measurements import center_of_mass
from scipy import stats
from scipy.optimize import curve_fit
import time

def alignwcs(wcsx,angles):
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
    rdirections=np.array([np.cos(angles),np.sin(angles)]).T
    enddirs=wcs.WCS(wcsx).wcs_pix2world(rdirections,0)
    wdirs=np.array([ enddirs[i]-wcs.WCS(wcsx).wcs_pix2world([[0,0]],0)[0] for i in range(len(enddirs))]).T
    wangs=np.angle(wdirs[0]+1j*wdirs[1])
    directions=np.array([np.sin(wangs),np.cos(wangs)]).T
    return directions,tags

def getcenter(signal):
    siggy=np.ma.masked_where(signal<0,signal)
    center=center_of_mass(siggy) #as usual is (y,x)
    scale=0.5*np.sqrt(len(signal)**2+len(signal[0])**2)
    mask=np.copy(siggy)
    for yy in range(len(mask)):
        for xx in range(len(mask[yy])):
            if np.sqrt((yy-center[0])**2+(xx-center[1])**2)>scale/3:
                mask[yy][xx]=0
    center=center_of_mass(mask)
    return center

def radmethod(signal,var,wcsx,show=False,output=None):
    
    center=getcenter(signal)
    numdirections=56
    numsteps=200
    argus=np.linspace(0,2*np.pi,numdirections+1)[:-1]
    directions,tags=alignwcs(wcsx,argus)
    scale=0.5*np.sqrt(len(signal)**2+len(signal[0])**2)
    steps=np.linspace(0,1,numsteps)*scale

    histvalues=np.full((numdirections,numsteps),np.nan)
    varvalues=np.full((numdirections,numsteps),np.nan)

    edgecs=np.full(numdirections,np.nan)
    edgecs2=np.full(numdirections,np.nan)

    

    for dire in range(numdirections):
        for st in range(numsteps):
            rawpoint=directions[dire]*steps[st]

            #print(str(rawpoint)+" "+str(center))

            yy=int(rawpoint[0]+0.5+center[0])
            xx=int(rawpoint[1]+0.5+center[1])
            if yy>=len(signal) or yy<0 or xx>=len(signal[0]) or xx<0:
                pass
            else:
                histvalues[dire][st]=signal[yy][xx]
                varvalues[dire][st]=np.abs(var[yy][xx])
            if np.isnan(edgecs[dire]) and (not np.isnan(histvalues[dire][st])) and (histvalues[dire][st]<=0):
                edgecs[dire]=steps[st]
            
            
        print(str(dire)+" angle out of "+str(numdirections))
        #np.ma.masked_where(vibearray==0,vibearray)

    print(edgecs)
    

    x=np.ma.masked_where(np.isnan(histvalues[0]),steps)
    y=np.ma.masked_where(np.isnan(histvalues[0]),histvalues[0])
    yv=np.ma.masked_where(np.isnan(histvalues[0]),np.sqrt(varvalues[0]))
    if show:
        fig,ax = plt.subplots()

    def smartav(limst):
        neg=np.count_nonzero(np.isnan(limst))
        return (np.nansum(limst))/(len(limst)-neg)

    averages=np.array(list(map(lambda x:smartav(x),histvalues.T)))
    xa=np.ma.masked_where(np.isnan(averages),steps)
    ya=np.ma.masked_where(np.isnan(averages),averages)
    
    combinedstd=np.array(list(map(lambda x:smartav(x),varvalues.T)))
    yav=np.ma.masked_where(np.isnan(averages),np.sqrt(combinedstd))
    #lineaa,=ax.plot(xa, ya+yav, color="darkcyan",alpha=0.5,linestyle="dashed")
    #lineab,=ax.plot(xa, ya-yav, color="darkcyan",alpha=0.5,linestyle="dashed")
    if show:
        linea,=ax.plot(xa, ya, color="black",alpha=0.5,linestyle="dashed")
        line1,=ax.plot(x, y, color="red")
        line1a,=ax.plot(x, y+yv, color="turquoise")
        line1b,=ax.plot(x, y-yv, color="turquoise")
        ax.set_xlim(0,np.nanmax(xa))
        ax.set_ylim(0,np.nanmax(histvalues)*1.1)
        ax.set_title("Radial profile in direction "+tags[0])
        dirax= plt.axes([0.6, 0.55, 0.25, 0.25])
        dirax.imshow(signal,cmap="cubehelix")
        dirax.set_xlim(0,len(signal[0]))
        dirax.get_xaxis().set_ticks([])
        dirax.set_ylim(0,len(signal))
        dirax.get_yaxis().set_ticks([])
        ang,=dirax.plot(center[1]+[0,scale*directions[0][1]],center[0]+[0,scale*directions[0][0]],color="blue")

    yedge=np.linspace(0,np.nanmax(averages)*1.1,10)
    edge=False
    for val in range(len(steps)):
        edgex=steps[val]
        if np.isnan(edgex):
            break
        edgey=(histvalues[0]-np.sqrt(varvalues[0]))[val]
        if edgey<=0:
            edge=True
            break
    if show:
        if edge:
            ledge,=ax.plot(0*yedge+edgex,yedge,color="green",linestyle="dashed")
        else:
            ledge,=ax.plot(0*yedge+edgex,yedge,color="red",linestyle="dashed")

        def update(val):
            amp = samp.val
            n=int(amp)%numdirections

            line1.set_xdata(np.ma.masked_where(np.isnan(histvalues[n]),steps))
            line1.set_ydata(np.ma.masked_where(np.isnan(histvalues[n]),histvalues[n]))

            line1a.set_xdata(np.ma.masked_where(np.isnan(histvalues[n]),steps))
            line1a.set_ydata(np.ma.masked_where(np.isnan(histvalues[n]),histvalues[n]+np.sqrt(varvalues[n])))

            line1b.set_xdata(np.ma.masked_where(np.isnan(histvalues[n]),steps))
            line1b.set_ydata(np.ma.masked_where(np.isnan(histvalues[n]),histvalues[n]-np.sqrt(varvalues[n])))

            edge=False
            for val in range(len(steps)):
                edgex=steps[val]
                if np.isnan(edgex):
                    break
                edgey=(histvalues[n]-np.sqrt(varvalues[n]))[val]
                if edgey<=0:
                    edge=True
                    break
            if edge:
                ledge.set_xdata(0*yedge+edgex)
                ledge.set_color("green")
            else:
                ledge.set_xdata(0*yedge+edgex)
                ledge.set_color("red")

            ax.set_title("Radial profile in direction "+tags[n])
            ang.set_xdata(center[1]+[0,scale*directions[n][1]])
            ang.set_ydata(center[0]+[0,scale*directions[n][0]])
            #plt.pause(0.1)
            fig.canvas.draw()

        ampax = plt.axes([0.2, 0.02, 0.65, 0.03], facecolor="beige")
        samp = Slider(ampax, 'Amp', 0, numdirections, valinit=0, valstep=1,orientation="horizontal")
        samp.on_changed(update)

        plt.show()


        fig,ax = plt.subplots()
        ax.set_title("Average radial profile")
        ax.set_xlabel("radial distance from weighted center (px)")
        ax.set_ylabel("counts")
        line1,=ax.plot(xa, ya, color="black")
        ax.set_xlim(0,np.nanmax(xa))
        ax.set_ylim(0,np.nanmax(averages)*1.1)
        line1a,=ax.plot(xa, ya+yav, color="darkcyan")
        line1b,=ax.plot(xa, ya-yav, color="darkcyan")

    def varav(limst,vlimst):
        return (np.nansum(limst/vlimst))/np.nansum(1/vlimst)
    
    #avrad=varav(edgecs,edgecs2)
    #vrad=varav(edgecs**2,edgecs2)-avrad**2
    avrad=smartav(edgecs)
    vrad=smartav(edgecs**2)-avrad**2
    return avrad,vrad

    

if __name__=="__main__":
    try:
        wcsx,signal,var,source,objname=bin_accretion.initialize(enternew=True)
        repeat=True
    except:
        repeat=False

    while repeat:
        output2=[]
        
        avrad,vrad=radmethod(signal,var,wcsx,show=True)
        print(avrad,vrad)
        try:
            wcsx,signal,var,source,objname=bin_accretion.initialize(enternew=True)
            repeat=True
        except:
            repeat=False

    print("Bye Bye!")