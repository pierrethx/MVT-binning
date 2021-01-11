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

def circularb(x,A,r,b):
    return A*(1+(x/r)**2)**(0.5-3*b)

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

    

    for dire in range(numdirections):
        for st in range(numsteps):
            rawpoint=directions[dire]*steps[st]

            print(str(rawpoint)+" "+str(center))

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

    edge=False
    yedge=np.linspace(0,np.nanmax(averages)*1.1,10)
    for val in range(len(ya)):
        edgex=xa[val]
        edgey=(ya-yav)[val]
        if edgey<=0:
            edge=True
            break
    try:
        med=np.nanmedian(edgecs)
    except:
        med=np.nan
    try:
        avl=np.nanmean(edgecs)
    except:
        avl=np.nan
    if show:
        if edge:
            ax.plot(0*yedge+edgex,yedge,color="green",linestyle="dashed")
            ax.plot([med,med],[np.nanmin(ya),np.nanmax(ya)],color="gold",linestyle="dashed")
            ax.annotate('potential edge', xy=(edgex, edgey), xytext=(edgex-np.max(xa)*0.02, np.nanmax(averages)*.2),horizontalalignment='right')
            print("edge "+str(edgex))
            try:
                slope,intercept,rv,pv,sterr=stats.linregress(xa[val+1:],ya[val+1:])
                print("slope "+str(slope))
                print("intercept "+str(intercept))
                print("rv "+str(rv))
                print("pv "+str(pv))
                print("sterr "+str(sterr))
            except:
                pass
        else:
            ax.plot(0*yedge+edgex,yedge,color="red",linestyle="dashed")
            ax.annotate('no edge detected', xy=(edgex, edgey), xytext=(edgex-np.max(xa)*0.02, np.nanmax(averages)*.2),horizontalalignment='right')
            print("no edge")

        plt.show()

        xl=functions.lerp(xa)
        yl=functions.lerp(ya)
        div1=functions.numdif(xl,yl)
        fig,ax = plt.subplots()
        ax.set_title("First derivative of radial profile")
        ax.set_xlabel("radial distance from weighted center (px)")
        ax.set_ylabel("counts")

        line1,=ax.plot(xl, div1, color="black")

        lines1,=ax.plot(xl, functions.smoother(div1,1), color="red")
        lines2,=ax.plot(xl, functions.smoother(div1,2), color="orange")
        lines3,=ax.plot(xl, functions.smoother(div1,3), color="yellow")
        lines4,=ax.plot(xl, functions.smoother(div1,4), color="green")
        lines5,=ax.plot(xl, functions.smoother(div1,5), color="blue")
        lines6,=ax.plot(xl, functions.smoother(div1,6), color="purple")
        ax.set_xlim(0,np.nanmax(xa))
        ax.set_ylim(np.nanmin(div1)*1.1,np.nanmax(div1)*1.1)
        

        if edge:
            ax.plot([edgex,edgex],[np.nanmin(div1),np.nanmax(div1)],color="green",linestyle="dashed")
            ax.annotate('potential edge', xy=(edgex, edgey), xytext=(edgex-np.max(xa)*0.02, np.nanmax(averages)*.2),horizontalalignment='right')
        else:
            ax.plot([edgex,edgex],[np.nanmin(div1),np.nanmax(div1)],color="red",linestyle="dashed")
            ax.annotate('no edge detected', xy=(edgex, edgey), xytext=(edgex-np.max(xa)*0.02, np.nanmax(averages)*.2),horizontalalignment='right')


        plt.show()

    I_0=np.nanmax(ya)
    smallest=1
    for l in range(len(ya)):
        if ya[l]>0.6065*I_0:
            if ya[l]<ya[smallest]:
                smallest=l
        else:
            break

    mini=0
    inflection=np.argmin(functions.smoother(functions.numdif(xa,ya),6)[mini:smallest+1])
    if inflection==mini:
        mini+=1
        inflection=np.argmin(functions.smoother(functions.numdif(xa,ya),6)[mini:smallest+1])+mini
    if inflection is list:
        inflection=inflection[-1]
    
    beta=0.07/(ya[inflection]/I_0 - 0.6065)
    if beta<0.168:
        beta=0.168
    r_c=xa[inflection]*np.sqrt(6*beta)

    print(inflection)

    p_0=[I_0,r_c,beta]

    print(p_0)
    try:
        if np.isnan(med):
            raise NameError("")
        val=np.argmin((xa-med)**2)
    except:
        val=len(xa)
    if val==0:
        popt=p_0
    else:
        try:
            print(ya[int(val/10):val])
            popt,pcov=curve_fit(circularb,xa[int(val/10):val],ya[int(val/10):val],p0=p_0,sigma=yav[int(val/10):val],absolute_sigma=True)
            perr = np.sqrt(np.diag(pcov))
        except:
            popt,pcov=p_0,np.sqrt(p_0)
            perr = np.sqrt(np.diag(pcov))

    try:
        for p in popt:
            output.append(p)
    except:
        pass

    print(output)
    yp=circularb(xa,popt[0],popt[1],popt[2])
    if show:
        fig,ax = plt.subplots()
        #ax.set_title("Average radial profile")
        #ax.set_xlabel("radial distance from weighted center (px)")
        #ax.set_ylabel("counts")
        line1,=ax.plot(xa, ya, color="black")
        ax.set_xlim(0,np.nanmax(xa))
        ax.set_ylim(0,np.nanmax(averages)*1.1)
        #linerf,=ax.plot(xa, circularb(xa,I_0,r_c,beta), color="red")

        
        lineff,=ax.plot(xa, yp, color="blue")

    for j in range(len(yp)):
        if yp[j]/ya[j]>2:
            break
    print(xa[j])
    
    
    print(med)
    if show:
        ax.plot([xa[j],xa[j]],[np.nanmin(ya),np.nanmax(ya)],color="green",linestyle="dashed")
        ax.plot([med,med],[np.nanmin(ya),np.nanmax(ya)],color="orange",linestyle="dashed")
        plt.show()

    print(avl)
    #edgex=0 edge from average profile
    #xa[j]=fit edge from profile
    #med=median of 0 edges from directions
    #avl=average of 0 edges from directions
    return edgex,xa[j],med,avl

if __name__=="__main__":
    try:
        wcsx,signal,var,source,objname=bin_accretion.initialize(enternew=True)
        repeat=True
    except:
        repeat=False

    while repeat:
        output2=[]
        edgex,edge2,edgec,edgea=radmethod(signal,var,wcsx,show=True,output=output2)
        print(output2)
        try:
            wcsx,signal,var,source,objname=bin_accretion.initialize(enternew=True)
            repeat=True
        except:
            repeat=False

    print("Bye Bye!")