import numpy as np
import matplotlib.pyplot as plt 
import tkinter
from tkinter.filedialog import askopenfilename
from matplotlib.widgets import Cursor,Slider,RadioButtons
import matplotlib as mpl
from matplotlib import cm
from astropy.io import fits
import functions,bin_accretion,wvt_iteration,main,radial_profile
from scipy.ndimage.measurements import center_of_mass
from scipy.optimize import curve_fit
import time

def endnumber(stri):
    beg=-1
    while stri[beg].isdigit():
        beg-=1
    if beg==-1:
        return 0.5*(np.sqrt((128**2+128**2)))
    else:
        print(int(stri[beg+1:]))
        return int(stri[beg+1:])

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
    rdirections=np.array([np.cos(angles),np.sin(angles)]).T
    enddirs=wcs.wcs_pix2world(rdirections,0)
    wdirs=np.array([ enddirs[i]-wcs.wcs_pix2world([[0,0]],0)[0] for i in range(len(enddirs))]).T
    wangs=np.angle(wdirs[0]+1j*wdirs[1])
    directions=np.array([np.sin(wangs),np.cos(wangs)]).T
    return directions,tags

def qradial(wcsxlist,siglist,varlist,sourcelist,objlist):
    datapoints=[]
    for i in range(len(sourcelist)):
        signal=siglist[i]
        var=varlist[i]
        wcsx=wcsxlist[i]

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
        '''
        y=np.ma.masked_where(np.isnan(histvalues[0]),histvalues[0])
        yv=np.ma.masked_where(np.isnan(histvalues[0]),np.sqrt(varvalues[0]))
        '''
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

        I_0=(ya[0]+np.nanmax(ya))/2
        smallest=1
        for l in range(len(ya)):
            if ya[l]>0.6065*I_0:
                if ya[l]<ya[smallest]:
                    smallest=l
            else:
                break

        inflection=np.argmin(functions.smoother(functions.numdif(xa,ya),6)[1:smallest+1])
        if inflection is list:
            inflection=inflection[-1]
        
        
        
        beta=0.07/(ya[inflection]/I_0 - 0.6065)
        r_c=xa[inflection]*np.sqrt(6*beta)

        print(xa[inflection])

        p_0=[I_0,r_c,beta]

        print(p_0)

        val=int(val*.9)

        popt,pcov=curve_fit(radial_profile.circularb,xa[:val+1],ya[:val+1],p0=p_0,sigma=yav[:val+1],absolute_sigma=True)
        perr = np.sqrt(np.diag(pcov))

        print(popt)
        yp=circularb(xa,popt[0],popt[1],popt[2])
        for j in range(len(yp)):
            if yp[j]/ya[j]>1.6:
                break
        print(xa[j])

        datapoint=[endnumber(sourcelist[i]),endnumber(objlist[i].split("_")[1]),edgex,endnumber(sourcelist[i].split("/")[-2]),xa[j]]
        datapoints.append(datapoint)

    data=np.array(datapoints)

    fig, ax = plt.subplots()
    cmap = cm.get_cmap('inferno')
    linearray=[]
    datam=data[:,0]
    ax.set_xlim(0.8*np.nanmin(datam),1.2*np.nanmax(datam))
    ax.set_ylim(0.8*np.nanmin(datam),1.2*np.nanmax(datam))
    
    xver=np.linspace(0,1000,2)
    edgelim,=ax.plot(xver,xver*0-10,color="red",linestyle="dashed")

    for m in datam:
        l, = ax.plot(m, m,markersize=8,marker="o",color=cmap(m/np.nanmax(datam)))
        linearray.append(l)
    plt.subplots_adjust(left=0.3)
    clbr=fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=np.nanmin(datam), vmax=np.nanmax(datam)), cmap=cmap))

    axcolor = 'lightgoldenrodyellow'
    xax = plt.axes([0.05, 0.4, 0.15, 0.15], facecolor=axcolor)
    radio = RadioButtons(xax, ('target', 'seed edge', 'rough edge','background','found edge'))
    xax.set_title("x axis")

    def xfunc(label):
        axisdict = {'target': data[:,0], 'seed edge': data[:,1], 'rough edge': data[:,2], 'background': data[:,3],'found edge':data[:,4]}
        axisdata = axisdict[label]
        for m in range(len(axisdata)):
            linearray[m].set_xdata(axisdata[m])
        ax.set_xlim(0.8*np.nanmin(axisdata),1.2*np.nanmax(axisdata))
        plt.draw()
    radio.on_clicked(xfunc)

    yax = plt.axes([0.05, 0.7, 0.15, 0.15], facecolor=axcolor)
    radio2 = RadioButtons(yax, ('target', 'seed edge', 'rough edge','background','found edge'))
    yax.set_title("y axis")


    def yfunc(label):
        axisdict = {'target': data[:,0], 'seed edge': data[:,1], 'rough edge': data[:,2], 'background': data[:,3],'found edge':data[:,4]}
        axisdata = axisdict[label]

        if label=='seed edge' or label=='found edge' or label=='rough edge':
            edgelim.set_ydata(xver*0+0.5*(np.sqrt((128**2+128**2))))
        else:
            edgelim.set_ydata(xver*0-10)
        
        for m in range(len(axisdata)):
            linearray[m].set_ydata(axisdata[m])
        ax.set_ylim(0.8*np.nanmin(axisdata),1.2*np.nanmax(axisdata))
        plt.draw()
    radio2.on_clicked(yfunc)

    zax = plt.axes([0.05, 0.1, 0.15, 0.15], facecolor=axcolor)
    radio3 = RadioButtons(zax, ('target', 'seed edge', 'rough edge','background','found edge'))
    zax.set_title("colors")


    def zfunc(label):
        axisdict = {'target': data[:,0], 'seed edge': data[:,1], 'rough edge': data[:,2], 'background': data[:,3],'found edge':data[:,4]}
        axisdata = axisdict[label]
        for m in range(len(axisdata)):
            linearray[m].set_color(cmap(axisdata[m]/np.nanmax(axisdata)))
        
        #clbr.remove()
        clbr.update_normal(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=np.nanmin(axisdata), vmax=np.nanmax(axisdata)), cmap=cmap))
        plt.draw()
    radio3.on_clicked(zfunc)

    

    plt.show()

if __name__ == "__main__":
    sourcelist=[]
    wcsxlist=[]
    siglist=[]
    varlist=[]
    objlist=[]
    contqueue=True
    while contqueue:
        try:
            wcsx,signal,var,sourcedir,objname=bin_accretion.initialize(enternew=True)
            sourcelist.append(sourcedir)
            objlist.append(objname)
            wcsxlist.append(wcsx)
            siglist.append(signal)
            varlist.append(var)
        except:
            contqueue=False
    print("Files loaded!")
    qradial(wcsxlist,siglist,varlist,sourcelist,objlist)
    print("Bye Bye!")