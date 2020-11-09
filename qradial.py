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
import time,os

def endnumber(stri):
    beg=-1
    while stri[beg].isdigit():
        beg-=1
    if beg==-1:
        return 0
    else:
        print("length: "+(stri[beg+1:]))
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

def qradial(wcsxlist,siglist,varlist,sourcelist,objlist,outputlist=None):
    datapoints=[]
    for i in range(len(sourcelist)):
        signal=siglist[i]
        var=varlist[i]
        wcsx=wcsxlist[i]
        output=[]

        edgex,edge2,edgec,edgea=radial_profile.radmethod(signal,var,wcsx,show=False,output=output)
        try:
            outputlist.append(output)
        except:
            pass
        print(sourcelist[i])
        if "block" in objlist[i]:
            datapoint=[endnumber(sourcelist[i]),endnumber(objlist[i].split("_")[2]),edgex,endnumber(sourcelist[i].split("/")[-2]),edge2,edgec,edgea]
        else:
            datapoint=[endnumber(sourcelist[i]),endnumber(objlist[i].split("_")[1]),edgex,endnumber(sourcelist[i].split("/")[-2]),edge2,edgec,edgea]

        print(datapoint)
        """[binning,seededge,firstpassedge,exptime,foundedge,stonedgem,stonedgea]"""
        datapoints.append(datapoint)

    data=np.array(datapoints)
    '''
    fig, ax = plt.subplots()
    cmap = cm.get_cmap('inferno')
    linearray=[]
    datam=data[:,0]
    ax.set_xlim(0.8*np.nanmin(datam),1.2*np.nanmax(datam))
    ax.set_ylim(0.8*np.nanmin(datam),1.2*np.nanmax(datam))
    
    xver=np.linspace(0,1000,2)
    edgelim,=ax.plot(xver,xver*0-10,color="red",linestyle="dashed")

    axisdict = {'target': data[:,0], 'seed edge': data[:,1], 'rough edge': data[:,2], 'background': data[:,3],'found edge':data[:,4],"stonedge":data[:,5],"stonedgea":data[:,6]}

    for m in datam:
        l, = ax.plot(m, m,markersize=8,marker="o",color=cmap(m/np.nanmax(datam)))
        linearray.append(l)
    plt.subplots_adjust(left=0.3)
    clbr=fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=np.nanmin(datam), vmax=np.nanmax(datam)), cmap=cmap))

    axcolor = 'lightgoldenrodyellow'
    xax = plt.axes([0.05, 0.7, 0.15, 0.15], facecolor=axcolor)
    radio = RadioButtons(xax, tuple(axisdict.keys()))
    xax.set_title("x axis")

    def xfunc(label):
        axisdata = axisdict[label]
        for m in range(len(axisdata)):
            linearray[m].set_xdata(axisdata[m])
        ax.set_xlim(0.8*np.nanmin(axisdata),1.2*np.nanmax(axisdata))
        plt.draw()
    radio.on_clicked(xfunc)

    yax = plt.axes([0.05, 0.4, 0.15, 0.15], facecolor=axcolor)
    radio2 = RadioButtons(yax, tuple(axisdict.keys()))
    yax.set_title("y axis")


    def yfunc(label):
        axisdata = axisdict[label]

        if label=='seed edge' or label=='found edge' or label=='rough edge' or label=="stonedge":
            edgelim.set_ydata(xver*0+0.5*(np.sqrt((len(siglist[0])**2+len(siglist[0][0])**2))))
        else:
            edgelim.set_ydata(xver*0-10)
        
        for m in range(len(axisdata)):
            linearray[m].set_ydata(axisdata[m])
        ax.set_ylim(0.8*np.nanmin(axisdata),1.2*np.nanmax(axisdata))
        plt.draw()
    radio2.on_clicked(yfunc)

    zax = plt.axes([0.05, 0.1, 0.15, 0.15], facecolor=axcolor)
    radio3 = RadioButtons(zax, tuple(axisdict.keys()))
    zax.set_title("colors")


    def zfunc(label):
        
        axisdata = axisdict[label]
        for m in range(len(axisdata)):
            linearray[m].set_color(cmap(axisdata[m]/np.nanmax(axisdata)))
        
        #clbr.remove()
        clbr.update_normal(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=np.nanmin(axisdata), vmax=np.nanmax(axisdata)), cmap=cmap))
        plt.draw()
    radio3.on_clicked(zfunc)

    

    plt.show()
    '''

    return data

def tabulate(pathfile,headers,rows,rel=False):
    f=open(pathfile,"w+")
    f.write("\\begin{tab"+"ular}{|c|"+"|c"*len(headers)+"|}\n \\hline\n")
    f.write("& ".join(headers)+"\\\\ \n \\hline \n")
    
    
    cols=[endnumber(stri) for stri in headers]
    print(cols)
    shortdata=np.zeros(shape=(len(rows),len(cols)))
    for r in range(len(rows)):
        shortdata[r][0]=rows[r]
    for point in data:
        if rel:
            shortdata[rows.index(point[3])][cols.index(point[1])]=round(point[6]/point[1],3)
        else:
            shortdata[rows.index(point[3])][cols.index(point[1])]=round(point[6],3)
    for line in shortdata:
        f.write(" & ".join([str(n) for n in line])+"\\\\ \n")
    f.write("\\hline\n\\end{tab"+"ular}")

def tabulate2(pathfile,headers,rows,row2,rel=False):
    f=open(pathfile,"w+")
    f.write("\\begin{tab"+"ular}{|c|"+"|c"*len(headers)+"|}\n \\hline\n")
    f.write("& ".join(headers)+"\\\\ \n \\hline \n")
    
    
    cols=[endnumber(stri) for stri in headers]
    print(cols)
    shortdata=np.zeros(shape=(len(rows),len(cols)))
    namelist=[]
    for r in range(len(rows)):
        namelist.append(str(rows[r])+", "+str(row2[r]))
        shortdata[r][0]=-1
    print(namelist)
    for point in data:
        if rel:
            shortdata[namelist.index(str(point[3])+", "+str(int(point[0])))][cols.index(point[1])]=round(point[6]/point[1],3)
        else:
            shortdata[namelist.index(str(point[3])+", "+str(int(point[0])))][cols.index(point[1])]=round(point[6],3)
    for r in range(len(shortdata)):
        f.write(namelist[r]+" & "+" & ".join([str(n) for n in shortdata[r] if np.isnan(n) or n>=0])+"\\\\ \n")
    f.write("\\hline\n\\end{tab"+"ular}")

if __name__ == "__main__":
    wcsxlist,siglist,varlist,sourcelist,objlist=bin_accretion.minitialize()
    fname=main.getname(sourcelist[0].split("/")[-1]+".txt")
    print("Files loaded!")
    data=qradial(wcsxlist,siglist,varlist,sourcelist,objlist)
    """[binning,seededge,firstpassedge,exptime,foundedge,stonedgem,stonedgea]"""
    
    path="/".join(sourcelist[0].split("/")[:-2])+"/chart/"
    
    headers=["expt, bin","edge50","edge75","(none) edge100"]
    rows=[1e2,1e2,1e2,1e2,1e4,1e4,1e4,1e4,1e6,1e6,1e6,1e6]
    row2=[0,3,5,10,0,3,5,10,0,3,5,10]

    tabulate2(path+fname,headers,rows,row2,rel=False)
    tabulate2(path+"r-"+fname,headers,rows,row2,rel=True)
    '''
    for fil in range(len(data)):
        print(objlist[fil],end="  ")
        print("fitedge",end=" ")
        print(data[fil][4],end="  ")
        print("stonedge",end=" ")
        print(data[fil][-1])
        '''
    print("Bye Bye!")