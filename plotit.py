import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.widgets import Cursor,Slider
import matplotlib as mpl
import tkinter
from tkinter.filedialog import askopenfilename
from astropy.io import fits
from astropy import wcs
from astropy.visualization.wcsaxes import WCSAxes
import functions,bin_accretion,wvt_iteration,radial_profile
import scipy.spatial as sp

wcsx,signal,var,sourcedir,objname=bin_accretion.initialize(enternew=True)
## objname would be testdata_edge##_wit_###.fits
## sourcedir would be stuff/stuff/128x128_peak100/bg##/target##

source="/".join(sourcedir.split("/")[:-2])
name="_".join(objname.split("_")[:-1])
sigfits=name+"_sig.fits"
varfits=name+"_var.fits"
stonfits="zston_"+name+".fits"

bg=[20,40,60,80,100,120]
target=[5,10,15,20]

fig=plt.figure()
scale=0.8
ax=plt.axes([(1-scale)/2,(1-scale)/2,scale,scale])
ax.set_title(name)
ax.set_xlim(0,len(bg))
ax.set_ylim(0,len(target))
ax.set_xticks([i+0.5 for i in range(len(bg))])
ax.set_xticklabels([bg[i] for i in range(len(bg))])
ax.set_yticks([i+0.5 for i in range(len(target))])
ax.set_yticklabels([target[i] for i in range(len(target))])
ax.set_xlabel(r"background noise level (ct$•$pix$^{-1}$)")
ax.set_ylabel(r"target signal-to-noise (pix$^{-1}$)")
lis=[]
#dirax= plt.axes([0.6, 0.55, 0.25, 0.25])
for b in range(len(bg)):
    lis.append([])
    for t in range(len(target)):
        signalname=source+"/bg"+str(bg[b])+"/target"+str(target[t])+"/"+stonfits
        #varname=source+"/bg"+bg[b]+"/target"+target[t]+"/"+varfits

        with fits.open(signalname) as hdul:
            signal=np.flipud(hdul[0].data)
            wcsx=wcs.WCS(hdul[0].header)
        #with fits.open(varname) as hdul:
        #    var=np.flipud(hdul[0].data)
        '''
        edgex,edge2,edgec=radial_profile.radmethod(signal,signal,wcsx,show=False)
        '''
        bum=plt.axes([(1-scale)/2+scale*b/len(bg),(1-scale)/2+scale*t/len(target),scale/len(bg),scale/len(target)])
        lis[b].append(bum)
        #lis[b][t].annotate(str(bg[b])+" at "+str(target[t]), xy=(0, 0), xytext=(0,0))
        bum.imshow(signal,cmap="cubehelix",vmin=0,vmax=30)
        '''
        t=np.linspace(0,2*np.pi,100)
        x=edgec*np.cos(t)+len(signal[0])/2
        y=edgec*np.sin(t)+len(signal)/2
        bum.plot(x,y,color='red')
        '''
        bum.set_xticks([])
        bum.set_yticks([])
cax=plt.axes([(9+7*scale)/16,(1-scale)/2,(1-scale)/8,scale])

clbr=plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=0, vmax=30), cmap=mpl.cm.get_cmap('cubehelix')),cax=cax)
plt.show()
