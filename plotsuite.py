import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.widgets import Cursor,Slider
import matplotlib as mpl
import tkinter
from tkinter.filedialog import askopenfilename
from astropy.io import fits
from astropy import wcs
from astropy.visualization.wcsaxes import WCSAxes
import functions,bin_accretion,wvt_iteration,radial_profile,main
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

wcsx,signal,var,sourcedir,objname=bin_accretion.minitialize()
## objname would be testdata_edge##_wit_###.fits
## sourcedir would be stuff/stuff/128x128_peak100/bg##/target##


for w in range(len(wcsx)):
    obj="_".join(objname[w].split("_")[:-1])
    SNR=signal[w]/np.sqrt(var[w])
    SNR[SNR==0]=np.NaN
    SNRf=np.flipud(SNR)

    levs=["target3","target5","target10"]
    max=15
    min=np.nanmin(SNR)

    fig=plt.figure()
    ax=plt.axes(projection=wcs.WCS(wcsx[w]))
    ax.set_xlabel('Right Ascension')
    ax.set_ylabel('Declination')
    g=ax.imshow(SNRf,cmap="coolwarm",norm=mpl.colors.TwoSlopeNorm(vmin=min, vcenter=0., vmax=max))
    plt.colorbar(g)
    plt.grid(color='black',ls='solid')
    fig.savefig(sourcedir[w]+"/"+obj+"_SNR_unbinned.png",bbox_inches="tight")

    for n in range(len(levs)):
        with fits.open(sourcedir[w]+"/"+levs[n]+"/zston_"+obj+".fits") as hdul:
            sigt=hdul[0].data
            sigt[sigt==0]=np.NaN
        fig=plt.figure()
        ax=plt.axes(projection=wcs.WCS(wcsx[w]))
        ax.set_xlabel('Right Ascension')
        ax.set_ylabel('Declination')
        g=ax.imshow(sigt,cmap="coolwarm",norm=mpl.colors.TwoSlopeNorm(vmin=min, vcenter=0., vmax=max))
        plt.colorbar(g)
        plt.grid(color='black',ls='solid')
        fig.savefig(sourcedir[w]+"/"+obj+"_SNR_"+levs[n]+".png",bbox_inches="tight")
