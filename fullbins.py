import numpy as np
import matplotlib.pyplot as plt
import tkinter
from tkinter.filedialog import askopenfilename
from astropy.io import fits
import bin_accretion
from astropy import wcs
from scipy.ndimage.measurements import center_of_mass

import scipy.spatial as sp 

scaleup=1

bloop=True
alist=[]
slist=[]
rslist=[]
vlist=[]
rvlist=[]
while (bloop):
    try:
        root=tkinter.Tk()
        root.withdraw()
        placeholder= askopenfilename(message="Select assign",multiple=False)
        root.update()
        root.destroy()
        with fits.open(placeholder) as hdul:
            alist.append(np.flipud(hdul[0].data).astype(int))
            wcsx=wcs.WCS(hdul[0].header)
        if placeholder.split("/")[-1][0]=='z':
            with fits.open("/".join(placeholder.split("/")[:-1])+"/"+"_".join(placeholder.split("_")[-3:-1])+"_wit_sig.fits") as hdul:
                slist.append(np.flipud(hdul[0].data))
                wcsx=wcs.WCS(hdul[0].header)
            with fits.open("/".join(placeholder.split("/")[:-1])+"/"+"_".join(placeholder.split("_")[-3:-1])+"_wit_var.fits") as hdul:
                vlist.append(np.flipud(hdul[0].data))
                wcsx=wcs.WCS(hdul[0].header)
            with fits.open("/".join(placeholder.split("/")[:-2])+"/unbinned/"+"_".join(placeholder.split("_")[-3:-1])+"_signal.fits") as hdul:
                rslist.append(np.flipud(hdul[0].data))
                wcsx=wcs.WCS(hdul[0].header)
            with fits.open("/".join(placeholder.split("/")[:-2])+"/unbinned/"+"_".join(placeholder.split("_")[-3:-1])+"_var.fits") as hdul:
                rvlist.append(np.flipud(hdul[0].data))
                wcsx=wcs.WCS(hdul[0].header)
        else:
            with fits.open("/".join(placeholder.split("/")[:-1])+"/"+"_".join(placeholder.split("_")[:-1])+"_wit_sig.fits") as hdul:
                slist.append(np.flipud(hdul[0].data))
                wcsx=wcs.WCS(hdul[0].header)
            with fits.open("/".join(placeholder.split("/")[:-1])+"/"+"_".join(placeholder.split("_")[:-1])+"_wit_var.fits") as hdul:
                vlist.append(np.flipud(hdul[0].data))
                wcsx=wcs.WCS(hdul[0].header)
            with fits.open("/".join(placeholder.split("/")[:-2])+"/unbinned/"+"_".join(placeholder.split("_")[:-1])+"_signal.fits") as hdul:
                rslist.append(np.flipud(hdul[0].data))
                wcsx=wcs.WCS(hdul[0].header)
            with fits.open("/".join(placeholder.split("/")[:-2])+"/unbinned/"+"_".join(placeholder.split("_")[:-1])+"_var.fits") as hdul:
                rvlist.append(np.flipud(hdul[0].data))
                wcsx=wcs.WCS(hdul[0].header)
        
    except:
        bloop=False


for a in range(len(alist)):
    assign=alist[a]

    bins=(np.nanmax(assign))+1

    binlist=[[] for i in range(bins)]
    for y in range(len(assign)):
        for x in range(len(assign[y])):
            binlist[assign[y][x]-1].append([y,x])
    binlist[-1]=[]
    sitem=np.copy(slist[a])
    vitem=np.copy(vlist[a])
    for binn in binlist:
        for point in binn:    
            sitem[point[0]][point[1]]*=len(binn)
            vitem[point[0]][point[1]]*=len(binn)
    fig,ax=plt.subplots(3,3)
    ax[0][2].imshow(rslist[a])
    ax[0][2].set_title("unbinned signal")
    ax[0][0].imshow(slist[a])
    ax[0][0].set_title("signal")
    ax[0][1].imshow(sitem)
    ax[0][1].set_title("pooled signal")
    ax[1][2].imshow(rvlist[a])
    ax[1][2].set_title("raw variance")
    ax[1][0].imshow(vlist[a])
    ax[1][0].set_title("variance")
    ax[1][1].imshow(vitem)
    ax[1][1].set_title("pooled variance")
    ax[2][1].imshow(sitem/np.sqrt(vitem))
    ax[2][1].set_title("StoN")
    plt.show()
    