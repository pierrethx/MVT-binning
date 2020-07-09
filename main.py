import numpy as np
import matplotlib.pyplot as plt 
import tkinter
from tkinter.filedialog import askopenfilename
from astropy.io import fits
import functions,bin_accretion,wvt_iteration
import scipy.spatial as sp

#from astropy.utils.data import download_file

sourcedir,signal,var=bin_accretion.initialize()
objname="J024815-081723"
target=3 #Target Signal To noise
wvt,init_generators, init_scalelengths=bin_accretion.cc_accretion(signal,var,target)
epsilon=10
wvt=wvt_iteration.iteration_func(target,signal,var,init_generators,init_scalelengths,epsilon)
hdu = fits.PrimaryHDU(wvt)
hdul = fits.HDUList([hdu])
hdul.writeto(sourcedir+"/iterated_"+objname+"_target_"+str(target)+".fits",overwrite=True)
fig,ax=plt.subplots()
image=ax.imshow(wvt,cmap="cubehelix")
fig.colorbar(image)
plt.show()
