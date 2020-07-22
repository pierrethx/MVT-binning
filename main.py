import numpy as np
import matplotlib.pyplot as plt 
import tkinter
from tkinter.filedialog import askopenfilename
from astropy.io import fits
import functions,bin_accretion,wvt_iteration
import scipy.spatial as sp
import time

#from astropy.utils.data import download_file

sourcedir,signal,var=bin_accretion.initialize(enternew=True)
#objname="J024815-081723"
objname="testdata"
subfolder=""
 ##empty for directly into sourcedir

start=time.time()

target=15 #Target Signal To noise
wvt,init_generators=bin_accretion.cc_accretion(signal,var,target)
init_scalelengths=np.full(len(init_generators),1)
epsilon=10
wvt=wvt_iteration.iteration_func(target,signal,var,init_generators,init_scalelengths,epsilon,displaywvt=False)

print("elapsed time: "+str(time.time()-start))

hdu = fits.PrimaryHDU(np.flipud(wvt))
hdul = fits.HDUList([hdu])
hdul.writeto(sourcedir+subfolder+"/iterated_"+objname+"_target_"+str(target)+".fits",overwrite=True)
fig,ax=plt.subplots()
image=ax.imshow(wvt,cmap="cubehelix")
fig.colorbar(image)
plt.show()
