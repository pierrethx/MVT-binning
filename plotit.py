import numpy as np
import matplotlib.pyplot as plt 
import tkinter
from tkinter.filedialog import askopenfilename
from astropy.io import fits
import functions,bin_accretion,wvt_iteration
import scipy.spatial as sp

sourcedir,signal,var=bin_accretion.initialize(enternew=True)
fig,ax=plt.subplots(1,2)
a=ax[0].imshow(signal,cmap="cubehelix")
ax[0].set_title("signal")
b=ax[1].imshow(var,cmap="cubehelix")
ax[1].set_title("variance")
plt.show()
plt.figure()
plt.imshow(signal/np.sqrt(var),cmap="cubehelix")
plt.show()