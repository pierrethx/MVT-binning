import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.widgets import Cursor,Slider
import matplotlib as mpl
import tkinter
from tkinter.filedialog import askopenfilename
from astropy.io import fits
from astropy import wcs
from astropy.visualization.wcsaxes import WCSAxes
import functions,bin_accretion,wvt_iteration
import scipy.spatial as sp

wcsx2,signal2,var2,sourcedir2,objname2=bin_accretion.initialize(enternew=True)
ston=signal2/var2

fig=plt.figure()

g=plt.imshow(ston,cmap="cubehelix")

fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=5, vmax=10), cmap=mpl.cm.cubehelix))
plt.show()

