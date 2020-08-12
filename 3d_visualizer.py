from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
import tkinter
from tkinter.filedialog import askopenfilename
from astropy.io import fits
import bin_accretion


wcsx2,signal2,var2,sourcedir2,objname2=bin_accretion.initialize(enternew=True)
ston=signal2/np.sqrt(var2)
ston=np.ma.masked_where(np.isnan(ston),ston)

x = np.arange(0, len(ston[0]), 1)
y = np.arange(0, len(ston),1)
X,Y=np.meshgrid(x,y)

fig = plt.figure()
ax = plt.axes(projection='3d')
#ax.contour3D(X, Y, signal, 5000, cmap='cubehelix')
ax.plot_surface(X, Y, ston, cmap='cubehelix')

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('brightness')
plt.show()