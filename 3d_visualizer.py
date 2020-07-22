from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
import tkinter
from tkinter.filedialog import askopenfilename
from astropy.io import fits


tkinter.Tk().withdraw()
signalname = askopenfilename()
with fits.open(signalname) as hdul:
    signal=np.abs(hdul[0].data)
x = np.arange(0, len(signal[0]), 1)
y = np.arange(0, len(signal),1)
X,Y=np.meshgrid(x,y)

fig = plt.figure()
ax = plt.axes(projection='3d')
#ax.contour3D(X, Y, signal, 5000, cmap='cubehelix')
ax.plot_surface(X, Y, signal, cmap='cubehelix')

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('brightness')
plt.show()