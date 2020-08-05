import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.widgets import Cursor,Slider
import tkinter
from tkinter.filedialog import askopenfilename
from astropy.io import fits
from astropy import wcs
from astropy.visualization.wcsaxes import WCSAxes
import functions,bin_accretion,wvt_iteration
import scipy.spatial as sp
'''
sourcedir,wcsx,signal,var=bin_accretion.initialize(enternew=True)
fig,ax=plt.subplots(1,2)
a=ax[0].imshow(signal,cmap="cubehelix")
ax[0].set_title("signal")
b=ax[1].imshow(var,cmap="cubehelix")
ax[1].set_title("variance")
plt.show()

fig,ax=plt.subplots()
burp=signal/np.sqrt(var)
plt.imshow(burp,cmap="cubehelix")
def onclick(event):
    x1,y1=event.xdata,event.ydata
    print(burp[int(y1+.5)][int(x1+.5)])
cursor=Cursor(ax, horizOn=False,vertOn=False,color='red',linewidth=2.0)
fig.canvas.mpl_connect('button_press_event',onclick)
plt.show()


pixcrd=np.array([[0,0],[1000000,0]],dtype=np.float64)
pixcrd2=wcsx.wcs_pix2world(pixcrd,0)
pixcrd3=wcsx.wcs_world2pix(pixcrd,0)
fig=plt.figure()
ax=WCSAxes(fig,[0.1,0.1,0.8,0.8],wcs=wcsx)
fig.add_axes(ax)
plt.imshow(np.flipud(signal),cmap="cubehelix")
print(pixcrd[0,:])
#plt.plot(pixcrd[:,0],pixcrd[:,1],color="red")
plt.plot(pixcrd2[:,0],pixcrd2[:,1],color="green")
#plt.plot(pixcrd3[:,0],pixcrd3[:,1],color="yellow")
def onclick2(event):
    x1,y1=event.xdata,event.ydata
    print(np.flipud(signal)[int(y1+.5)][int(x1+.5)])
cursor=Cursor(ax, horizOn=False,vertOn=False,color='red',linewidth=2.0)
fig.canvas.mpl_connect('button_press_event',onclick2)
plt.show()

'''
x = np.linspace(0, 10*np.pi, 100)
y = np.sin(x)

#plt.ion()
fig = plt.figure()
ax = fig.add_subplot(111)
line1, = ax.plot(x, y, 'b-')
print("oop1")
def update(val):
    amp = samp.val
    line1.set_ydata(amp*np.sin(x))
    #plt.pause(0.1)
    fig.canvas.draw()

ampax = plt.axes([0.2, 0.02, 0.65, 0.03], facecolor="beige")
samp = Slider(ampax, 'Amp', 0.1, 1, valinit=0.1, valstep=0.1,orientation="horizontal")
samp.on_changed(update)
print("oop2")
plt.show()
print("oop3")


    