import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.widgets import Cursor,Slider,Button
#import tkinter
from tkinter.filedialog import askopenfilename
from astropy.io import fits
from astropy import wcs
from astropy.visualization.wcsaxes import WCSAxes
#import functions,bin_accretion,wvt_iteration
import scipy.spatial as sp

## The Sarazin circular beta model supposes that the brightness profile is like
## I=I_0 [1+(r/r_c)^2]^(0.5-3b)

I_0=1
r_c=1
b=.34

fix=.5
fix2=3

r=np.linspace(-fix2*2,fix2*2,10000)
I=I_0*(1+(r/r_c)**2)**(0.5-3*b)

#aI=I_0*(1+(0.5-3*b)*(r/r_c)**2)
eI=I_0*(np.exp((0.5-3*b)*(r/r_c)**2))
fI=(1+1/(6*b))**(0.5-3*b)

v=np.linspace(0,I_0,5)
stored=[]
fig = plt.figure()
ax = plt.axes([0.1, 0.3, .8,.6])
ax.set_ylim(0,I_0*1.05)
ax.set_xlim(-I_0,I_0)
#line0, = ax.plot(r, 0*r+fix, 'r',linestyle="dashed",alpha=0.5)
linev1, = ax.plot(v*0+fix2, v, 'r',linestyle="dotted",alpha=0.5)
linev2, = ax.plot(v*0-fix2, v, 'r',linestyle="dotted",alpha=0.5)
line1, = ax.plot(r, I, 'b-')
#linea1, = ax.plot(r, aI, 'k-')
#linee1, = ax.plot(r, eI, 'c-')
#linef, = ax.plot(r, 0*r+fI, 'k',linestyle="dashed")
print("oop1")

def update(val):
    r_c=srad.val
    b=sbet.val
    #I_0=samp.val
    I_0=1
    #r=np.linspace(-I_0,I_0,10000)
    #line1.set_xdata(r)

    I=I_0*(1+(r/r_c)**2)**(0.5-3*b)

    #aI=I_0*(1+(0.5-3*b)*(r/r_c)**2)
    eI=I_0*(np.exp((0.5-3*b)*(r/r_c)**2))
    fI=(1+1/(6*b))**(0.5-3*b)

    line1.set_ydata(I)
    #linea1.set_ydata(I_0*(1+(0.5-3*b)*(r/r_c)**2))
    #linee1.set_ydata(eI)

    #line0.set_xdata(r)
    #line0.set_ydata(0*r+fix)
    #linef.set_ydata(0*r+fI)
    linev1.set_xdata(v*0+fix2)
    linev2.set_xdata(v*0-fix2)
    #plt.pause(0.1)
    fig.canvas.draw()

def storit(val):
    r_c=srad.val
    b=sbet.val
    #I_0=samp.val
    I_0=1
    #r=np.linspace(-I_0,I_0,10000)
    #line1.set_xdata(r)

    I=I_0*(1+(r/r_c)**2)**(0.5-3*b)

    #aI=I_0*(1+(0.5-3*b)*(r/r_c)**2)
    eI=I_0*(np.exp((0.5-3*b)*(r/r_c)**2))
    fI=(1+1/(6*b))**(0.5-3*b)

    ax.plot(r,I,color="green",alpha=0.5)
    #ax.plot(r,I_0*(1+(0.5-3*b)*(r/r_c)**2),color="black",alpha=0.5)
    #ax.plot(r,eI,color="cyan",alpha=0.5)
    #ax.plot(r,0*r+fI,color="black",alpha=0.5)
    
"""
ampax = plt.axes([0.2, 0.16, 0.6, 0.03], facecolor="beige")
samp = Slider(ampax, 'Amp', 1, 20, valinit=2, valstep=0.1,orientation="horizontal")
samp.on_changed(update)"""

radax = plt.axes([0.2, 0.10, 0.6, 0.03], facecolor="beige")
srad = Slider(radax, 'Core Rad', 0.05, 1, valinit=1, valstep=0.01,orientation="horizontal")
srad.on_changed(update)

betax = plt.axes([0.2, 0.04, 0.6, 0.03], facecolor="beige")
sbet= Slider(betax, 'Beta', .34, 2.5, valinit=.34, valstep=0.01,orientation="horizontal")
sbet.on_changed(update)

buttax=plt.axes([0.4,0.15,.2,.05],facecolor="beige")
sbutt=Button(buttax,"Store")
sbutt.on_clicked(storit)

print("oop2")
plt.show()
print("oop3")