import numpy as np
import matplotlib.pyplot as plt
import tkinter
from tkinter.filedialog import askopenfilename
from astropy.io import fits
import bin_accretion
from astropy import wcs
from scipy.ndimage.measurements import center_of_mass

def trim(lis,targ):
    other=[]
    for i in range(len(lis)):
        if lis[i]<=targ:
            other.append(lis[i])
    return other

# make sig noise -> ston
root=tkinter.Tk()
root.withdraw()
signalname= askopenfilename(message="Select signal")
root.update()
root.destroy()


fig,ax=plt.subplots()
xun=[]
yun=[]


with fits.open(signalname) as hdul:
    ston=np.flipud(hdul[0].data)
    wcsx=hdul[0].header

center=center_of_mass(np.ma.masked_where(np.isnan(ston),ston)) #as usual is (y,x)
scale=0.5*np.sqrt(len(ston)**2+len(ston[0])**2)
mask=np.copy(ston)
for yy in range(len(mask)):
    for xx in range(len(mask[yy])):
        if np.sqrt((yy-center[0])**2+(xx-center[1])**2)>scale/3:
            mask[yy][xx]=0
center=center_of_mass(mask)

for y in range(len(ston)):
    for x in range(len(ston[0])):
        if not np.isnan(ston[y][x]):
            xun.append(np.sqrt((y-center[0])**2+(x-center[1])**2))
            yun.append(ston[y][x])
ax.plot(xun,np.full(len(xun),0),linestyle="dashed",color="red")
ax.plot(xun,yun,linewidth=0,marker=".",label="unbinned")
ax.set_xlabel("Radius from center")
ax.set_ylabel("Signal")

plt.show()