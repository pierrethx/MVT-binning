import numpy as np
import matplotlib.pyplot as plt
import tkinter
from tkinter.filedialog import askopenfilename
from astropy.io import fits
import bin_accretion
from astropy import wcs
from scipy.ndimage.measurements import center_of_mass

import scipy.spatial as sp 


root=tkinter.Tk()
root.withdraw()
placeholder= askopenfilename(message="Select assign",multiple=False)
root.update()
root.destroy()
with fits.open(placeholder) as hdul:
    assign=np.flipud(hdul[0].data).astype(int)
    wcsx=wcs.WCS(hdul[0].header)

center=[len(assign)/2,len(assign[0])/2]


bins=(np.nanmax(assign))

bcenters=[[0,0]]*bins

binlist=[[]]*(bins+1)
for y in range(len(assign)):
    for x in range(len(assign[y])):
        binlist[assign[y][x]].append([y,x])
        bcenters[assign[y][x]-1][0]+=y
        bcenters[assign[y][x]-1][1]+=x
    print(str(y/len(assign))+" through stage 1")

blue=[0]*bins
brads=np.array(blue,copy=True)
#reff=[0]*bins
bareas=np.array(blue,copy=True)
ratios=[0]*bins

radius=[0]*bins
for m in range(1, len(binlist)):
    bareas[m-1]=len(binlist[m])
    bcenters[m-1][0]/=bareas[m-1]
    bcenters[m-1][1]/=bareas[m-1]
    
    brads[m-1]=sum(sp.distance.cdist([bcenters[m-1]],binlist[m])[0])/bareas[m-1]

    print(str(m/len(binlist))+" through stage 2")
reff=np.sqrt(bareas/3.14159)
ratios=brads/bareas
radius=sp.distance.cdist([center],bcenters)[0]
print("done with stage 3")
print(radius)
'''
fig,ax=plt.subplots()
ax.set_xlabel("radius from center")
ax.set_ylabel("Rav/Reff")
plt.plot(radius,ratios,linewidth=0,marker="o")
plt.show()
'''
