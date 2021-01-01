import numpy as np
import matplotlib.pyplot as plt
import tkinter
from tkinter.filedialog import askopenfilename
from astropy.io import fits
import bin_accretion
from astropy import wcs
from scipy.ndimage.measurements import center_of_mass

import scipy.spatial as sp 

scaleup=40

bloop=True
alist=[]
while (bloop):
    try:
        root=tkinter.Tk()
        root.withdraw()
        placeholder= askopenfilename(message="Select assign",multiple=False)
        root.update()
        root.destroy()
        with fits.open(placeholder) as hdul:
            assign=np.flipud(hdul[0].data).astype(int)
            wcsx=wcs.WCS(hdul[0].header)
        #alist.append(assign)
        aa=np.zeros((scaleup*len(assign),scaleup*len(assign[0])),dtype=int)
        for l in range(len(aa)):
            for m in range(len(aa[l])):
                aa[l][m]=assign[int(l/scaleup)][int(m/scaleup)]
        alist.append(aa)
    except:
        bloop=False

    

ratlist=[]
radlist=[]

for assign in alist:

    center=[len(assign)/2,len(assign[0])/2]


    bins=(np.nanmax(assign))

    bcenters=[[0,0] for i in range(bins)]

    binlist=[[] for i in range(bins)]
    for y in range(len(assign)):
        for x in range(len(assign[y])):
            binlist[assign[y][x]-1].append([y,x])
            bcenters[assign[y][x]-1][0]+=y
            bcenters[assign[y][x]-1][1]+=x
        print(str(y/len(assign))+" through stage 1")

    blue=[0]*bins
    brads=np.array(blue,copy=True)
    #reff=[0]*bins
    bareas=np.array(blue,copy=True)
    ratios=[0]*bins

    radius=[0]*bins
    for m in range(0, len(binlist)):
        bareas[m]=len(binlist[m])
        bcenters[m][0]/=bareas[m]
        bcenters[m][1]/=bareas[m]
        
        brads[m]=sum(sp.distance.cdist([bcenters[m]],binlist[m])[0])/bareas[m]

        print(str(m/len(binlist))+" through stage 2")
    reff=np.sqrt(bareas/3.14159)
    ratios=(brads)/reff
    radius=sp.distance.cdist([center],bcenters)[0]
    print("done with stage 3")
    ratlist.append(ratios)
    radlist.append(radius/scaleup)
fig,ax=plt.subplots()
ax.set_xlabel("radius from center")
ax.set_ylabel("Rav/Reff")
labelist=["target 3","target 5"]
for i in range(len(ratlist)):
    bot,top=plt.ylim()
    ax.set_ylim(2/3,top)
    plt.plot(radlist[i],ratlist[i],linewidth=0,marker="o",label=labelist[i])
fig.legend()
plt.show()
