import numpy as np
import matplotlib.pyplot as plt
import tkinter
from tkinter.filedialog import askopenfilename
from astropy.io import fits
import bin_accretion
from astropy import wcs
from scipy.ndimage.measurements import center_of_mass


wcsx,signal,var,sourcedir,objname=bin_accretion.initialize(enternew=True)




ssum=0
svar=0
for y in range(len(signal)):
    for x in range(len(signal[0])):
        ssum+=signal[y][x]
        svar+=var[y][x]



print("objname: "+objname)
print("sig sum: "+str(ssum))
print("var sum: "+str(svar))
print("all ston "+str(ssum/np.sqrt(svar)))

target=[3,5]


for t in target:
    source="/".join(sourcedir.split("/")[:-1])
    name="_".join(objname.split("_")[:-1])
    sigfits=name+"_wit_sig.fits"
    varfits=name+"_wit_var.fits"
    stonfits="zston_"+name+"_wit.fits"
    signalname=source+"/target"+str(t)+"/"+sigfits
    with fits.open(signalname) as hdul:
            signal=np.flipud(hdul[0].data)
            wcsx=wcs.WCS(hdul[0].header)
    signalname=source+"/target"+str(t)+"/"+varfits
    with fits.open(signalname) as hdul:
            var=np.flipud(hdul[0].data)
            wcsx=wcs.WCS(hdul[0].header)

    ssum=0
    svar=0

    for y in range(len(signal)):
        for x in range(len(signal[0])):
            ssum+=signal[y][x]
            svar+=var[y][x]
    print("objname: /target"+str(t)+"/"+name)
    print("sig sum: "+str(ssum))
    print("var sum: "+str(svar))
    print("all ston "+str(ssum/np.sqrt(svar)))
    