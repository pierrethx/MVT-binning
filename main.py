import numpy as np
import matplotlib.pyplot as plt 
import tkinter
from tkinter.filedialog import askopenfilename
from astropy.io import fits
import functions,bin_accretion,wvt_iteration
import scipy.spatial as sp
import time

import tkinter as tk
from tkinter import simpledialog

def getname(initial=""):
    application_window = tk.Tk()
    application_window.withdraw()
    target = simpledialog.askstring("Input", "What is the objectname?",parent=application_window,initialvalue=initial)
    if target is None:
        raise NameError("No name, goodbye.")
    return target

def gettarget(initial=0):
    application_window = tk.Tk()
    application_window.withdraw()
    target = simpledialog.askinteger("Input", "What is the target Signal-to-Noise? \n For cvt instead of wvt make this number negative :]",parent=application_window,initialvalue=initial)
    if target is None:
        raise NameError("No target, goodbye.")
    return target

def mainfunc(signal,var,target,weighting=True,displayWVT=True,epsilon=10):
    binlist,init_generators=bin_accretion.cc_accretion(signal,var,target)
    init_scalelengths=np.full(len(init_generators),1)
    binlist=wvt_iteration.iteration_func(target,signal,var,init_generators,init_scalelengths,epsilon,weighting=weighting,displaywvt=False)
    wvt=functions.generate_wvt(binlist,signal,displayWVT)
    vwvt=functions.generate_wvt(binlist,var,displayWVT)
    return wvt,vwvt

def saveiteratedfits(target,wcsx,wvt,vwvt,objname,sourcedir,subfolder="",weighting=True):
    header=wcsx.to_header()
    hdu = fits.PrimaryHDU(np.flipud(wvt),header=header)
    hdul = fits.HDUList([hdu])
    if weighting:
        hdul.writeto(sourcedir+subfolder+"/witerated_"+objname+"_target_"+str(target)+".fits",overwrite=True)
    else:
        hdul.writeto(sourcedir+subfolder+"/citerated_"+objname+"_target_"+str(target)+".fits",overwrite=True)

    hdu2 = fits.PrimaryHDU(np.flipud(vwvt),header=header)
    hdul2 = fits.HDUList([hdu2])
    if weighting:
        hdul2.writeto(sourcedir+subfolder+"/wit_var_"+objname+"_target_"+str(target)+".fits",overwrite=True)
    else:
        hdul2.writeto(sourcedir+subfolder+"/cit_var_"+objname+"_target_"+str(target)+".fits",overwrite=True)


if __name__ == "__main__":
    wcsx,signal,var,sourcedir,objname=bin_accretion.initialize(enternew=True)
    #objname="J024815-081723"
    ##empty for directly into sourcedir
    objname=getname(objname)
    target=gettarget()
    
    if target<0:
        weighting=False
        target=-target
    else:
        weighting=True
    wvt,vwvt=mainfunc(signal,var,target,displayWVT=True)
    saveiteratedfits(target,wcsx,wvt,vwvt,objname,sourcedir)
