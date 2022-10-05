import time, os,sys
sys.path.append(("/".join(os.path.realpath(__file__).split("/")[:-2])))


import numpy as np
import matplotlib.pyplot as plt 
import tkinter as tk
from tkinter.filedialog import askopenfilename
from astropy.io import fits
from astropy import wcs
import functions,bin_accretion,wvt_iteration,main
import scipy.spatial as sp
from scipy import ndimage

from collections import Counter

SMALL_SIZE=14
MEDIUM_SIZE=16
BIGGER_SIZE=20

## this is meant to be a faithful implementation of the original CVT and WVT but using the
## tech of the new implementation. Specifically, we use a minimum input SNR  and thus
## do not need an output SNR cutoff
## Cappellari and Copin 2003, Diehl and Statler 2006


if __name__ == "__main__":

    ## negative epsilon means iterate that number of times
    ## positive epsilon gives us a tolerance our iterations must get below
    eps=-30

    ## minimum input SNR=1 means we use all pixels with more signal than noise
    minSNR=1

    modetypes=[]

    #modetypes.append("WVT")
    modetypes.append("CVT")
    
    ## here for easy testing we allow you to select multiple options but 
    ## unless you were actually using multiple options you would just comment the rest out

    targhold="3,5,10"
    wcsxlist,siglist,varlist,sourcelist,objlist=bin_accretion.minitialize()
    if len(siglist)==0:
        print("No files loaded!")
    else:
        print("Files loaded!")
        targlist=bin_accretion.getmulttarget(targhold)
        for i in range(len(sourcelist)):
            for modetype in modetypes:
                epslist=[]
                diflistlist=[]
                wcsx=wcsxlist[i]
                signal=siglist[i]
                var=varlist[i]
                if "_sig.fits" in objlist[i].lower():
                    objname="_".join(objlist[i].split("_")[:-1])
                else:
                    objname=".".join(objlist[i].split(".")[:-1])
                sourcedir=sourcelist[i]
                print(sourcedir)
                for m in range(len(targlist)):

                    target=targlist[m]
                    
                    subfolder=bin_accretion.makesubfolder(sourcedir,target,modetype+"0")

                    ## this is us applying the data mask. We block out negative values and then run the binning algorithm on it
                    signal2=np.copy(signal)
                    var2=np.copy(var)
                    var2[signal2<=0]=1e10
                    var2[var2<=0]=1e10
                    ston2=signal2/np.sqrt(var2)

                    mask=np.full_like(signal2,1)
                    mask[ston2<minSNR]=0
                    
                    ## this is the most important function
                    ## first we do bin accretion. then iteration. the wvt and vwvt is not necessary here but im scared to remove it
                    binlist,init_generators,init_scalelengths=bin_accretion.cc_accretion(signal,var,target,minsize=1,mode=modetype,mask=mask,display=False)
                    
                    binlist,diflist=wvt_iteration.iteration_func(target,signal,var,init_generators,init_scalelengths,eps,modetype,display=False,mask=mask)

                    ## now to generate the binned signal and variance
                    wvt,ston=functions.generate_wvt3(binlist,signal*mask,var*mask,np.full(len(binlist),1))
                    vwvt=functions.generate_wvt(binlist,var*mask) 

                    main.saveunblockedfits(wcsx,wvt,vwvt,objname,sourcedir,modetype+"0",subfolder)
                    main.savestonimg(target,ston,sourcedir,objname,subfolder,modetype+"0")
                    main.saveston(wcsx,ston,sourcedir,objname,modetype+"0",subfolder)
                    assig=main.assign(binlist,target,ston)
                    main.saveassign(wcsx,assig,sourcedir,objname,modetype+"0",subfolder)
                    epslist.append(eps)
                    diflistlist.append(diflist)
                if(len(diflistlist[0])>0):
                    main.convergencelist(epslist,diflistlist,targlist,sourcelist[i],objname,modetype+"0")
    print("Bye bye")