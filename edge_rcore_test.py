import edge_detect
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt



for rcore in ["rcore4","rcore8","rcore12","rcore16","rcore20"]:
    
    edgenum=75
    edge=str(int(edgenum))
    #targ_levels=["target3","target5","target8","target10","target15"]
    #targ_levelsi=[3,5,8,10,15]
    targ_levels=["target3","target5","target8","target10","target15"]
    param_levels=["beta0.67","beta1.17","beta1.67","beta2.17","beta2.67","beta3.17"]
    param_levelsi=[0.67,1.17,1.67,2.17,2.67,3.17]
    targ_levelsi=[3,5,8,10,15]
    f,a=plt.subplots()

    uedge=np.zeros(len(param_levels))
    uedgev=np.zeros(len(param_levels))
    redge=np.zeros((len(targ_levelsi),len(param_levels)))
    redgev=np.zeros((len(targ_levelsi),len(param_levels)))

    for i in range(len(param_levels)):
        home_directory="/Volumes/TOSHIBA EXT/MARTIN2/simdata/betatest"
        sourcefilenames="testdata0"+"_"+rcore+"_signal.fits"
        sourcefilenamev="testdata0"+"_"+rcore+"_var.fits"
        blockfiles="block_testdata0"+"_"+rcore+"_wit_sig.fits"
        blockfilev="block_testdata0"+"_"+rcore+"_wit_var.fits"

        with fits.open(home_directory+"/"+param_levels[i]+"/unbinned/"+sourcefilenames) as hdul:
            ubs=hdul[0].data
            wcsx=hdul[0].header
        with fits.open(home_directory+"/"+param_levels[i]+"/unbinned/"+sourcefilenamev) as hdul:
            ubv=hdul[0].data
        uedge[i],uedgev[i]=edge_detect.radmethod(ubs,ubv,wcsx)
        for targi in range(len(targ_levels)):
            with fits.open(home_directory+"/"+param_levels[i]+"/"+targ_levels[targi]+"/"+blockfiles) as hdul:
                bs=hdul[0].data
                wcsx=hdul[0].header
            with fits.open(home_directory+"/"+param_levels[i]+"/"+targ_levels[targi]+"/"+blockfilev) as hdul:
                bv=hdul[0].data
            redge[targi,i],redgev[targi,i]=edge_detect.radmethod(bs,bv,wcsx)




    #a.set_xlabel("unbinned recovered radius")
    a.set_ylabel("binned recovered radius")
    a.set_xlabel("beta")
    for targi in range(len(targ_levelsi)):
        a.plot(param_levelsi,redge[targi,:],linewidth=1,marker="^",label="target "+str(targ_levelsi[targi]))
    a.plot([0,4],[75,75],"k",linestyle="dashed",label="true edge")
    a.plot(param_levelsi,uedge,"k",linestyle="dotted",marker="o",label="unbinned")
        
    a.legend(title=rcore)
    #plt.legend()
    #plt.show()
    plt.savefig(home_directory+"/edge_recovery_rcore"+rcore+".png")