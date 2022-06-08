import edge_detect
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt


home_directory="/Volumes/TOSHIBA EXT/MARTIN2/simdata/betatest/"
edgenum=75
edge=str(int(edgenum))

for rcore in ["rcore4","rcore8","rcore12","rcore16","rcore20"]:

    #targ_levels=["target3","target5","target8","target10","target15"]
    #targ_levelsi=[3,5,8,10,15]
    targ_levels=["target3","target5","target8","target10","target15"]
    param_levels=["beta0.67","beta1.17","beta1.67","beta2.17","beta2.67","beta3.17"]
    param_levelsi=[0.67,1.17,1.67,2.17,2.67,3.17]
    targ_levelsi=[3,5,8,10,15]
    f,a=plt.subplots()

    #cols=["tab:blue","tab:orange","tab:green","tab:red","tab:purple"]

    unbinnedsig_m=np.zeros(len(param_levels))
    unbinnedvar_m=np.zeros(len(param_levels))
    blockedsig_m=np.zeros((len(targ_levels),len(param_levels)))
    blockedvar_m=np.zeros((len(targ_levels),len(param_levels)))
    for q in range(len(param_levels)):
        with fits.open(home_directory+"/"+param_levels[q]+"/unbinned/testdata_"+rcore+"_model.fits") as hdul:
            model=hdul[0].data

        mask=np.full_like(model,1)
        with fits.open(home_directory+param_levels[q]+"/unbinned//testdata0_"+rcore+"_signal.fits") as hdul:
            unbinnedsig_m[q]=np.sum(mask*(hdul[0].data-model)**2)
  #    unbinnedvar_m[n]=np.sum((hdul[0].data-model)**2)
        for targi in range(len(targ_levels)):
            with fits.open(home_directory+param_levels[q]+"/"+targ_levels[targi]+"/block_testdata0_"+rcore+"_wit_sig.fits") as hdul:
                blockedsig_m[targi][q]=np.sum(mask*(hdul[0].data-model)**2)

    #a.set_xlabel("unbinned recovered radius")
    a.set_ylabel("summed squared residuals")
    a.set_xlabel("beta")
    for targi in range(len(targ_levelsi)):
        a.plot(param_levelsi,blockedsig_m[targi,:],linewidth=1,marker="^",label="target "+str(targ_levelsi[targi]))
    a.plot(param_levelsi,unbinnedsig_m,"k",marker="o",linestyle="dashed",label="unbinned")  
    a.legend(title=rcore)
    #plt.legend()
    #plt.show()
    plt.savefig(home_directory+"/model_recovery_"+rcore+".png",bbox_inches="tight")