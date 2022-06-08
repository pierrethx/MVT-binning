import numpy as np
import matplotlib.pyplot as plt
from tkinter.filedialog import askopenfilename
from astropy.io import fits
import bin_accretion,functions
from astropy import wcs
import os

home_directory="/Volumes/TOSHIBA EXT/MARTIN2/simdata/123test/"
bintype="wit"
edgenum=75
edge=str(int(edgenum))
blockt=False
sig_levels=["sig1","sig10","sig100","sig1000","sig10000"]
targ_levels=["target3","target5","target10"]
cols=["tab:blue","tab:orange","tab:green"]

sourcefilenames="testdata0_edge"+edge+"_signal.fits"
sourcefilenamev="testdata0_edge"+edge+"_var.fits"
if blockt:
    blockfiles="block_testdata0_edge"+edge+"_"+bintype+"_sig.fits"
    blockfilev="block_testdata0_edge"+edge+"_"+bintype+"_var.fits"
else:
    unblockfiles="testdata0_edge"+edge+"_"+bintype+"_sig.fits"
    unblockfilev="testdata0_edge"+edge+"_"+bintype+"_var.fits"

unbinnedsig=np.zeros(len(sig_levels))
unbinnedvar=np.zeros(len(sig_levels))
if blockt:
    blockedsig=np.zeros((len(targ_levels),len(sig_levels)))
    blockedvar=np.zeros((len(targ_levels),len(sig_levels)))
else:
    binnedsig=np.zeros((len(targ_levels),len(sig_levels)))
    binnedvar=np.zeros((len(targ_levels),len(sig_levels)))


for sigi in range(len(sig_levels)):
    with fits.open(home_directory+sig_levels[sigi]+"/unbinned/"+sourcefilenames) as hdul:
        unbinnedsig[sigi]=np.sum(hdul[0].data)
    with fits.open(home_directory+sig_levels[sigi]+"/unbinned/"+sourcefilenamev) as hdul:
        unbinnedvar[sigi]=np.sum(hdul[0].data)
    for targi in range(len(targ_levels)):
        if blockt:
            with fits.open(home_directory+sig_levels[sigi]+"/"+targ_levels[targi]+"/"+blockfiles) as hdul:
                blockedsig[targi,sigi]=np.sum(hdul[0].data)
            with fits.open(home_directory+sig_levels[sigi]+"/"+targ_levels[targi]+"/"+blockfilev) as hdul:
                blockedvar[targi,sigi]=np.sum(hdul[0].data)
        else:
            with fits.open(home_directory+sig_levels[sigi]+"/"+targ_levels[targi]+"/"+unblockfiles) as hdul:
                binnedsig[targi,sigi]=np.sum(hdul[0].data)
            with fits.open(home_directory+sig_levels[sigi]+"/"+targ_levels[targi]+"/"+unblockfilev) as hdul:
                binnedvar[targi,sigi]=np.sum(hdul[0].data)


fig,ax=plt.subplots(2,1,gridspec_kw={'height_ratios': [3, 1]})
fig2,ax2=plt.subplots(2,1,gridspec_kw={'height_ratios': [3, 1]})
for t in range(len(targ_levels)):
    
    
    if blockt:
        ax[0].set_title(blockfiles)
        ax2[0].set_title(blockfilev)
        ax[0].plot(unbinnedsig,blockedsig[t,:],marker="o",linewidth=0,color=cols[t],label=targ_levels[t])
        ax[1].plot(unbinnedsig,blockedsig[t,:]/unbinnedsig-1,marker="o",linewidth=0,label="block")
        ax[0].set_ylabel("binned, blocked sig")
        ax2[0].plot(unbinnedvar,blockedvar[t,:],marker="o",linewidth=0,color=cols[t],label=targ_levels[t])
        ax2[1].plot(unbinnedvar,blockedvar[t,:]/unbinnedvar-1,marker="o",color=cols[t],linewidth=1)
        ax[0].set_ylabel("binned, blocked var")
    else:
        ax[0].set_title(unblockfiles)
        ax2[0].set_title(unblockfilev)
        ax[0].plot(unbinnedsig,binnedsig[t,:],marker="o",linewidth=0,color=cols[t],label=targ_levels[t])
        ax[1].plot(unbinnedsig,binnedsig[t,:]/unbinnedsig-1,marker="o",color=cols[t],linewidth=1)
        ax[0].set_ylabel("binned, unblocked sig")
        ax2[0].plot(unbinnedvar,binnedvar[t,:],marker="o",linewidth=0,color=cols[t],label=targ_levels[t])
        ax2[1].plot(unbinnedvar,binnedvar[t,:]/unbinnedvar-1,marker="o",color=cols[t],linewidth=1)
        ax[0].set_ylabel("binned, unblocked var")
    ax[0].plot([np.min(unbinnedsig),np.max(unbinnedsig)],[np.min(unbinnedsig),np.max(unbinnedsig)],color="k",linestyle="dashed")
    ax[1].plot([np.min(unbinnedsig),np.max(unbinnedsig)],[0,0],color="k",linestyle="dashed")
    
    ax[0].get_xaxis().set_visible(False)
    ax[0].set_xscale("log")
    ax[0].set_yscale("log")
    ax[1].set_xscale("log")
    ax[1].set_xlabel("unbinned signal")
    ax[1].set_ylabel("normalized residual")
    #ax[0,t].set_yscale("log")
    ax[0].legend()

    
    
    #ax[t].plot(unbinnedsig,blockedsig[t,:]/unbinnedsig-1,marker="o",linewidth=0,label="block")
    ax2[0].plot([np.min(unbinnedvar),np.max(unbinnedvar)],[np.min(unbinnedvar),np.max(unbinnedvar)],color="k",linestyle="dashed")
    ax2[1].plot([np.min(unbinnedvar),np.max(unbinnedvar)],[0,0],color="k",linestyle="dashed")
    ax2[0].get_xaxis().set_visible(False)
    ax2[0].set_xscale("log")
    ax2[0].set_yscale("log")
    ax2[1].set_xscale("log")
    ax2[1].set_xlabel("unbinned var")
    ax2[1].set_ylabel("normalized residual")
    #ax[0,t].set_yscale("log")
    ax2[0].legend()

fig.tight_layout()
fig2.tight_layout()
if blockt:
    fig.savefig(home_directory+"b_testdata_edge"+edge+"_"+bintype+"_sig_conservation.png")
    fig2.savefig(home_directory+"b_testdata_edge"+edge+"_"+bintype+"_var_conservation.png")
else:
    fig.savefig(home_directory+"bn_testdata_edge"+edge+"_"+bintype+"_sig_conservation.png")
    fig2.savefig(home_directory+"bn_testdata_edge"+edge+"_"+bintype+"_var_conservation.png")
plt.show()
