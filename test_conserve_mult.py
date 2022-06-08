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
numrange=10
edge=str(int(edgenum))
blockt=False
sig_levels=["sig1","sig10","sig100","sig1000","sig10000"]
targ_levels=["target3","target5","target10"]
cols=["tab:blue","tab:orange","tab:green"]
markers=["s","v","*"]

unbinnedsig_m=np.zeros((len(sig_levels),numrange))
unbinnedvar_m=np.zeros((len(sig_levels),numrange))
if blockt:
    blockedsig_m=np.zeros((len(targ_levels),len(sig_levels),numrange))
    blockedvar_m=np.zeros((len(targ_levels),len(sig_levels),numrange))
else:
    unblockedsig_m=np.zeros((len(targ_levels),len(sig_levels),numrange))
    unblockedvar_m=np.zeros((len(targ_levels),len(sig_levels),numrange))

for n in range(numrange):
    for sigi in range(len(sig_levels)):
        with fits.open(home_directory+sig_levels[sigi]+"/unbinned/"+"testdata"+str(n)+"_edge"+edge+"_signal.fits") as hdul:
            unbinnedsig_m[sigi,n]=np.sum(hdul[0].data)
        with fits.open(home_directory+sig_levels[sigi]+"/unbinned/"+"testdata"+str(n)+"_edge"+edge+"_var.fits") as hdul:
            unbinnedvar_m[sigi,n]=np.sum(hdul[0].data)
        for targi in range(len(targ_levels)):
            if blockt:
                with fits.open(home_directory+sig_levels[sigi]+"/"+targ_levels[targi]+"/"+"block_testdata"+str(n)+"_edge"+edge+"_"+bintype+"_sig.fits") as hdul:
                    blockedsig_m[targi,sigi,n]=np.sum(hdul[0].data)
                with fits.open(home_directory+sig_levels[sigi]+"/"+targ_levels[targi]+"/"+"block_testdata"+str(n)+"_edge"+edge+"_"+bintype+"_var.fits") as hdul:
                    blockedvar_m[targi,sigi,n]=np.sum(hdul[0].data)
            else:
                with fits.open(home_directory+sig_levels[sigi]+"/"+targ_levels[targi]+"/"+"testdata"+str(n)+"_edge"+edge+"_"+bintype+"_sig.fits") as hdul:
                    unblockedsig_m[targi,sigi,n]=np.sum(hdul[0].data)
                with fits.open(home_directory+sig_levels[sigi]+"/"+targ_levels[targi]+"/"+"testdata"+str(n)+"_edge"+edge+"_"+bintype+"_var.fits") as hdul:
                    unblockedvar_m[targi,sigi,n]=np.sum(hdul[0].data)

unbinnedsig=np.mean(unbinnedsig_m,axis=1)
unbinnedvar=np.mean(unbinnedvar_m,axis=1)
unbinnedsige=np.std(unbinnedsig_m,axis=1)
unbinnedvare=np.std(unbinnedvar_m,axis=1)

if blockt:
    blockedsig=np.mean(blockedsig_m,axis=2)
    blockedvar=np.mean(blockedvar_m,axis=2)
    blockedsige=np.std(blockedsig_m,axis=2)
    blockedvare=np.std(blockedvar_m,axis=2)
else:
    unblockedsig=np.mean(unblockedsig_m,axis=2)
    unblockedvar=np.mean(unblockedvar_m,axis=2)
    unblockedsige=np.std(unblockedsig_m,axis=2)
    unblockedvare=np.std(unblockedvar_m,axis=2)




fig,ax=plt.subplots(2,1,gridspec_kw={'height_ratios': [3, 1]})
fig2,ax2=plt.subplots(2,1,gridspec_kw={'height_ratios': [3, 1]})
for t in range(len(targ_levels)):
    if blockt:
        ax[0].set_title("blocked_testdata_edge"+edge+"_"+bintype+"_sig.fits")
        ax2[0].set_title("blocked_testdata_edge"+edge+"_"+bintype+"_var.fits")
        #ax[0].plot(unbinnedsig,blockedsig[t,:],marker="o",linewidth=0,color=cols[t],label=targ_levels[t])
        ax[0].errorbar(unbinnedsig,blockedsig[t,:],blockedsige[t,:],unbinnedsige,marker="o",capsize=2,linewidth=0,color=cols[t],label=targ_levels[t])
        #ems=(blockedsig[t,:]/unbinnedsig)*np.sqrt( (blockedsige[t,:]/blockedsig[t,:])**2+(unbinnedsige/unbinnedsig)**2)
        #ax[1].errorbar(unbinnedsig,blockedsig[t,:]/unbinnedsig-1,ems,unbinnedsige,marker="o",linewidth=0)
        ax[1].plot(unbinnedsig,blockedsig[t,:]/unbinnedsig-1,marker="o",linewidth=0)
        ax[0].set_ylabel("binned, blocked sig")
        ax2[0].errorbar(unbinnedvar,blockedvar[t,:],blockedvare[t,:],unbinnedvare,marker="o",linewidth=0,color=cols[t],label=targ_levels[t])
        #emv=(blockedvar[t,:]/unbinnedvar)*np.sqrt( (blockedvare[t,:]/blockedvar[t,:])**2+(unbinnedvare/unbinnedvar)**2)
        #ax2[1].errorbar(unbinnedvar,blockedvar[t,:]/unbinnedvar-1,emv,unbinnedvare,marker="o",color=cols[t],linewidth=1)
        ax2[1].plot(unbinnedvar,blockedvar[t,:]/unbinnedvar-1,marker="o",color=cols[t],linewidth=1)
        ax[0].set_ylabel("binned, blocked var")
    else:
        ax[0].set_title("testdata_edge"+edge+"_"+bintype+"_sig.fits")
        ax2[0].set_title("testdata_edge"+edge+"_"+bintype+"_var.fits")
        #ax[0].plot(unbinnedsig,blockedsig[t,:],marker="o",linewidth=0,color=cols[t],label=targ_levels[t])
        ax[0].errorbar(unbinnedsig,unblockedsig[t,:],unblockedsige[t,:],unbinnedsige,marker=markers[t],capsize=2,linewidth=0,color=cols[t],label=targ_levels[t])
        #ems=(blockedsig[t,:]/unbinnedsig)*np.sqrt( (blockedsige[t,:]/blockedsig[t,:])**2+(unbinnedsige/unbinnedsig)**2)
        #ax[1].errorbar(unbinnedsig,blockedsig[t,:]/unbinnedsig-1,ems,unbinnedsige,marker="o",linewidth=0)
        ax[1].plot(unbinnedsig,unblockedsig[t,:]/unbinnedsig-1,marker=markers[t],linewidth=1)
        ax[0].set_ylabel("binned sig")
        ax2[0].errorbar(unbinnedvar,unblockedvar[t,:],unblockedvare[t,:],unbinnedvare,marker=markers[t],linewidth=0,color=cols[t],label=targ_levels[t])
        #emv=(blockedvar[t,:]/unbinnedvar)*np.sqrt( (blockedvare[t,:]/blockedvar[t,:])**2+(unbinnedvare/unbinnedvar)**2)
        #ax2[1].errorbar(unbinnedvar,blockedvar[t,:]/unbinnedvar-1,emv,unbinnedvare,marker="o",color=cols[t],linewidth=1)
        ax2[1].plot(unbinnedvar,unblockedvar[t,:]/unbinnedvar-1,marker=markers[t],color=cols[t],linewidth=1)
        ax2[0].set_ylabel("binned var")
    
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
