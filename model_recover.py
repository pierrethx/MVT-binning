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
siglevel="sig10"
edge=str(int(edgenum))

targ_levels=["target3","target5","target8","target10","target15"]
cols=["tab:blue","tab:orange","tab:green","tab:red","tab:purple"]

unbinnedsig_m=np.zeros(numrange)
unbinnedvar_m=np.zeros(numrange)
blockedsig_m=np.zeros((len(targ_levels),numrange))
blockedvar_m=np.zeros((len(targ_levels),numrange))

with fits.open(home_directory+siglevel+"/unbinned/"+"testdata_edge"+edge+"_model.fits") as hdul:
    model=hdul[0].data

mask=np.full_like(model,1)
'''
mask=np.zeros_like(model)
xcent=int(160/2)
ycent=int(160/2)
for j in range(len(mask)):
    for i in range(len(mask[0])):
        if ((i-xcent)**2+(j-ycent)**2)<75**2:
            mask[j][i]=1

fig0,ax0=plt.subplots(1,1+len(targ_levels))
fig0.subplots_adjust(right=0.82)
th=np.linspace(0,6.2831,1000)
xc=75*np.cos(th)+80
yx=75*np.sin(th)+80
'''

for n in range(numrange):
        with fits.open(home_directory+siglevel+"/unbinned/"+"testdata"+str(n)+"_edge"+edge+"_signal.fits") as hdul:
            unbinnedsig_m[n]=np.sum(mask*(hdul[0].data-model)**2)
            '''
            if n==0:
                mi=np.min(mask*(hdul[0].data-model)**2)
                ma=np.max(mask*(hdul[0].data-model)**2)
                q=ax0[0].imshow(mask*(hdul[0].data-model)**2,vmin=mi,vmax=ma,cmap="cubehelix")
                ax0[0].set_xlabel("unbinned diff")
                ax0[0].set_xticklabels([])
                ax0[0].set_yticklabels([])
                ax0[0].plot(xc,yx,'r',linewidth=2)
            '''
                
        #with fits.open(home_directory+siglevel+"/unbinned/"+"testdata"+str(n)+"_edge"+edge+"_var.fits") as hdul:
        #    unbinnedvar_m[n]=np.sum((hdul[0].data-model)**2)
for targi in range(len(targ_levels)):
    for n in range(numrange):
        with fits.open(home_directory+siglevel+"/"+targ_levels[targi]+"/"+"block_testdata"+str(n)+"_edge"+edge+"_"+bintype+"_sig.fits") as hdul:
            blockedsig_m[targi,n]=np.sum(mask*(hdul[0].data-model)**2)
            '''
            print(blockedsig_m[targi,n])
            if n==0:
                ax0[targi+1].imshow(mask*(hdul[0].data-model)**2,vmin=mi,vmax=ma,cmap="cubehelix")
                ax0[targi+1].set_xlabel(targ_levels[targi]+" diff")
                ax0[targi+1].set_xticklabels([])
                ax0[targi+1].set_yticklabels([])
                ax0[targi+1].plot(xc,yx,'r',linewidth=2)
            '''
        #with fits.open(home_directory+siglevel+"/"+targ_levels[targi]+"/"+"block_testdata"+str(n)+"_edge"+edge+"_"+bintype+"_var.fits") as hdul:
        #    blockedvar_m[targi,n]=np.sum((hdul[0].data-model)**2)
'''
cb=fig0.add_axes([0.85, 0.35, 0.025, 0.30])
fig0.colorbar(q,cax=cb,orientation="vertical")
'''

unbinnedsig=np.mean(unbinnedsig_m,axis=0)
#unbinnedvar=np.mean(unbinnedvar_m,axis=1)
unbinnedsige=np.std(unbinnedsig_m,axis=0)
#unbinnedvare=np.std(unbinnedvar_m,axis=1)

blockedsig=np.mean(blockedsig_m,axis=1)

#blockedvar=np.mean(blockedvar_m,axis=2)
blockedsige=np.std(blockedsig_m,axis=1)
#blockedvare=np.std(blockedvar_m,axis=2)

fig,ax=plt.subplots()
goop=np.linspace(0,len(targ_levels)+1,10)
#ax.plot(goop,0*goop,linestyle="dashed",color="black")
#a.set_xlabel("unbinned recovered radius")
ax.set_ylabel("sum of square residuals")
print([unbinnedsig]+[blockedsig[targi] for targi in range(len(targ_levels))])
print(["unbinned"]+[targ_levels[targi] for targi in range(len(targ_levels))])
ax.boxplot([unbinnedsig_m]+[blockedsig_m[targi] for targi in range(len(targ_levels))],labels=["unbinned"]+[targ_levels[targi] for targi in range(len(targ_levels))])



plt.show()

     
