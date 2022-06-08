
import edge_detect
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt



home_directory="/Volumes/TOSHIBA EXT/MARTIN2/simdata/123test/"

edgenum=75
edge=str(int(edgenum))
sig_levels=["sig1","sig10","sig100","sig1000","sig10000"]
targ_levels=["target3","target5","target10"]
targ_levelsi=[3,5,8,10,15]
sourcefilenames="testdata0_edge"+edge+"_signal.fits"
sourcefilenamev="testdata0_edge"+edge+"_var.fits"
blockfiles="block_testdata0_edge"+edge+"_wit_sig.fits"
blockfilev="block_testdata0_edge"+edge+"_wit_var.fits"

uedge=np.zeros(len(sig_levels))
uedgev=np.zeros(len(sig_levels))
redge=np.zeros((len(sig_levels),len(targ_levels)))
redgev=np.zeros((len(sig_levels),len(targ_levels)))

for sigi in range(len(sig_levels)):
    with fits.open(home_directory+sig_levels[sigi]+"/unbinned/"+sourcefilenames) as hdul:
        ubs=hdul[0].data
        wcsx=hdul[0].header
    with fits.open(home_directory+sig_levels[sigi]+"/unbinned/"+sourcefilenamev) as hdul:
        ubv=hdul[0].data
    uedge[sigi],uedgev[sigi]=edge_detect.radmethod(ubs,ubv,wcsx)
    for targi in range(len(targ_levels)):
        with fits.open(home_directory+sig_levels[sigi]+"/"+targ_levels[targi]+"/"+blockfiles) as hdul:
            bs=hdul[0].data
            wcsx=hdul[0].header
        with fits.open(home_directory+sig_levels[sigi]+"/"+targ_levels[targi]+"/"+blockfilev) as hdul:
            bv=hdul[0].data
        redge[sigi,targi],redgev[sigi,targi]=edge_detect.radmethod(bs,bv,wcsx)

f,a=plt.subplots(len(sig_levels),1)
for ai in range(len(sig_levels)):
    a[ai].set_title(sig_levels[ai])
    a[ai].errorbar(0,uedge[ai],uedgev[ai],linestyle='dashed',marker='o',color="red",label="unbinned")
    a[ai].plot([0,15],[edgenum,edgenum],linestyle='dashed',color="blue",label="true")
    a[ai].errorbar(targ_levelsi,redge[ai,:],redgev[ai,:],marker="o")
    #a[ai].legend()
plt.show()