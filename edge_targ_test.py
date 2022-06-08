import edge_detect
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt



home_directory="/Volumes/TOSHIBA EXT/MARTIN2/simdata/123test/sig10/"

numm=10
edgenum=75
edge=str(int(edgenum))
#targ_levels=["target3","target5","target8","target10","target15"]
#targ_levelsi=[3,5,8,10,15]
targ_levels=["target3","target5","target8","target10","target15"]
targ_levelsi=[3,5,8,10,15]
f,a=plt.subplots()

uedge=np.zeros(numm)
uedgev=np.zeros(numm)
redge=np.zeros((len(targ_levelsi),numm))
redgev=np.zeros((len(targ_levelsi),numm))

for num in range(numm):
    sourcefilenames="testdata"+str(num)+"_edge"+edge+"_signal.fits"
    sourcefilenamev="testdata"+str(num)+"_edge"+edge+"_var.fits"
    blockfiles="block_testdata"+str(num)+"_edge"+edge+"_wit_sig.fits"
    blockfilev="block_testdata"+str(num)+"_edge"+edge+"_wit_var.fits"

    with fits.open(home_directory+"/unbinned/"+sourcefilenames) as hdul:
        ubs=hdul[0].data
        wcsx=hdul[0].header
    with fits.open(home_directory+"/unbinned/"+sourcefilenamev) as hdul:
        ubv=hdul[0].data
    uedge[num],uedgev[num]=edge_detect.radmethod(ubs,ubv,wcsx)
    for targi in range(len(targ_levels)):
        with fits.open(home_directory+"/"+targ_levels[targi]+"/"+blockfiles) as hdul:
            bs=hdul[0].data
            wcsx=hdul[0].header
        with fits.open(home_directory+"/"+targ_levels[targi]+"/"+blockfilev) as hdul:
            bv=hdul[0].data
        redge[targi,num],redgev[targi,num]=edge_detect.radmethod(bs,bv,wcsx)



goop=np.linspace(0,len(targ_levels)+1,10)
a.plot(goop,edgenum+0*goop,linestyle="dashed",color="black")
#a.set_xlabel("unbinned recovered radius")
a.set_ylabel("recovered radius")
a.boxplot([uedge]+[redge[targi,:] for targi in range(len(targ_levels))],labels=["unbinned"]+[targ_levels[targi] for targi in range(len(targ_levels))])

fi,ax=plt.subplots()
ies=np.argsort(uedge)
for targi in range(len(targ_levels)):
    
    ax.errorbar(uedge[ies],redge[targi,:][ies],np.sqrt(redgev[targi,:])[ies],np.sqrt(uedgev[ies]),marker=".",linewidth=0,elinewidth=1,label=targ_levels[targi])
    #ax.plot(uedge[ies],redge[targi,:][ies],marker=".",linewidth=2,label=targ_levels[targi])
    #ax.plot(uedge[ies],edgenum+0*uedge,linestyle="dashed",color="black")
    #a.boxplot(redge[targi,:],positions=targi+1,labels=targ_levels[targi])
    #pass
ax.legend()
#plt.legend()
plt.show()