import bin_accretion,main,functions
import numpy as np
import matplotlib.pyplot as plt

SMALL_SIZE=14
MEDIUM_SIZE=16
BIGGER_SIZE=20

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

def convergencelist(contarglist,diflistlist,targlist,sourcedir,objname,subfolder="unbinned"):
    ## generates a chart to see how a function converges over many iterations. 
    # This is to check condition 4 of the function converging over
    fpath=sourcedir+"/"+subfolder+"/y_"+objname+"_convergence.png"

    fig,ax=plt.subplots()
    ax.set_xlabel("iteration number")
    ax.set_ylabel("normalized difference")

    for i in range(len(contarglist)):
        n=len(diflistlist[i])
        diflist=np.where(np.isnan(diflistlist[i]),np.zeros_like(diflistlist[i]),diflistlist[i])[1:]
        
        xes=range(1,len(diflist)+1)
        if(contarglist[i]>0):
            ax.plot([1,len(diflist)],[contarglist[i],contarglist[i]],linestyle="dashed")
            
        ax.plot(xes,diflist,marker="o",label=targlist[i])
    ax.set_yscale('log')
    fig.legend(title="Target SNR")
    plt.tight_layout()
    
    plt.savefig(fpath)
    plt.close()

if __name__ == "__main__":

    targhold=0
    targlist=[]
    wcsxlist,siglist,varlist,sourcelist,objlist=bin_accretion.minitialize()
    contqueue=True
    while contqueue:
        try:
            target=main.gettarget(targhold)
            targlist.append(target)
        except:
            contqueue=False
    print("Files loaded!")
    for i in range(len(sourcelist)):
        epslist=[]
        diflistlist=[]
        for m in range(len(targlist)):
            
            wcsx=wcsxlist[i]
            signal=siglist[i]
            var=varlist[i]
            objname="_".join(objlist[i].split("_")[:-1])
            target=targlist[m]
            sourcedir=sourcelist[i]
            print(sourcedir)
            subfolder=main.makesubfolder(sourcedir,target)
            #main.saveston(wscxlist[i],siglist[i],varlist[i],sourcelist[i],objlist[i],subfolder="unbinned")

            ## this is us applying the data mask. We block out negative values and then run the binning algorithm on it
            signal2=np.copy(signal)
            var2=np.copy(var)
            var2[signal2<=0]=1e10
            signal2[signal2<=0]=0
            eps=-10

            binlist,diflist=main.mainfunc(signal2,var2,target,displayWVT=False,epsilon=eps)
            
            #main.saveblockoutfits(targlist[m],binlist,wscxlist[i],siglist[i],varlist[i],objlist[i],sourcelist[i],subfolder=subfolder)
            #wvt,ston=functions.generate_wvt2(binlist,siglist[i],varlist[i])
            ## then we apply the bins to the actual data and save all our files.
            wvt,ston=functions.generate_wvt3(binlist,signal,var,np.full(len(binlist),1))
            vwvt=functions.generate_wvt(binlist,var)
            main.saveiteratedfits(target,wcsx,wvt,vwvt,objname,sourcedir,subfolder=subfolder)
            main.saveblockoutfits(target,ston,wcsx,wvt,vwvt,objname,sourcedir,subfolder=subfolder,check=1)
            #main.saveblockoutoldfits(target,ston,wcsx,wvt,vwvt,objname,sourcedir,subfolder=subfolder)
            main.saveston(wcsx,ston,sourcedir,objname,subfolder=subfolder)
            assign=functions.assign(binlist,target,ston)
            main.saveassign(wcsx,assign,sourcedir,objname,subfolder=subfolder)
            epslist.append(eps)
            diflistlist.append(diflist)
        convergencelist(epslist,diflistlist,targlist,sourcedir,objname)
    print("Bye bye")