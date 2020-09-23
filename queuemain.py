import bin_accretion,main,functions
import numpy as np

if __name__ == "__main__":
    sourcelist=[]
    wscxlist=[]
    siglist=[]
    varlist=[]
    objlist=[]
    targlist=[]
    contqueue=True
    targhold=0
    
    while contqueue:
        try:
            wcsx,signal,var,sourcedir,objname=bin_accretion.initialize(enternew=True)
            objname=main.getname("_".join(objname.split("_")[:-1]))
            sourcelist.append("/".join(sourcedir.split("/")[:-1]))
            wscxlist.append(wcsx)
            siglist.append(signal)
            varlist.append(var)
            objlist.append(objname)
        except:
            contqueue=False
    contqueue=True
    while contqueue:
        try:
            target=main.gettarget(targhold)
            targlist.append(target)
        except:
            contqueue=False
    print("Files loaded!")
    for i in range(len(sourcelist)):
        for m in range(len(targlist)):
            if targlist[m]<0:
                #weighting=False 
                """Would be false except i havent fixed cvt"""
                weighting=True
                targlist[m]=-targlist[m]
            else:
                weighting=True
            subfolder="target"+str(targlist[m])
            #main.saveston(wscxlist[i],siglist[i],varlist[i],sourcelist[i],objlist[i],subfolder="unbinned")

            signal2=np.copy(siglist[i])
            var2=np.copy(varlist[i])
            #var2[signal2<=0]=1e10
            signal2[signal2<=0]=0

            binlist=main.mainfunc(signal2,var2,targlist[m],displayWVT=False,epsilon=-10)
            
            #main.saveblockoutfits(targlist[m],binlist,wscxlist[i],siglist[i],varlist[i],objlist[i],sourcelist[i],subfolder=subfolder)
            #wvt,ston=functions.generate_wvt2(binlist,siglist[i],varlist[i])
            wvt,ston=functions.generate_wvt3(binlist,siglist[i],varlist[i],np.full(len(binlist),1))
            vwvt=functions.generate_wvt(binlist,varlist[i])
            #main.saveiteratedfits(targlist[m],binlist,wscxlist[i],wvt,vwvt,objlist[i],sourcelist[i],subfolder=subfolder)
            main.saveblockoutfits(targlist[m],ston,wscxlist[i],wvt,vwvt,objlist[i],sourcelist[i],subfolder=subfolder)
            main.saveston(wscxlist[i],ston,sourcelist[i],objlist[i],subfolder=subfolder)
            assign=functions.assign(binlist,siglist[i])
            main.saveassign(wscxlist[i],assign,sourcelist[i],objlist[i],subfolder=subfolder)

print("Bye Bye!")