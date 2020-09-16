import bin_accretion,main,functions

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
            binlist=main.mainfunc(siglist[i],varlist[i],targlist[m],weighting=weighting,displayWVT=False,epsilon=-10)
            
            #main.saveblockoutfits(targlist[m],binlist,wscxlist[i],siglist[i],varlist[i],objlist[i],sourcelist[i],subfolder=subfolder)
            wvt,ston=functions.generate_wvt2(binlist,siglist[i],varlist[i])
            vwvt=functions.generate_wvt(binlist,varlist[i])
            main.saveiteratedfits(targlist[m],binlist,wscxlist[i],wvt,vwvt,objlist[i],sourcelist[i],subfolder=subfolder)
            main.saveston(wscxlist[i],ston,sourcelist[i],objlist[i],subfolder=subfolder)
            assign=functions.assign(binlist,siglist[i])
            main.saveassign(wscxlist[i],assign,sourcelist[i],objlist[i],subfolder=subfolder)

print("Bye Bye!")