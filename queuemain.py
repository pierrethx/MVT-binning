import bin_accretion,main

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
            objname=main.getname(objname)
            target=main.gettarget(targhold)
            sourcelist.append(sourcedir)
            wscxlist.append(wcsx)
            siglist.append(signal)
            varlist.append(var)
            objlist.append(objname)
            targlist.append(target)
            targhold=target
        except:
            contqueue=False
    print("Files loaded!")
    for i in range(len(sourcelist)):
        if targlist[i]<0:
            weighting=False
            targlist[i]=-targlist[i]
        else:
            weighting=True
        wvt,vwvt=main.mainfunc(siglist[i],varlist[i],targlist[i],weighting=weighting,displayWVT=False,epsilon=-10)
        main.saveiteratedfits(targlist[i],wscxlist[i],wvt,vwvt,objlist[i],sourcelist[i],weighting=weighting)