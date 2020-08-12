import bin_accretion,main,qradial,queuemain

if __name__ == "__main__":
    superwcsxlist=[]
    supersiglist=[]
    supervarlist=[]
    supersourcelist=[]
    superobjlist=[]

    wcsxlist=[]
    siglist=[]
    varlist=[]
    sourcelist=[]    
    objlist=[]
    targlist=[]
    contqueue=True
    targhold=0
    
    while contqueue:
        try:
            wcsx,signal,var,sourcedir,objname=bin_accretion.initialize(enternew=True)
            objname=main.getname("_".join(objname.split("_")[:-1]))
            sourcelist.append(sourcedir)
            wcsxlist.append(wcsx)
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
            wvt,vwvt=main.mainfunc(siglist[i],varlist[i],targlist[m],displayWVT=False,epsilon=-10)
            a,b,c,d,e=main.saveiteratedfits(targlist[m],wcsxlist[i],wvt,vwvt,objlist[i],sourcelist[i],weighting=weighting,subfolder=subfolder)
            superwcsxlist.append(a)
            supersiglist.append(b)
            supervarlist.append(c)
            supersourcelist.append(d+"/"+subfolder)
            superobjlist.append(e)
    qradial.qradial(superwcsxlist,supersiglist,supervarlist,supersourcelist,superobjlist)
print("Bye Bye!")