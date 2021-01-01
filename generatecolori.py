
import os

def findfilewith(directory,xxx,yyy):
    files=os.listdir(directory)
    for m in files:
        if xxx in m and yyy in m:
            print(m)
            return m
    raise NameError

def colorat(directory):
    try:
        os.mkdir(directory+"/colorimages")
        f=open("../../Downloads/"+subfolder+"/colorimages/prescript.sh",'w+')
        f.write("alias ds9='/Applications/ds9.darwinsierra.8.1/ds9'")
    except:
        print("exists")

nebulae=['j0248-0817','j0823+0313','j1044+0353','j1238+1009']
nebnames=['j0248','j0823','j1044','j1238']

targs=[0,3,5,10]


for n in range(len(nebulae)):
    for t in range(len(targs)):
        if targs[t]==0:
            f=open("../../Downloads/binned_Sep29/colorimages/"+str(nebnames[n])+"_ub_r4862_g4686_b4472.sh",'w+')
            f.writelines(". prescript.sh\n")
            r="../"+nebulae[n]+"/unbinned/"+findfilewith("../../Downloads/binned_Sep29/"+nebulae[n]+"/unbinned","486","icubes")
            g="../"+nebulae[n]+"/unbinned/"+findfilewith("../../Downloads/binned_Sep29/"+nebulae[n]+"/unbinned","468","icubes")
            b="../"+nebulae[n]+"/unbinned/"+findfilewith("../../Downloads/binned_Sep29/"+nebulae[n]+"/unbinned","447","icubes")
            f.writelines("ds9 -frame delete -rgb -red "+r+" -scale log -scalelimits 0 2 -green "+g+" -scale log -scalelimits 0 2 -blue "+b+" -scale log -scalelimits 0 2\n")
        else:
            f=open("../../Downloads/binned_Sep29/colorimages/"+nebnames[n]+"_t"+str(targs[t])+"_r4862_g4686_b4472.sh",'w+')
            f.write(". prescript.sh\n")
            r="../"+nebulae[n]+"/target"+str(targs[t])+"/"+findfilewith("../../Downloads/binned_Sep29/"+nebulae[n]+"/target"+str(targs[t]),"486","block")
            g="../"+nebulae[n]+"/target"+str(targs[t])+"/"+findfilewith("../../Downloads/binned_Sep29/"+nebulae[n]+"/target"+str(targs[t]),"468","block")
            b="../"+nebulae[n]+"/target"+str(targs[t])+"/"+findfilewith("../../Downloads/binned_Sep29/"+nebulae[n]+"/target"+str(targs[t]),"447","block")
            f.write("ds9 -frame delete -rgb -red "+r+" -scale log -scalelimits 0 2 -green "+g+" -scale log -scalelimits 0 2 -blue "+b+" -scale log -scalelimits 0 2\n")

for n in range(len(nebulae)):
    f=open("../../Downloads/binned_Sep29/colorimages/"+str(nebnames[n])+"_all_r4862_g4686_b4472.sh",'w+')
    f.writelines(". prescript.sh\n")
    for t in range(len(targs)):
        if targs[t]==0:
            r="../"+nebulae[n]+"/unbinned/"+findfilewith("../../Downloads/binned_Sep29/"+nebulae[n]+"/unbinned","486","icubes")
            g="../"+nebulae[n]+"/unbinned/"+findfilewith("../../Downloads/binned_Sep29/"+nebulae[n]+"/unbinned","468","icubes")
            b="../"+nebulae[n]+"/unbinned/"+findfilewith("../../Downloads/binned_Sep29/"+nebulae[n]+"/unbinned","447","icubes")
            f.writelines("ds9 -frame delete -frame lock image -rgb -red "+r+" -scale log -scalelimits 0 2 -green "+g+" -scale log -scalelimits 0 2 -blue "+b+" -scale log -scalelimits 0 2 ")
        else:
            r="../"+nebulae[n]+"/target"+str(targs[t])+"/"+findfilewith("../../Downloads/binned_Sep29/"+nebulae[n]+"/target"+str(targs[t]),"486","block")
            g="../"+nebulae[n]+"/target"+str(targs[t])+"/"+findfilewith("../../Downloads/binned_Sep29/"+nebulae[n]+"/target"+str(targs[t]),"468","block")
            b="../"+nebulae[n]+"/target"+str(targs[t])+"/"+findfilewith("../../Downloads/binned_Sep29/"+nebulae[n]+"/target"+str(targs[t]),"447","block")
            f.write("-rgb -red "+r+" -scale log -scalelimits 0 2 -green "+g+" -scale log -scalelimits 0 2 -blue "+b+" -scale log -scalelimits 0 2 ")


'''
nebulae=['exp10000']
nebnames=nebulae

subfolder="Nov6/beta0.9"
colorat("../../Downloads/"+subfolder)

targs=[0,3,5,10]
edges=["edge50","edge75","edge100"]


for n in range(len(nebulae)):
    for t in range(len(targs)):
        for e in range(len(edges)):
            if targs[t]==0:
                f=open("../../Downloads/Sep4/128x128_peak30/colorimages/"+str(nebnames[n])+"_ub_"+edges[e]+".sh",'w+')
                f.writelines(". prescript.sh\n")
                r="../"+nebulae[n]+"/unbinned/"+findfilewith("../../Downloads/Sep4/128x128_peak30/"+nebulae[n]+"/unbinned","testdata_",edges[e])
                f.writelines("ds9 -frame delete -fits "+r+" -scale log -scalelimits 0 2\n")
            else:
                f=open("../../Downloads/Sep4/128x128_peak30/colorimages/"+nebnames[n]+"_t"+str(targs[t])+"_"+edges[e]+".sh",'w+')
                f.write(". prescript.sh\n")
                r="../"+nebulae[n]+"/target"+str(targs[t])+"/"+findfilewith("../../Downloads/Sep4/128x128_peak30/"+nebulae[n]+"/target"+str(targs[t]),"testdata_",edges[e])
                f.writelines("ds9 -frame delete -fits "+r+" -scale log -scalelimits 0 2\n")

for n in range(len(nebulae)):
    for e in range(len(edges)):
        f=open("../../Downloads/"+subfolder+"/colorimages/"+str(nebnames[n])+"_all_"+edges[e]+".sh",'w+')
        f.writelines(". prescript.sh\n")
        for t in range(len(targs)):
            if targs[t]==0:
                r="../"+nebulae[n]+"/unbinned/"+findfilewith("../../Downloads/"+subfolder+"/"+nebulae[n]+"/unbinned","testdata_",edges[e])
                f.writelines("ds9 -frame delete -frame lock image -fits "+r+" -scale log -scalelimits 0 40 ")
            else:
                r="../"+nebulae[n]+"/target"+str(targs[t])+"/"+findfilewith("../../Downloads/"+subfolder+"/"+nebulae[n]+"/target"+str(targs[t]),"testdata_",edges[e])
                f.write("-fits "+r+" -scale log -scalelimits 0 40 ")
        f.write("-single")

for e in range(len(edges)):
    for t in range(len(targs)):
        if targs[t]==0:
            f=open("../../Downloads/"+subfolder+"/colorimages/"+edges[e]+"_all_ub.sh",'w+')
        else:
            f=open("../../Downloads/"+subfolder+"/colorimages/"+edges[e]+"_all_t"+str(targs[t])+".sh",'w+')
        f.writelines(". prescript.sh\n")
        for n in range(len(nebulae)):
            if targs[t]==0:
                r="../"+nebulae[n]+"/unbinned/"+findfilewith("../../Downloads/"+subfolder+"/"+nebulae[n]+"/unbinned","testdata_",edges[e])
            else:
                r="../"+nebulae[n]+"/target"+str(targs[t])+"/"+findfilewith("../../Downloads/"+subfolder+"/"+nebulae[n]+"/target"+str(targs[t]),"testdata_",edges[e])

            if n==0:
                f.writelines("ds9 -frame delete -frame lock image -fits "+r+" -scale log -scalelimits 0 40 ")
            else:
                f.write("-fits "+r+" -scale log -scalelimits 0 40 ")
        f.write("-single")
'''