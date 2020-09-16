import numpy as np
import matplotlib.pyplot as plt 
import tkinter
from tkinter.filedialog import askopenfilename
from astropy.io import fits
import bin_accretion,wvt_iteration,functions
from matplotlib.animation import FuncAnimation,PillowWriter
import time

def iteration_func(target,signal,var,geocarray,scalearray,displaywvt=False):
    start=time.time()
    ## have to manually kill terminal is does not converge
    supremebinlist=[]
    supremevibelist=[]
    freaky=[] 
    
    print("another iteration")
    geocarray2=np.copy(geocarray)
    scalearray2=np.copy(scalearray)
    wvt2=np.copy(wvt)
    assign=np.full_like(signal,-1, dtype=int)
    viable=[]
    for g in range(len(geocarray)):
        point=(int(geocarray[g][0]),int(geocarray[g][1]))
        try:
            assign[point[0]][point[1]]=g
        except:
            print(point)
            print(geocarray[g])
            print(g)
            raise NameError("ouchi")
        viable.append([])
        wvt_iteration.append_validate((point[0]+1,point[1]),viable[g],assign)
        wvt_iteration.append_validate((point[0]-1,point[1]),viable[g],assign)
        wvt_iteration.append_validate((point[0],point[1]+1),viable[g],assign)
        wvt_iteration.append_validate((point[0],point[1]-1),viable[g],assign)
        if len(viable[g])==0:
            if(geocarray[g][0]<0 or geocarray[g][1]<0):
                print("Disapprove but like. ok")
                freaky.append(geocarray[g])
            else:
                raise NameError("uh ruh roh")
        #print(str(int(g*100/len(geocarray)))+" percent done with init pass")
    while wvt_iteration.checkneg(assign) or wvt_iteration.viabempty(viable):
        for g in range(len(geocarray)):
            prune=True
            while prune and len(viable[g])>0:
                point=viable[g].pop(0)
                if assign[point[0]][point[1]]==g:
                    prune=True
                else:
                    prune=False
            if len(viable[g])>0:
                if assign[point[0]][point[1]]==-1:
                    assign[point[0]][point[1]]=g
                    wvt_iteration.append_validate((point[0]+1,point[1]),viable[g],assign)
                    wvt_iteration.append_validate((point[0]-1,point[1]),viable[g],assign)
                    wvt_iteration.append_validate((point[0],point[1]+1),viable[g],assign)
                    wvt_iteration.append_validate((point[0],point[1]-1),viable[g],assign)
                else:
                    if ((geocarray[g][0]-point[0])**2+(geocarray[g][1]-point[1])**2)/(scalearray[g]**2)<((geocarray[assign[point[0]][point[1]]][0]-point[0])**2+(geocarray[assign[point[0]][point[1]]][1]-point[1])**2)/(scalearray[assign[point[0]][point[1]]]**2):
                        assign[point[0]][point[1]]=g
                        wvt_iteration.append_validate((point[0]+1,point[1]),viable[g],assign)
                        wvt_iteration.append_validate((point[0]-1,point[1]),viable[g],assign)
                        wvt_iteration.append_validate((point[0],point[1]+1),viable[g],assign)
                        wvt_iteration.append_validate((point[0],point[1]-1),viable[g],assign)
        
    binlist=[ [] for _ in range(len(geocarray)) ]
    for j in range(len(assign)):
        for i in range(len(assign[0])):
            binlist[assign[j][i]].append((j,i))
    geocarray2,scalearray2=functions.calculate_scales(target,binlist,signal,var)
    for r in range(len(binlist)):
        if len(binlist[r])==0:
            print("empty index"+str(r))
            geocarray2[r]=geocarray[r]
            scalearray2[r]=scalearray[r]
        
    
    print("elapsed time "+str(time.time()-start))

    return geocarray2,scalearray2

def iteration_funca(target,signal,var,geocarray,scalearray,displaywvt=False):
    start=time.time()
    ## have to manually kill terminal is does not converge
    supremebinlist=[]
    supremevibelist=[]
    freaky=[] 
    contin=True
    print("another iteration")
    geocarray2=np.copy(geocarray)
    scalearray2=np.copy(scalearray)
    wvt2=np.copy(wvt)
    assign=np.full_like(signal,-1, dtype=int)
    viable=[]
    for g in range(len(geocarray)):
        point=(int(geocarray[g][0]),int(geocarray[g][1]))
        try:
            assign[point[0]][point[1]]=g
        except:
            print(point)
            print(geocarray[g])
            print(g)
            raise NameError("ouchi")
        viable.append([])
        wvt_iteration.append_validate((point[0]+1,point[1]),viable[g],assign)
        wvt_iteration.append_validate((point[0]-1,point[1]),viable[g],assign)
        wvt_iteration.append_validate((point[0],point[1]+1),viable[g],assign)
        wvt_iteration.append_validate((point[0],point[1]-1),viable[g],assign)
        if len(viable[g])==0:
            if(geocarray[g][0]<0 or geocarray[g][1]<0):
                print("Disapprove but like. ok")
                freaky.append(geocarray[g])
            else:
                raise NameError("uh ruh roh")
        #print(str(int(g*100/len(geocarray)))+" percent done with init pass")
    while wvt_iteration.checkneg(assign) or wvt_iteration.viabempty(viable):
        for g in range(len(geocarray)):
            prune=True
            while prune and len(viable[g])>0:
                point=viable[g].pop(0)
                if assign[point[0]][point[1]]==g:
                    prune=True
                else:
                    prune=False
            if len(viable[g])>0:
                if assign[point[0]][point[1]]==-1:
                    assign[point[0]][point[1]]=g
                    wvt_iteration.append_validate((point[0]+1,point[1]),viable[g],assign)
                    wvt_iteration.append_validate((point[0]-1,point[1]),viable[g],assign)
                    wvt_iteration.append_validate((point[0],point[1]+1),viable[g],assign)
                    wvt_iteration.append_validate((point[0],point[1]-1),viable[g],assign)
                else:
                    try:
                        if ((geocarray[g][0]-point[0])**2+(geocarray[g][1]-point[1])**2)/(scalearray[g]**2)<((geocarray[assign[point[0]][point[1]]][0]-point[0])**2+(geocarray[assign[point[0]][point[1]]][1]-point[1])**2)/(scalearray[assign[point[0]][point[1]]]**2):
                            assign[point[0]][point[1]]=g
                            wvt_iteration.append_validate((point[0]+1,point[1]),viable[g],assign)
                            wvt_iteration.append_validate((point[0]-1,point[1]),viable[g],assign)
                            wvt_iteration.append_validate((point[0],point[1]+1),viable[g],assign)
                            wvt_iteration.append_validate((point[0],point[1]-1),viable[g],assign)
                    except:
                        print(point)
                        print(assign[point[0]][point[1]])
                        print(geocarray[assign[point[0]][point[1]]])
                        print(geocarray[g])
                        raise NameError("ooch")
        if contin:
            vibearray=np.zeros_like(signal)
            binnarray=np.zeros_like(signal)
            for thing in range(len(viable)):
                for point in viable[thing]:
                    vibearray[point[0]][point[1]]=100
            for j in range(len(assign)):
                for i in range(len(assign)):
                    if assign[j][i]!=-1:
                        binnarray[j][i]=assign[j][i]
            supremebinlist.append(binnarray)
            supremevibelist.append(np.ma.masked_where(vibearray==0,vibearray))
            if len(supremevibelist)>(100+offset):
                contin=False
        #update2(vibearray,binnarray)
        #fig.canvas.draw()
        
    binlist=[ [] for _ in range(len(geocarray)) ]
    for j in range(len(assign)):
        for i in range(len(assign[0])):
            binlist[assign[j][i]].append((j,i))
    geocarray2,scalearray2=functions.calculate_scales(target,binlist,signal,var)
    for r in range(len(binlist)):
        if len(binlist[r])==0:
            print("empty index"+str(r))
            geocarray2[r]=geocarray[r]
            scalearray2[r]=scalearray[r]
        
    
    print("elapsed time "+str(time.time()-start))

    return supremebinlist,supremevibelist,np.array(freaky)

def update(num,contourb,contourf):
    ax1.clear
    g=ax1.imshow(contourb[num],cmap="plasma")
    ax1.imshow(contourf[num],cmap="binary")
    print(str(num)+" out of "+str(len(contourb)))

def update2(contourb,contourf):
    ax1.clear
    g=ax1.imshow(contourb,cmap="plasma")
    ax1.imshow(contourf,cmap="binary")

wscx,signal,var,sourcedir,objname=bin_accretion.initialize(enternew=True)

off=[800]


#plt.ion()
fig,ax1=plt.subplots()

target=5
epsilon=0

binlist,geocarray=bin_accretion.cc_accretion(signal,var,target)

print(geocarray)

scalearray=np.full(len(geocarray),1)
#wvt=functions.generate_wvt(binlist,signal)

for i in range(epsilon):
    geocarray,scalearray=iteration_func(target,signal,var,geocarray,scalearray)

np.random.shuffle(geocarray)
for offset in off:
    sbl,svl,freakya=iteration_funca(target,signal,var,geocarray,scalearray)



    print("iteration over")

    anim=FuncAnimation(fig,update,frames=range(offset,len(sbl)),fargs=(sbl,svl),interval=50)
    print("setup animation")
    objname="s_testdata"

    anim.save(sourcedir+"/"+objname+"_"+str(offset)+"_vid.gif",writer=PillowWriter(fps=24))

print("exit")