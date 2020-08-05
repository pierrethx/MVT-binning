import numpy as np
import matplotlib.pyplot as plt 
import tkinter
from tkinter.filedialog import askopenfilename
from astropy.io import fits
import bin_accretion,wvt_iteration,functions
from matplotlib.animation import FuncAnimation,PillowWriter
import time



def iteration_funca(target,signal,var,geocarray,scalearray,epsilon,displaywvt=False):
    target=5
    start=time.time()
    ## have to manually kill terminal is does not converge
    difference=2*epsilon
    supremebinlist=[]
    supremevibelist=[]
    freaky=[]
    

    while difference>epsilon:
        print("another iteration")
        geocarray2=np.copy(geocarray)
        scalearray2=np.copy(scalearray)
        wvt2=np.copy(wvt)
        assign=np.full_like(signal,-1)
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
            if len(supremebinlist)>100:
                break
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
        difference=0
        
    
    print("elapsed time "+str(time.time()-start))

    return supremebinlist,supremevibelist,np.array(freaky)

def update(num,contourb,contourf):
    ax1.clear
    g=ax1.imshow(contourb[num],cmap="plasma")
    ax1.imshow(contourf[num],cmap="binary")
    print(str(num)+" out of "+str(len(contourb)))

wscx,signal,var,sourcedir,objname=bin_accretion.initialize(enternew=True)


target=20
epsilon=4000
binlist,geocarray=bin_accretion.cc_accretion(signal,var,target)

print(geocarray)

scalearray=np.full(len(geocarray),1)
wvt=functions.generate_wvt(binlist,signal)



#np.random.shuffle(geocarray)
sbl,svl,freakya=iteration_funca(target,signal,var,geocarray,scalearray,epsilon)

fig,ax=plt.subplots()

ax.imshow(wvt,cmap="cubehelix")
for binn in binlist:
    bim=np.array(binn)
    ax.plot(bim[:,1],bim[:,0])
plt.show()

raise NameError(":0")

print("iteration over")
fig,ax1=plt.subplots()

anim=FuncAnimation(fig,update,frames=range(len(sbl)),fargs=(sbl,svl),interval=50)
print("setup animation")
objname="testdata"
anim.save(sourcedir+"/"+objname+"_vid.gif",writer=PillowWriter(fps=24))
print("exit")