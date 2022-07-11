import os
import numpy as np
from astropy.io import fits

def makefolder(destination,subfolder,text,expt,roundnum=2):
    slist=subfolder.split("/")
    for s in range(len(slist)):
        try:
            os.mkdir(destination+"/".join(slist[:s+1]))
        except:
            pass
    for ex in expt:
        try:
            os.mkdir(destination+subfolder+"/"+text+str(round(ex,roundnum)))
        except:
            pass
        try:
            os.mkdir(destination+subfolder+"/"+text+str(round(ex,roundnum))+"/unbinned")
        except:
            pass
    return [destination+subfolder+"/"+text+str(round(ex,roundnum)) for ex in expt]

def profile(widpix,heipix,A,r,b,edge):
    xcent=int(widpix/2)
    ycent=int(heipix/2)

    xlist=np.arange(0,widpix,1)
    ylist=np.arange(0,heipix,1)
    xx,yy=np.meshgrid(xlist,ylist)

    model=A*((1+((xx-xcent)**2+(yy-ycent)**2)/r**2))**(0.5-3*b)
    if np.isnan(edge):
        pass
    else:
        for j in range(len(xx)):
            for i in range(len(xx[0])):
                if ((xx[j][i]-xcent)**2+(yy[j][i]-ycent)**2)>edge**2:
                    model[j][i]=0
    return model


def generator(widpix,heipix,func,arguments,I_b,frames):

    model=func(widpix,heipix,*arguments)
    
    noise=np.sqrt((model+I_b)/frames)*(2*np.random.rand(heipix,widpix)-1)

    var=(model+noise+I_b)/frames
    signal=model+noise

    return signal,var,model


def genprof(num,wp,hp,frames,cases,destination,subfolder,objname,subtitle,filetitle,mode=True,pairings=None):
    ## this is a function for generating Circular Beta Profiles (Sarazin 1988)
    rn=2 #number to round decimal place to
    makefolder(destination,subfolder,subtitle,cases[subtitle],rn)

    keyo=[]
    keyp=[]
    keyt=[]
    for key in cases:
        if pairings==None or not any([key in p for p in pairings]):
            keyo.append(key)
            keyp.append(-1)
            keyt.append(len(cases[key]))
        else:
            for p in range(len(pairings)):
                if key==pairings[p][0]:
                    keyo.append(key)
                    keyp.append(p)
                    keyt.append(len(cases[key]))
                    break
    for i in range(np.product(keyt)):
        temp={"sig":-1,"rcore":-1,"beta":-1,"bg":-1,"edge":-1}
        t=np.unravel_index(i,tuple(keyt))
        for i in range(len(t)):
            if keyp[i]==-1:
                temp[keyo[i]]=cases[keyo[i]][t[i]]
            else:
                for title in pairings[keyp[i]]:
                    temp[title]=cases[title][t[i]]   
        for n in range(num):
            sig,var,model=generator(wp,hp,profile,(temp["sig"],temp["rcore"],temp["beta"],temp["edge"]),temp["bg"],frames)
            hdu = fits.PrimaryHDU(sig)
            hdul = fits.HDUList([hdu])
            hdul.writeto(destination+subfolder+"/"+subtitle+str(round(temp[subtitle],rn))+"/unbinned/"+objname+str(n)+"_"+filetitle+str(round(temp[filetitle],rn))+"_sig.fits",overwrite=True)
            hdu = fits.PrimaryHDU(var)
            hdul = fits.HDUList([hdu])
            hdul.writeto(destination+subfolder+"/"+subtitle+str(round(temp[subtitle],rn))+"/unbinned/"+objname+str(n)+"_"+filetitle+str(round(temp[filetitle],rn))+"_var.fits",overwrite=True)
        if mode:
            hdu0 = fits.PrimaryHDU(model)
            hdul0 = fits.HDUList([hdu0])
            hdul0.writeto(destination+subfolder+"/"+subtitle+str(round(temp[subtitle],rn))+"/unbinned/"+objname+"_"+filetitle+str(round(temp[filetitle],rn))+"_model.fits",overwrite=True)

def grad(widpix,heipix,left,right):
    xlist=np.arange(0,widpix,1)
    ylist=np.arange(0,heipix,1)
    xx,yy=np.meshgrid(xlist,ylist)
    lazx=widpix-1
    lgrad=right*xx/lazx+left*(1-xx/lazx)
    return lgrad

def hgrad(widpix,heipix,right,power):
    xlist=np.arange(0,widpix,1)
    ylist=np.arange(0,heipix,1)
    xx,yy=np.meshgrid(xlist,ylist)
    lazx=widpix-1
    pgrad=right*((2*xx/lazx-1)**power)
    pgrad[xx<(widpix)/2]=0
    return pgrad

def island(widpix,heipix,island,ocean,radius):
    xlist=np.arange(0,widpix,1)
    ylist=np.arange(0,heipix,1)
    xx,yy=np.meshgrid(xlist,ylist)
    xcent=int(widpix/2.0)
    ycent=int(heipix/2.0)
    ocean=0*xx+ocean
    island=0*xx+island
    
    rrsq=(xx-xcent)**2+(yy-ycent)**2
    for y in range(len(rrsq)):
        for x in range(len(rrsq[0])):
            if rrsq[y][x]<=radius**2:
                rrsq[y][x]=island[y][x]
            else:
                rrsq[y][x]=ocean[y][x]
    return rrsq

def gentestcases(num,wp,hp,frames,params,destination,subfolder,objname,mode=True):
    
    #params of form [A,r,bg]
    A,r,bg=params
    ## this is a function for generating Circular Beta Profiles (Sarazin 1988)
    rn=2 #number to round decimal place to
    makefolder(destination,subfolder,"sig",[A],rn)

    tests={"constp":{"func":grad,"params":[A,A],},
    "const0":{"func":grad,"params":[0,0],},
    "constn":{"func":grad,"params":[-A,-A],},
    "lgradpp":{"func":grad,"params":[A/10,A],},
    "lgradpn":{"func":grad,"params":[A,-A],},
    "islandp0":{"func":island,"params":[A,0,r],},
    "islandppn":{"func":island,"params":[A,-0.5*A,r],},
    "islandpnn":{"func":island,"params":[0.5*A,-A,r],},
    "island0p":{"func":island,"params":[0,A,r],},
    "islandnpp":{"func":island,"params":[-A,2*A,r],},
    "islandnnp":{"func":island,"params":[-A,0.5*A,r],},
    "lgradppn":{"func":grad,"params":[-0.3*A,A],},
    "lgradnnp":{"func":grad,"params":[-A,0.3*A],},
    "hLgrad":{"func":hgrad,"params":[A,1],},
    "h2grad":{"func":hgrad,"params":[A,2],},
    "h3grad":{"func":hgrad,"params":[A,3],},}

    for key in tests:
        for n in range(num):
            sig,var,model=generator(wp,hp,tests[key]["func"],tests[key]["params"],bg,frames)
            hdu = fits.PrimaryHDU(sig)
            hdul = fits.HDUList([hdu])
            hdul.writeto(destination+subfolder+"/"+"sig"+str(round(A,rn))+"/unbinned/"+objname+str(n)+"_"+key+"_signal.fits",overwrite=True)
            hdu = fits.PrimaryHDU(var)
            hdul = fits.HDUList([hdu])
            hdul.writeto(destination+subfolder+"/"+"sig"+str(round(A,rn))+"/unbinned/"+objname+str(n)+"_"+key+"_var.fits",overwrite=True)
        if mode:
            hdu0 = fits.PrimaryHDU(model)
            hdul0 = fits.HDUList([hdu0])
            hdul0.writeto(destination+subfolder+"/"+"sig"+str(round(A,rn))+"/unbinned/"+objname+"_"+key+"_model.fits",overwrite=True)


if __name__ == "__main__":
    frames=10000
    widpix=160
    heipix=160
    destination="/Volumes/TOSHIBA EXT/MARTIN2/tupletest/"
    
    cases = { "sig": [1,10,100,1000,10000], "rcore": [16], "beta":[0.67], "bg":[0.01,0.1,1,10,100], "edge":[75]}
    genprof(10,widpix,heipix,frames,cases,destination,"123test","testdata","sig","edge",pairings=[("sig","bg")])
    cases = { "sig": [10], "rcore": [4,8,12,16,20], "beta":[0.67,1.17,1.67,2.17,2.67,3.17], "bg":[0.01], "edge":[75]}
    genprof(1,widpix,heipix,frames,cases,destination,"betatest","testdata","beta","rcore")

    gentestcases(1,widpix,heipix,frames,[5,48,10],destination,"testcases","testcase",mode=False)