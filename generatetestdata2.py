import os
import numpy as np
from astropy.io import fits

def makeitfolder(destination,subfolder,text,expt,roundnum=2):
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
        targs=["unbinned","target3","target5","target8","target10","target15"]
        for t in targs:
            try:
                os.mkdir(destination+subfolder+"/"+text+str(round(ex,roundnum))  +"/"+t)
            except:
                pass
                print(destination+subfolder+"/"+text+str(round(ex,roundnum))  +"/"+t," already exists")
    return [destination+subfolder+"/"+text+str(round(ex,roundnum)) for ex in expt]


def generator(widpix,heipix,A,r,b,I_b,edge,frames):
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
    noise=np.sqrt((model+I_b)/frames)*(2*np.random.rand(heipix,widpix)-1)

    var=(model+noise+I_b)/frames
    signal=model+noise

    return signal,var,model

def gendata(num,wp,hp,frames,cases,destination,subfolder,objname,subtitle,filetitle,model=False):
    rn=2 #number to round decimal place to
    makeitfolder(destination,subfolder,subtitle,cases[subtitle],rn)
    for As in cases["sig"]:
        for rs in cases["rcore"]:
            for be in cases["beta"]:
                for bg in cases["bg"]:
                    for edge in cases["edge"]:
                        temp = { "sig": As, "rcore": rs, "beta":be, "bg":bg, "edge":edge}
                        for n in range(num):
                            sig,var,model=generator(wp,hp,As,rs,be,bg,edge,frames)
                            hdu = fits.PrimaryHDU(sig)
                            hdul = fits.HDUList([hdu])
                            hdul.writeto(destination+subfolder+"/"+subtitle+str(round(temp[subtitle],rn))+"/unbinned/"+objname+str(n)+"_"+filetitle+str(round(temp[filetitle],rn))+"_signal.fits",overwrite=True)
                            hdu = fits.PrimaryHDU(var)
                            hdul = fits.HDUList([hdu])
                            hdul.writeto(destination+subfolder+"/"+subtitle+str(round(temp[subtitle],rn))+"/unbinned/"+objname+str(n)+"_"+filetitle+str(round(temp[filetitle],rn))+"_var.fits",overwrite=True)

def genmodel(num,wp,hp,frames,cases,destination,subfolder,objname,subtitle,filetitle):
    rn=2 #number to round decimal place to
    makeitfolder(destination,subfolder,subtitle,cases[subtitle],rn)
    for As in cases["sig"]:
        for rs in cases["rcore"]:
            for be in cases["beta"]:
                for bg in cases["bg"]:
                    for edge in cases["edge"]:
                        temp = { "sig": As, "rcore": rs, "beta":be, "bg":bg, "edge":edge}
                        sig,var,model=generator(wp,hp,As,rs,be,bg,edge,frames)
                        hdu0 = fits.PrimaryHDU(model)
                        hdul0 = fits.HDUList([hdu0])
                        hdul0.writeto(destination+subfolder+"/"+subtitle+str(round(temp[subtitle],rn))+"/unbinned/"+objname+"_"+filetitle+str(round(temp[filetitle],rn))+"_model.fits",overwrite=True)

def gendataAB(num,wp,hp,frames,cases,destination,subfolder,objname,subtitle,filetitle,model=False):
    rn=2 #number to round decimal place to
    makeitfolder(destination,subfolder,subtitle,cases[subtitle],rn)
    for i in range(len(cases["sig"])):
        for rs in cases["rcore"]:
            for be in cases["beta"]:
                for edge in cases["edge"]:
                    temp = { "sig": cases["sig"][i], "rcore": rs, "beta":be, "bg":cases["bg"][i], "edge":edge}
                    for n in range(num):
                        sig,var,model=generator(wp,hp,temp["sig"],rs,be,temp["bg"],edge,frames)
                        hdu = fits.PrimaryHDU(sig)
                        hdul = fits.HDUList([hdu])
                        hdu2 = fits.PrimaryHDU(var)
                        hdul2 = fits.HDUList([hdu2])
                        hdul.writeto(destination+subfolder+"/"+subtitle+str(round(temp[subtitle],rn))+"/unbinned/"+objname+str(n)+"_"+filetitle+str(round(temp[filetitle],rn))+"_signal.fits",overwrite=True)
                        hdul2.writeto(destination+subfolder+"/"+subtitle+str(round(temp[subtitle],rn))+"/unbinned/"+objname+str(n)+"_"+filetitle+str(round(temp[filetitle],rn))+"_var.fits",overwrite=True)
                    hdu0 = fits.PrimaryHDU(model)
                    hdul0 = fits.HDUList([hdu0])
                    hdul0.writeto(destination+subfolder+"/"+subtitle+str(round(temp[subtitle],rn))+"/unbinned/"+objname+"_"+filetitle+str(round(temp[filetitle],rn))+"_model.fits",overwrite=True)
                          

if __name__ == "__main__":
    '''
    A=[1,10,100,1000,10000]
    r=[16]
    b=[0.67]
    I_b=[0.01*i for i in A]
    edge=[75]
    frames=10000
    widpix=160
    heipix=160
    destination="/Volumes/TOSHIBA EXT/MARTIN2/"
    subfolder="simdata/123test"
    cases = { "sig": A, "rcore": r, "beta":b, "bg":I_b, "edge":edge}
    gendataAB(10,widpix,heipix,frames,cases,destination,subfolder,"testdata","sig","edge",model=True)
    '''
    A=[10]
    r=[4,8,12,16,20]
    b=[0.67,1.17,1.67,2.17,2.67,3.17]
    I_b=[0.01]
    edge=[75]
    frames=10000

    widpix=160
    heipix=160

    destination="/Volumes/TOSHIBA EXT/MARTIN2/"
  
    subfolder="simdata/betatest"
    cases = { "sig": A, "rcore": r, "beta":b, "bg":I_b, "edge":edge}
    #gendata(1,widpix,heipix,frames,cases,destination,subfolder,"testdata","beta","rcore",model=True)
    genmodel(1,widpix,heipix,frames,cases,destination,subfolder,"testdata","beta","rcore")
    
    