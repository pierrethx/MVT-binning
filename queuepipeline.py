import bin_accretion,main,qradial,queuemain,generatetestdata,radial_profile,functions
import numpy as np
import matplotlib.pyplot as plt

def tabulate(pathfile,slist,obj,first,sim,maxi,nedges,conc):
    f=open(pathfile,"w+")
    f.write("\\begin{tab"+"ular}{|c||c|c|c|}\n \\hline\n")
    f.write("file & recovered edge & simulated edge & sim maximum edge & num nan edges& conclusion\\\\ \n \\hline \n")
    
    for i in range(len(nedges)):
        f.write(slist[i].split("/")[-1]+"/"+obj[i]+" & "+first[i]+" & "+sim[i]+" & "+maxi[i]+" & "+nedges[i]+" & "+conc[i]+" \\\\ \n")
    f.write("\\hline\n\\end{tab"+"ular}")

if __name__ == "__main__":

    name=main.getname("_.txt")
    
    pathfile="/Users/pierre/Downloads/"+name

    wcsxlist,siglist,varlist,sourcelist,objlist=bin_accretion.minitialize()
    outputlist=[]


    data=qradial.qradial(wcsxlist,siglist,varlist,sourcelist,objlist,outputlist)
    """[binning,seededge,firstpassedge,exptime,foundedge,stonedgem,stonedgea]"""
    # outputlist is a list of popts, so [popt1,popt2,...]

    print(outputlist)

    firstlist=[]
    maxlist=[]
    simlist=[]
    numnanedges=[]
    conclist=[]
    numtimes=5

    for fil in range(len(data)):
        first=data[fil][6]
        firstlist.append(str(round(first,3)))
        if np.isnan(first) or first<1:
            maxlist.append("nan")
            simlist.append("nan")
            numnanedges.append("nan")
            conclist.append("no edge")
        else:
            edges=[]
            edges2=[]
            nedge=0
            for too in range(numtimes):
                expt=np.nanmax(siglist[fil]/varlist[fil])
                Bg=np.abs(np.average(expt*varlist[fil]-siglist[fil],weights=np.ma.masked_where(np.isnan(siglist[fil]/np.sqrt(varlist[fil])),siglist[fil]/np.sqrt(varlist[fil]))))

                popt=outputlist[fil]
                print(first)
                center=radial_profile.getcenter(siglist[fil])
                factor=1.4
                si,va=generatetestdata.ogenerator(int(factor*len(siglist[fil])),int(factor*len(siglist[fil][0])),(factor*center[0],factor*center[1]),np.abs(popt[0]),popt[1],popt[2],Bg,expt,np.nan)

                if data[fil][0]==0:
                    wvt=si
                    vwvt=va
                else:
                    binlist=main.mainfunc(si,va,data[fil][0],displayWVT=False,epsilon=-10)
                    wvt,ston=functions.generate_wvt3(binlist,si,va,np.full(len(binlist),1))
                    vwvt=functions.generate_wvt(binlist,va)
                edgex,edge2,edgec,edgea=radial_profile.radmethod(wvt,vwvt,wcsxlist[fil],show=False)
                #these are first pass, fit, ston median, and ston average. we want the last one.
                edges.append(edgea)
                if np.isnan(edgea):
                    nedge+=1

                si,va=generatetestdata.ogenerator(len(siglist[fil]),len(siglist[fil][0]),center,popt[0],popt[1],popt[2],Bg,expt,first)
                if data[fil][0]==0:
                    wvt=si
                    vwvt=va
                else:
                    binlist=main.mainfunc(si,va,data[fil][0],displayWVT=False,epsilon=-10)
                    wvt,ston=functions.generate_wvt3(binlist,si,va,np.full(len(binlist),1))
                    vwvt=functions.generate_wvt(binlist,va)
                edgex,edge2,edgec,edgea=radial_profile.radmethod(wvt,vwvt,wcsxlist[fil],show=False)
                #these are first pass, fit, ston median, and ston average. we want the last one.
                edges2.append(edgea)
            
            eem=np.nanmean(edges)
            if len(edges)>1:
                ees=np.nanstd(edges)
            else:
                ees=0

            eem2=np.nanmean(edges2)
            if len(edges2)>1:
                ees2=np.nanstd(edges2)
            else:
                ees2=0

            maxlist.append(str(round(eem,3))+"$\pm$"+str(round(ees,3)))
            simlist.append(str(round(eem2,3))+"$\pm$"+str(round(ees2,3)))
            numnanedges.append(str(int(nedge))+"/"+str(numtimes))
            if first<=eem2+ees2 and first<eem-ees and ((first-eem2)/ees2 < (first-eem)/ees or (nedge>3)):
                conclist.append("likely edge")
            else:
                conclist.append("inconclusive")
        tabulate(pathfile,sourcelist,objlist,firstlist,simlist,maxlist,numnanedges,conclist)
    
    print("Bye Bye!")