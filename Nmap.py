import bin_accretion
import numpy as np
import matplotlib.pyplot as plt

wcsxlist,siglist,varlist,sourcelist,objlist=bin_accretion.minitialize()
stonlist=[siglist[i]/varlist[i] for i in range(len(siglist))]
snmax=[np.nanmax(stonlist[i]) for i in range(len(siglist))]
simax=[np.nanmax(siglist[i]) for i in range(len(siglist))]
vmax=[np.nanmax(varlist[i]) for i in range(len(siglist))]

snmed=[np.nanmedian(stonlist[i]) for i in range(len(siglist))]
simed=[np.nanmedian(siglist[i]) for i in range(len(siglist))]
vmed=[np.nanmedian(varlist[i]) for i in range(len(siglist))]

snmen=[np.nanmean(stonlist[i]) for i in range(len(siglist))]
simen=[np.nanmean(siglist[i]) for i in range(len(siglist))]
vmen=[np.nanmean(varlist[i]) for i in range(len(siglist))]
f,a=plt.subplots()
a.plot(range(0,len(snmax)),np.sort(snmax),label="SNR max")
a.plot(range(0,len(snmed)),np.sort(snmed),label="SNR median")
a.plot(range(0,len(snmen)),np.sort(snmen),label="SNR mean")
a.set_yscale("log")
a.set_title("SNR")
f.legend()
#print(snmax)

f2,a2=plt.subplots()
f3,a3=plt.subplots()
a2.plot(range(0,len(simax)),np.sort(simax),label="sig max")
a3.plot(range(0,len(simax)),np.sort(vmax),label="var max")
a2.plot(range(0,len(simed)),np.sort(simed),label="sig median")
a3.plot(range(0,len(simed)),np.sort(vmed),label="var median")
a2.plot(range(0,len(simen)),np.sort(simen),label="sig mean")
a3.plot(range(0,len(simen)),np.sort(vmen),label="var mean")
a2.set_yscale("log")
a2.set_title("sig")
a3.set_yscale("log")
a3.set_title("var")
#print(simax)
#print(vmax)
f2.legend()
f3.legend()
plt.show()
