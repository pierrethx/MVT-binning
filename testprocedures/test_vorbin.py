#!/usr/bin/env python

"""
Copyright (C) 2003-2014, Michele Cappellari
E-mail: michele.cappellari_at_physics.ox.ac.uk

    V1.0.0: Michele Cappellari, Vicenza, 13 February 2003
    V1.0.1: Use astro library routines to read and write files.
        MC, Leiden, 24 July 2003
    V2.0.0: Translated from IDL into Python. MC, London, 19 March 2014
    V2.0.1: Support both Python 2.6/2.7 and Python 3.x. MC, Oxford, 25 May 2014
    V2.0.2: Make files paths relative to this file, to run the example from
        any directory. MC, Oxford, 23 January 2017
    V2.0.3: Changed imports for vorbin as package. 
        Make file paths relative to the vorbin package to be able to run the
        example unchanged from any directory. MC, Oxford, 17 April 2018    
    V2.0.4: Dropped legacy Python 2.7 support. MC, Oxford, 10 May 2018

"""

import numpy as np
from astropy.io import fits
import sys,os
sys.path.append(("/".join(os.path.realpath(__file__).split("/")[:-2])))
import bin_accretion
from vorbin.voronoi_2d_binning import voronoi_2d_binning

#-----------------------------------------------------------------------------

def make_input(signal,variance,file_dir,objname):

    with open(file_dir+'/'+objname+'_input.txt', 'w') as f:
        f.write("##########################################################\n")
        f.write("#          X             Y          Signal       Noise\n")
        f.write("##########################################################\n")
        for y in range(len(signal)):
            for x in range(len(signal[0])):
                if np.sqrt(variance[y][x])>0:
                    f.write("       ")
                    f.write(str(x+1))
                    f.write("       ")
                    f.write(str(y+1))
                    f.write("       ")
                    f.write(str(signal[y][x]))
                    f.write("       ")
                    f.write(str(np.sqrt(variance[y][x])))
                    f.write('\n')

def make_output(signal,variance,file_dir,objname):
    x, y, binnum = np.loadtxt(file_dir + '/'+objname+'_output.txt').T
    
    out=np.zeros((int(y[-1]),int(x[-1])))

    ss=np.zeros(int(np.nanmax(binnum)+1))
    vs=np.zeros(int(np.nanmax(binnum)+1))
    for i in range(len(binnum)):
        ss[int(binnum[i])]+=signal[int(y[i]-1)][int(x[i]-1)]
        vs[int(binnum[i])]+=variance[int(y[i]-1)][int(x[i]-1)]
    
    for i in range(len(binnum)):
        out[int(y[i]-1)][int(x[i]-1)]=ss[int(binnum[i])]/np.sqrt(vs[int(binnum[i])])


    hdu3 = fits.PrimaryHDU(out)
    hdul3 = fits.HDUList([hdu3])
    #manipulate(hdul3)
    hdul3.writeto(file_dir+'/vorbin_output_'+objname+'.fits',overwrite=True,checksum=True)

def voronoi_binning_example():
    """
    Usage example for the procedure VORONOI_2D_BINNING.

    It is assumed below that the file voronoi_2d_binning_example.txt
    resides in the current directory. Here columns 1-4 of the text file
    contain respectively the x, y coordinates of each SAURON lens
    and the corresponding Signal and Noise.

    """
    ##file_dir = "/Users/pierre/Downloads/vorbin-3.1.4/example/unbinned"  # path of vorbin

    wcsx,signal,var,sourcedir,objname=bin_accretion.initialize()
    targetSN=bin_accretion.gettarget()
    #print(sourcedir)
    #print(os.path.exists(sourcedir))
    subfolder=bin_accretion.makesubfolder(sourcedir,targetSN)

    make_input(signal,var,sourcedir+"/"+subfolder,objname)
    x, y, signal, noise = np.loadtxt(sourcedir+"/"+subfolder + '/'+objname+'_input.txt').T

    # Perform the actual computation. The vectors
    # (binNum, xNode, yNode, xBar, yBar, sn, nPixels, scale)
    # are all generated in *output*
    #
    binNum, xNode, yNode, xBar, yBar, sn, nPixels, scale = voronoi_2d_binning(
        x, y, signal, noise, targetSN, plot=1, quiet=0,pixelsize=1)

    # Save to a text file the initial coordinates of each pixel together
    # with the corresponding bin number computed by this procedure.
    # binNum uniquely specifies the bins and for this reason it is the only
    # number required for any subsequent calculation on the bins.
    #
    np.savetxt(sourcedir+"/"+subfolder + '/'+objname+'_output.txt', np.column_stack([x, y, binNum]),
               fmt=b'%10.6f %10.6f %8i')
    make_output(signal,var,sourcedir+"/"+subfolder,objname)

#-----------------------------------------------------------------------------

if __name__ == '__main__':
    voronoi_binning_example()