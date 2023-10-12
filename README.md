## pierrethx/MVT-binning
# MVT-binning: Modified Adaptive Binning Methods
Pierre Thibodeaux at University of California, Santa Barbara
under the guidance of Dr. Crystal Martin

We present a set of codes for adaptively binning 2D intensity maps to form less-resolved bins with greater signal-to-noise ratio (SNR). These methods are based off of Cappellari and Copin's Centroidal Voronoi Tessellation (2003) and Diehl and Statler's Weighted Voronoi Tessellation (2006) algorithms. These modifications allow the methods to be run on data that include negative values, which may result from procedures like continuum subtraction in emission line images.

Here, we treat the negative signal like zero signal during the binning, and then apply the formulated bins to the real data. This is based on the assumption that the positive and negative noise will be distributed approximately evenly in the bins and should cancel out when large enough bins are formed.

There are two stages of Adaptive Binning: Bin Accretion, which was created in Cappellari and Copin (2003), and Iteration, where the bins formed are regularized, the SNR is equalized between bins. Our method also masks bins which do not sufficiently form SNR, and bases the cutoff value on the shape of the SNR histogram over the images.

## Versions
July 2022

Version 1.0: Cleanup and consolodation of codes. First official version of code. 
5 October 2022

Version 1.1: Change in the cutoff method from defaultly finding minimum to finding maximum or secondary maximum, then tracking left to find minimum. Change in the 2-stage WVT method to mark the transition from high SNR to low SNR using "SNR density" as opposed to SNR. Addition of original.py script, minor modifications to testsuite.py, addition of license.

## Quickstart
Place all unbinned signal and variance files in a directory titled "unbinned" inside the desired file location. Run main.py on these files, and their outputs will be deposited in directories that are parallel to "unbinned", according to their target SNR. This program is best run in an astroconda (https://astroconda.readthedocs.io/en/latest/) virtual environment of Python 3.


## PRIMARY FILES

### MAIN.PY
**main.py** is the script to run to perform adaptive binning. When run, a tkinter GUI will appear, allowing the user to select multiple files to serve as the signal and variance input images. This will repeat until the user cancels the window, confirming all input sets of files to be used. It will also ask for target SNR levels, seperated by commas. The program will process EVERY set of files at EVERY target SNR levels, so it is better to break up long jobs that might get interrupted.

The other toggleable parameters of this program must be edited in the code, which consist of:
- the number of iterations OR convergence tolerance
- the minimum size of a bin (which depends on the minimum resolution element of the instrument which produced the images)
- which binning modes to use. There are four:
    1. 'CVT': Centroidal Voronoi Tessellation. Aside from our negative masking procedure, this is similar to Cappellari and Copin's formulation. This forms a Voronoi Tessellation then iteratively improves upon it by placing the generators of the tessellation at the centroids of each bin.
    2. 'WVT': Weighted Voronoi Tessellation. Aside from our negative masking procedure, this is similar to Diehl and Statler's formulation. This forms a Weighted Voronoi Tessellation (where the bin sizes are determined by an additional scalelength parameter) then iteratively improves upon it by placing recalculating the scalelengths (which replaces the process of moving the bin generators according to the centroids). It is advised one does not choose a minsize too big or else it will calculate poor initial scalelengths and bins.
    3. 'VT': Voronoi Tessellation. This can be interpreted as a CVT where we use geometric centers rather than centroids, or a WVT where we never calculate scalelength. This method was devised to preserve minimum bin size, since the previous two methods have difficulty achieving that, though does not attempt to tend towards better bin configurations, as the previous do. 
    4. 'WVT2s': Two Stage Voronoi Tessellation. This method was devised to preserve minimum bin size while also iterating (some) bins. It performs bin accretion, preserves bins above a certain SNR, then performs iteration on the unpreserved bins.
- The option to choose a new cutoff for the end bin mask. This is done by changing the check parameter of *saveblockoutfits*:
    - check=0 uses the code-calculated cut and outputs no image of the histogram
    - check=1 uses the code-calculated cut and outputs an image of the histogram
    - check=2 allows the user to override the cut and outputs an image of the histogram

This program outputs a variety of files:
1. the binned, masked signal fits file (file name is prefaced: "block_" and ends with "_sig.fits")
2. the binned, masked variance fits file (file name is prefaced: "block_" and ends with "_var.fits")
    - both of the above may be replaced with the unmasked versions by uncommenting *saveunblockedfits*
3. the binned SNR map fits file (file name is prefaced "zston_" and ends with ".fits")
    - may be removed by commenting out *saveston*
4. the binned assigned map fits file (file name is prefaced "z_" and ends with ".fits")
    - may be removed by commenting out *saveassign*. These files assign every bin a unique integer, with masked bins 0, and are used for the reconstruction of the bins outside of this program
5. The png image displaying the SNR histogram used to make the bin mask (file name is prefaced "block_" and ends with "_hist.png")
6. The png image of the binned SNR map (file name is prefaced "block_" and ends with "_ston.png")
    - both of the above toggled by changing the check parameter in *saveblockoutfits*. 0 outputs no image, 1 or 2 outputs both
7. The png image of the convergence over iterations (file name is prefaced "y_" and ends with "_convergence.png")
    - may be removed by commenting out *convergencelist*. This is to visualize how well the binning method converges, and is placed in "unbinned" if multiple targets are used.

### OTHER PRIMARY FILES
**bin_accretion.py** handles bin accretion and most of the tkinter/file entry.

**wvt_iteration.py** handles iteration, which differs depending on the selected mode(s)

**functions.py** holds functions that are used in the other files.

## AUXILIARY FILES
### DATA GENERATION
**circular_beta_modeler.py** is a little script that displays an editable Circular Beta profile so that one can gain an intuition for how the parameters modify the shape of the profile.

**generate_testdata.py** contains functions for generating simulated intensity data. It can generate data according to a Circular Beta profile, and also contains some simple cases for which the adaptive binnings performance can be tested.

### TEST PROCEDURES
**edge_detect.py** is a simple test that looks for the cessation of signal in a circularly symmetric image (such as in one of our simulated profiles) and outputs the radius at which that occurs. This would not apply to images with arbitrary shapes.

**test_vorbin.py** is a script that is to be used with Cappellari's VorBin package (https://pypi.org/project/vorbin/). It converts .fits file inputs into txt files compatible with the input files of the VorBin code and converts its output back into .fits files.

**original.py** implements main.py but disabling the mask and applying a minimum input SNR, in order to replicate Cappellari and Copin's CVT and Diehl and Statler's WVT methods but using the machinery of our code. We use this to directly compare our method to the "originals".

**testsuite.py** generates a pair of graphs which show SNR, bin roundness (as defined in Diehl and Statler 2006), and bin size as a function of radius. This is used for the visualization of how these parameters change over the image.

## License

MIT License

Copyright (c) 2022 Pierre Thibodeaux

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
