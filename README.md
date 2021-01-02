# Adaptive-Binning
Testing out adaptive binning procedures.
Pierre Thibodeaux at University of California, Santa Barbara
under the guidance of Dr. Crystal Martin

bin_accretion.py is a program that performs the Bin Accretion algorithm, as defined by Cappellari and Copin (2003).
wvt_iteration.py is a program that performs the iterative WVT algorithm, as defined by Diehl and Statler (2006).
functions.py is a reference for functions used in the above two programs.
main.py is a program that runs bin_accretion and wvt_iteration on a file to output a binned file and its effects. This specifically runs the data mask method, one of the two devised ways to handle negative data. Here, we treat the negative signal like zero signal during the binning, and then apply the formulated bins to the real data, with the assumption that the negative signal will be distributed approximately evenly in the bins.
queuemain.py runs the main.py sequence for a set of files and a set of target S/N. It uses the functions from main, but is better equipped to quickly process multiple files.
queuemain2.py runs a procedure modified from main.py. It instead uses the expansion method, which expands bins that are 0 or negative, so that they might find more signal elsewhere.

generatecolori.py generates a suite of sh scripts that can be used to look at various line images after binning them. ds9 doesnt allow you to save RGB images, so in order to access them with whatever parameters you want, you need to access them through a sh script using line commands. It is kind of tedious to set them all up and name them, so this does that. prescript.sh is meant to be easy to modify to direct to the users instance of ds9 so that they dont have to edit each colorimage script.
plotit.py has done many various things but it is currently configured to make a StoN file when given an unbinned signal and variance. You can technically apply this to binned images but the StoN for a bin is not simply the S of one of its pixels over sqrt(V) for that pixel, the entire bin needs to be considered.

generatetestdata.py generates a set of Circular Beta Profiles (Sarazin 1988) for the purpose of binning them to see how the code responds. Generates profiles with varying parameters. The variance is modified with an extra parameter to better match the slim variance of real data.
generatetestcases.py generates a set of test cases to check that the code works for simple geometric situations.
showiteration.py constructs a gif to visualize the wvt construction process
circular_beta_modeler.py displays a circular beta profile with mutable parameters to visualize the curve and develop a fitting procedure
3d_visualizer.py displays a .fits file as a 3-d plot.
fullbins.py generates a panel of images to compare various levels of binning.

radial_profile.py generates a radial profile for a binned or unbinned signal file and fits it to a Circular Beta Profile (Sarazin 1988)
qradial.py performs the function of radial_profile on multiple files and plots recovered edge data on a plot.
queuepipeline takes unbinned files and performs queuemain and qradial on them to generate a plot of recovered edge data.

radtest.py is the precursor to testcheck2.py
scalestest.py runs through the binning with 4 different formulations of the calculate_scales function (which is the step in the WVT that is most directly hindered by the presence of negative data)
testcheck1.py checks that signal and variance are conserved from an unbinned image to its binned product.
testcheck2.py looks at StoN as a function of radius from center for binned images and also looks at histograms to show overall brightness/StoN distributions for binned/unbinned images.
testcheck3.py looks at the roundness of the bins, following the procedure in DS06.
testcheck4.py is a modification of testcheck2.py to generate a more informative figure that compares the results of different binning levels alongside each other. Following the writeup, 4 should correspond to convergence of the program, but that must be constructed and saved at runtime unlike the rest of the tests, which are meant to be applied after the binning has been completed. So I just used this empty file for something else



