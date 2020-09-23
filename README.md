# Adaptive-Binning
Testing out adaptive binning procedures.
Pierre Thibodeaux at University of California, Santa Barbara

bin_accretion.py is a program that performs the Bin Accretion algorithm, as defined by Cappellari and Copin (2003).
wvt_iteration.py is a program that performs the iterative WVT algorithm, as defined by Diehl and Statler (2006).
functions.py is a reference for functions used in the above two programs.
main.py is a program that runs bin_accretion and wvt_iteration on a file to output a binned file and its effects.
queuemain.py runs the main.py sequence for a set of files and a set of target S/N.

radial_profile.py generates a radial profile for a binned or unbinned signal file and fits it to a Circular Beta Profile (Sarazin 1988)
qradial.py performs the function of radial_profile on multiple files and plots recovered edge data on a plot.
queuepipeline takes unbinned files and performs queuemain and qradial on them to generate a plot of recovered edge data.

showiteration.py constructs a gif to visualize the wvt construction process
circular_beta_modeler.py displays a circular beta profile with mutable parameters to visualize the curve and develop a fitting procedure
3d_visualizer.py displays a .fits file as a 3-d plot.
generatetestdata.py generates test data using the circulat beta profile to run through the binning procedure.

