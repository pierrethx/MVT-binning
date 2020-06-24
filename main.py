import numpy as np
import matplotlib.pyplot as plt 
import tkinter
from tkinter.filedialog import askopenfilename
from astropy.io import fits
import functions,bin_accretion
import scipy.spatial as sp

#from astropy.utils.data import download_file

sourcedir,signal,var=bin_accretion.initialize()
target=5 #Target Signal To noise
init_generators, init_scalelengths=bin_accretion.cc_accretion(sourcedir,signal,var,target)

fig=sp.voronoi_plot_2d(sp.Voronoi(init_generators),show_points=True,show_vertices=False)
plt.show()