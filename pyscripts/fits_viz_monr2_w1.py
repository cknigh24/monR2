# Fits Visualization  MonR2 W1
#+
#-----------------------------------------------------------------------------
# These templates will look at an input fits file and output some quick statistics 
# and a histogram distribution to find an optimal scale to output an image  
# highlighting relevant emission within the FOV  
#-----------------------------------------------------------------------------
#-

from fits_viz_template import fits_viz as viz

viz(filename = '../fits/monr2_w1.fits',a1=6,a2=2400,b1 =1600,b2 = 4095)