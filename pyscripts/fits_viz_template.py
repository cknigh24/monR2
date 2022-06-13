# Fits Visualization templates
#+
#-----------------------------------------------------------------------------
# These templates will look at an input fits file and output some quick statistics 
# and a histogram distribution to find an optimal scale to output an image  
# highlighting relevant emission within the FOV  
#-----------------------------------------------------------------------------
#-


# Import libraries

import matplotlib.pyplot as plt
import numpy as np
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
np.set_printoptions(precision=4, suppress= True)
from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS
import matplotlib.ticker

#-----------------------------------------------------------------------------

# formatter for tickmarks
   
#see https://stackoverflow.com/questions/43324152/python-matplotlib-colorbar-scientific-notation-base
class OOMFormatter(matplotlib.ticker.ScalarFormatter):
    def __init__(self, order=0, fformat="%1.1f", offset=True, mathText=True):
        self.oom = order
        self.fformat = fformat
        matplotlib.ticker.ScalarFormatter.__init__(self,useOffset=offset,useMathText=mathText)
    def _set_order_of_magnitude(self):
        self.orderOfMagnitude = self.oom
    def _set_format(self, vmin=None, vmax=None):
        self.format = self.fformat
        if self._useMathText:
             self.format = r'$\mathdefault{%s}$' % self.format

#-----------------------------------------------------------------------------

# default code for histogram

def make_hist(image, pcent, units, plotpath, nametag, save=False):
    """
    NAME:
    make_hist
    
    DESCRIPTION:
    This function takes a 2-D array and outputs a histogram distribution
    for a given input range.
    
    INPUT:
    image = 2-D array to plot.
    pcent = Lower and upper percentiles of array. 
    units = Data units of image array. 
    plotpath = Path to output save plot. 
    nametag = Unique string to add to filename. 
    
    INPUT KEYWORD PARAMETERS:
    save = Boolean value to save plot.
    
    
    """
       
    image_vec=image.reshape(image.size) #reshape subimage into a 1D vector
    plt.rcParams['axes.facecolor'] = 'whitesmoke'
    fig1 = plt.figure(figsize=(12, 12))
    plt.hist(image_vec,bins = 50,range=(pcent[0],pcent[1] ),color='tab:blue')
    plt.xlabel(str(units))
    plt.ylabel('Counts')
    plt.grid(color='tab:grey')
    plt.show()
    if save == True:
        fig1.savefig(plotpath + nametag + "_hist.png", dpi=300, bbox_inches='tight') # add more file output options
        print("Saving histogram titled: '" + plotpath + nametag + "_hist.png'" ) 

#-----------------------------------------------------------------------------

# default code for cumulative histogram
def make_cumu_hist(image, pcent, units, plotpath, nametag, save=False):
    """
    NAME:
    make_cumu_hist
    
    DESCRIPTION:
    This function takes a 2-D array and outputs a cumulative histogram 
    distribution for a given input range.
    
    INPUT:
    image = 2-D array to plot.
    pcent = Lower and upper percentiles of array.
    units = Data units of image array.
    plotpath = Path to output save plot.
    nametag = Unique string to add to filename.
    
    INPUT KEYWORD PARAMETERS:
    save = Boolean value to save plot. 
    
    """
    image_vec=image.reshape(image.size) #reshape subimage into a 1D vector
    plt.rcParams['axes.facecolor'] = 'whitesmoke'
    fig2 = plt.figure(figsize=(12, 12))
    plt.hist(image_vec,bins = 50,range=(pcent[0],pcent[1] ),color='tab:blue', cumulative = True, histtype= 'stepfilled' )
    plt.xlabel (str(units))
    plt.ylabel('Counts')
    plt.grid(color='tab:grey')
    plt.show()
    if save == True:
        fig2.savefig(plotpath + nametag + "_cumu_hist.png", dpi=300, bbox_inches='tight') # add more file output options
        print("Saving Cumulative histogram titled: '" + plotpath + nametag + "_cumu_hist.png'" ) 

#-----------------------------------------------------------------------------

# default code to plot a 2d image




def make_image_plot(image, figx, figy, subwcs, 
                    tickspacing, scaling, cmap, pcent, units,
                    plotpath, nametag, fontsize, hdu,
                    contourname1 = None, contour_levels1 = None,
                    contourname2= None, contour_levels2 =None, 
                    coordx = [], coordy=[], grid =False, bad_pix_mark = False, 
                    wcs_axes = False, rot_wcs = False, save=False):
    """
    NAME:
    make_image_plot
    
    DESCRIPTION:
    This function takes a 2-D array and outputs an image plot.
    
    INPUT:
    image = 2-D array to plot.
    figx = x dimension of figure window.
    figy = y dimension of figure window.
    subwcs = wcs corresponding to 2-D array.
    tickspacing = 2 component array to specify spacing between major tick marks 
        in x and y axes (or wcs axes). 
    scaling = scaling of color bar [supports linear, sqrt, or log normalizations] 
    cmap =  color map used to visualize image
    pcent = Lower and upper percentiles of array.
    units = Data units of image array.
    plotpath = Path to output save plot.
    nametag = Unique string to add to filename.
    fontsize = fontsize of axes labels.
    hdu = header data unit of fits files to read.
    contourname1 = filename of fits file to use as first set of overplotted contours.
    contourname2 = filename of fits file to use as second set of overplotted contours.
    contour_levels1 = array of numeric values to use in first set of overplotted contours.
    contour_levels2 = array of numeric values to use in second set of overplotted contours
    coordx = array of x coordinates of pixels to mask out.
    coordy = array of y coordinates of pixels to mask out (must match coordx).
    
    INPUT KEYWORD PARAMETERS:
    grid =   Boolean value to include overplotted coordinate grid. 
    bad_pix_mark =  Boolean value to include bad pixel masking. 
    wcs_axes = Boolean value to use astronomical coordinates in axes.
    rot_wcs = Boolean value to rotate wcs axes in cases where image is rotated by 90 degrees
    save = Boolean value to save plot.
    
    """
    import matplotlib.colors as colors
    fig3 =plt.figure(figsize=(figx,figy),tight_layout=True)
    if wcs_axes == False:
        #ax =plt.gca()
        ax = fig3.add_subplot(111)
        plt.xlabel(r'X (Pixels)',fontsize=fontsize)
        plt.ylabel(r'Y (Pixels)',fontsize=fontsize)
        ax.tick_params(labelsize =fontsize)
    else:
        ax = fig3.add_subplot(111,projection=subwcs)
        ra =ax.coords['ra']
        dec= ax.coords['dec']
        dec.set_axislabel(r'Declination',fontsize=fontsize)
        #for cases where image grid is rotated wrt to wcs
        if rot_wcs == True: 
            ra.set_axislabel_position('l')
            dec.set_axislabel_position('b')
            ra.set_ticklabel_position('l')
            dec.set_ticklabel_position('b')
            ra.set_ticks_position('l')
            dec.set_ticks_position('b')
        else:
            ra.set_ticklabel_position('bl')
            dec.set_ticklabel_position('bl')
            ra.set_ticks_position('bl')
            dec.set_ticks_position('bl')
        ra.set_axislabel(r'Right Ascension',fontsize=fontsize)
        ra.set_ticks(spacing=tickspacing[0] * u.hourangle, color='black', size = 8,
                     width =1, direction ='out')
        dec.set_ticks(spacing=tickspacing[1] * u.arcsec, color='black', size =8,
                      width =1, direction ='out')
        ra.set_ticklabel(exclude_overlapping=True, size = 8)
        dec.set_ticklabel(exclude_overlapping=True, size = 8)
        
    if scaling == "sqrt":
        norm = colors.PowerNorm(gamma=0.5,  vmax= pcent[1], vmin =pcent[0])
        im=ax.imshow(image, cmap=cmap,origin='lower',norm=norm) 
    elif scaling == "log":
        norm = colors.LogNorm(vmax= pcent[1], vmin =pcent[0])
        im=ax.imshow(image, cmap=cmap,origin='lower',norm=norm)
    elif scaling == "linear":
        norm = None
        im=ax.imshow(image, cmap=cmap,origin='lower',
                     vmax= pcent[1], vmin =pcent[0])
              
             
    im_ratio = image.shape[0]/image.shape[1]        
        
    if scaling != 'log':
        v1 = np.floor(np.log10(np.nanmax(image)))
        cbar=fig3.colorbar(im, pad=0.02,fraction=0.046*im_ratio,
                        format=OOMFormatter(v1, mathText=False))
    if scaling == 'log':
     cbar=fig3.colorbar(im, pad=0.02,fraction=0.046*im_ratio)
         
    cbar.ax.tick_params(length =8, width =2,direction ='out',which = 'major')   
    
    # read in fits file for contours if desired (must be same size as image file)
    if contourname1 != None:
       header_data_unit_list_c1 = fits.open(contourname1)
       contours1 = header_data_unit_list_c1[hdu].data
       header_c1 = header_data_unit_list_c1[hdu].header
       if contour_levels1 == None:
           contour_levels1 = [0.25*i*np.nanmax(contours1) for i in range(1,5)]
           
    if contourname2 != None:
       header_data_unit_list_c2 = fits.open(contourname2)
       contours2 = header_data_unit_list_c2[hdu].data
       header_c2 = header_data_unit_list_c2[hdu].header
       if contour_levels2 == None:
           contour_levels2 = [0.25*i*np.nanmax(contours2) for i in range(1,5)] 
     
    if wcs_axes != True:
        ax.set_autoscale_on(False)
        if contourname1 != None:
            ax.contour(contours1, levels = contour_levels1, colors=('green'),
                   linestyles=('solid'))
        if contourname2 != None:
            ax.contour(contours2, levels = contour_levels2, colors=('cyan'),
                   linestyles=('solid'))
       
            
             
    if wcs_axes == True:    
        ax.set_autoscale_on(False)
        if contourname1 != None:
            ax.contour(contours1, levels = contour_levels1, colors=('green'),
                   linestyles=('solid'),transform=ax.get_transform(WCS(header_c1)))
           
        if contourname2 != None:
            ax.contour(contours2, levels = contour_levels2, colors=('cyan'),
                   linestyles=('solid'),transform=ax.get_transform(WCS(header_c2)))
           
    # Overplot markers for bad pixels
    if bad_pix_mark == True:
        from matplotlib.patches import Rectangle
        for i in range(0, len(coordx)):
            r = Rectangle((coordx[i], coordy[i]), 1., 1., edgecolor='red', facecolor='lightgrey')
            ax.add_patch(r)
        
        
         
    if units != None:
       cbar.set_label(units)
         
    if grid == True:
        ax.grid(color='black', ls='dotted', lw = 2)
    plt.show()
         
    if save == True:
        fig3.tight_layout()
        fig3.savefig(plotpath + nametag + "_im.png", dpi=300, bbox_inches='tight')
        print("Saving image titled: '" + plotpath + nametag + "_im.png'"  )

#-----------------------------------------------------------------------------

# default code to save modified fits file

def output_to_fits(image, header, subwcs, units, filename, fitspath, nametag):
    """
    NAME:
    output_to_fits
    
    DESCRIPTION:
    This function takes a 2-D array and outputs it to a fits file.
    
    INPUT:
    image = 2-D array to  output as data in fits file.
    header = Header structure to use as template header in output fits file.
    subwcs = wcs corresponding to 2-D array.
    units = Data units of image array.
    filename = fits file to read in as template header data unit list.
    fitspath = Path to output save fits.
    nametag = Unique string to add to filename.
    
    
    """
    header_data_unit_list = fits.open(filename)
    header_new=header.copy() # Need to update header
    subshape=subwcs.array_shape
    c_x= subshape[1]//2
    c_y= subshape[0]//2
    c_ra,c_dec=subwcs.pixel_to_world_values(c_x, c_y)
    c_ra,c_dec=float(c_ra),float(c_dec)
    header_new["BUNIT"] = str(units)
    header_new["CRVAL1"] = c_ra
    header_new["CRVAL2"] = c_dec
    header_new["CRPIX1"] = c_x
    header_new["CRPIX2"] = c_y
    header_new["NAXIS1"] = subshape[0]
    header_new["NAXIS2"] = subshape[1]
    header_data_unit_list[0].header = header_new
    header_data_unit_list[0].data = image
    if fitspath != "":
        header_data_unit_list.writeto(fitspath + nametag + "_cut.fits", overwrite=True)
        print("Saving fits titled: '" + fitspath + nametag + "_cut.fits'"  )
    if fitspath == "":
        header_data_unit_list.writeto(nametag + "_cut.fits", overwrite=True)
        print("Saving fits titled: '" + nametag + "_cut.fits'"  )
    print("\n=================================================================\n")    

#-----------------------------------------------------------------------------

# Define function that can output simple visualization of fits files
def fits_viz(filename,errorname= None, a1 = 0, b1 = 0, a2= 0, b2= 0,  
            scaling ="linear", lower_p =0.25, upper_p = 99.75, cmap='plasma', 
             plotpath = '', nametag = "", fitspath = "",  hdu = 0, 
            convert_factor = None, pixel_size = None, units = None, sigma_cut = [3,False],
            fontsize = 18, figx = 10, figy = 10,tickspacing = [2.7777777e-4, 15],
            contourname1 = None , contour_levels1=None,
            contourname2 = None , contour_levels2=None,
            coordx = [], coordy=[], grid = False, bad_pix_mark = False,  
            wcs_axes = False,rot_wcs = False, hist_display = False, 
            cumhist_display = False, im_display = True, show_stats = False, 
            save = False, save_fits = False, print_hdr = False, 
            extract_units = False):
    """
    NAME:
    fits_viz
    
    DESCRIPTION:
    This function is an overlay for multiple visualization 
    i/o procedures involving fits files, and wcs manipulation.
    
    INPUT:
    filename = string filename pertaining to 2-d fits array to visualize.
    errorname = string filename pertaining error plane of 2-d fits array
        to visualize. Only relevant when using sigma cut procedure. 
    a1 = starting row index of 2-D array to plot.
    a2 = ending row index of 2-D array to plot.
    b1 = starting column index of 2-D array to plot.
    b2 = ending column index of 2-D array to plot.
    scaling = scaling of color bar [supports linear, sqrt, or log normalizations].
    lower_p =  Lower percentile of array to consider.
    upper_p =  Upper percentile of array to consider.
    cmap =  color map used to visualize image
    plotpath = Path to output save plot.
    nametag = Unique string to add to filename.
    fitspath = Path to output save fits.
    hdu = header data unit of fits files to read.
    convert_factor = convert from one set of units to another (currently limited).
    pixel_size = pixel size of input array (used in unit_conversions)
    units = Data units of image array.
    sigma_cut = 2 component list to implement masking based on signal-to-noise limits.
        First component specifies the desired signal-to-noise value to mask image.
        Second component is a boolean value to implement this mask if true.
        Note: file for error plane must be included to use this procedure.
    fontsize = fontsize of axes labels.
    figx = x dimension of figure window.
    figy = y dimension of figure window.
    tickspacing = 2 component list to specify spacing between major tick marks 
        in x and y axes (or wcs axes) respectively. 
    contourname1 = filename of fits file to use as first set of overplotted contours.
    contourname2 = filename of fits file to use as second set of overplotted contours.
    contour_levels1 = array of numeric values to use in first set of overplotted contours.
    contour_levels2 = array of numeric values to use in second set of overplotted contours.
    coordx = array of x coordinates of pixels to mask out.
    coordy = array of y coordinates of pixels to mask out (must match coordx).
    
    INPUT KEYWORD PARAMETERS:
    grid =   Boolean value to include overplotted coordinate grid.     
    bad_pix_mark =  Boolean value to include bad pixel masking. 
    wcs_axes = Boolean value to use astronomical coordinates in axes.
    rot_wcs = Boolean value to rotate wcs axes in cases where image is rotated by 90 degrees
    save = Boolean value to save plot.
    hist_display = Boolean value to make histogram plot.
    cumu_hist_display = Boolean value to make cumulative histogram plot.
    im_display = Boolean value to make image plot.
    show_stats =  Boolean value to output simple statistics on imported fits image.
    save = Boolean value to save plot.
    save_fits = Boolean value to save image to new fits file.
    print_hdr = Boolean value to output header information.
    extract_units = Boolean value to extract units for fits header
    
    """
    

   
    
    # read in designated hdu of fits file to visualize
    header_data_unit_list = fits.open(filename)
    header_data_unit_list.info()
    image = header_data_unit_list[hdu].data
    header = header_data_unit_list[hdu].header
    
    
    # print out formatted header info if requested
    if print_hdr == True:
        print("\n=================================================================\n")
        print(header.index)
    print("\n=================================================================\n")
    
    
    # Extract nametag from filename for output file 
    if nametag == "":
        nametag = [word for line in filename.split(".") for word in line.split("/")][-2]
        print("nametag: " + nametag)
        plotpath +nametag
        
    print("\n=================================================================\n")
    
    # apply image array cut (default = None)
    if (a2 == 0 or b2 == 0):
        a2 = header["NAXIS1"]
        b2 = header["NAXIS2"]
        image = image[b1:b2,a1:a2]
    else:
        image = image[b1:b2,a1:a2]
      
    
    # read in error fits if sigma cut is desired
    if sigma_cut[1] == True and errorname != None:
        header_data_unit_list_e = fits.open(errorname)
        header_data_unit_list.info()
        errors = header_data_unit_list_e[hdu].data
        errors = errors[b1:b2,a1:a2]
        #header_e = header_data_unit_list_e[hdu].header
        print, "reading in associated error file: "+errorname+" to perform a " +str(sigma_cut[0])+ " sigma cut."
        snr = image/errors
        for i in range(b1,b2):  
            for j in range(a1,a2):
                if snr[i,j] <= sigma_cut[0]:
                    image[i,j] = np.nan
                    
     
    # Get units from header
    if extract_units == True:
        units = header["BUNIT"]
    
    # Optional call to unit converter here
    # default unit conversion: DN to MJy/sr
    if convert_factor != None and pixel_size != None:
        if units != None:
            from unit_converter import convert
            image =convert(image,cf = convert_factor, pixel_size = pixel_size,
                           in_unit = "DN/pixel", out_unit = "MJy/sr", verbose=False)
            units =  "MJy/sr"
        
    # Select colorbar range based on percentiles    
    pcent = np.nanpercentile(image, [lower_p, upper_p]) 
    
    
    #Do some quick statistics excluding NaNs
    if show_stats == True:
        print("Image units =  {}".format(units))
        print('max = ',np.nanmax(image), " ",str(units))
        print('min = ',np.nanmin(image), " ",str(units))
        print('range =',np.nanmax(image) -np.nanmin(image), " ",str(units))
        print('mean = ',np.nanmean(image), " ",str(units) )
        print('median = ',np.nanmedian(image), " ",str(units) )
        print('mean = ',np.nanmean(image), " ",str(units) )
        print('std = ', np.nanstd(image), " ",str(units))
        print("{}th percentile = {} {} & {}th percentile = {} {}"
              .format(lower_p, pcent[0],units,upper_p, pcent[1],units))
        print("\n=================================================================\n")      
        
    #Make histogram of emission distribution
    plt.rcParams["axes.grid"] = False
    if hist_display == True:
        make_hist(image, pcent, units,plotpath,nametag, save)  
        
    #Make a cumulative histogram of emission distribution
    if cumhist_display == True:      
       make_cumu_hist(image, pcent, units,plotpath,nametag, save)  
       
    # Get WCS information and update if subimage is desired
    im_wcs = WCS(header)    
    #Slice wcs to match latest subimage
    from astropy.wcs.wcsapi import SlicedLowLevelWCS
    slices=[slice(a1,a2), slice(b1,b2)]
    subwcs = SlicedLowLevelWCS(im_wcs,slices=slices )

        
    #plot image of fits array 
        
        
    # Can define different colormapping options, ranges, and colortables from function call 
    if im_display == True:
       make_image_plot(image, figx, figy, subwcs, 
                           tickspacing, scaling, cmap, pcent, units,
                           plotpath, nametag, fontsize, hdu,
                           contourname1, contour_levels1, contourname2, contour_levels2, 
                           coordx, coordy, grid, bad_pix_mark, 
                           wcs_axes, rot_wcs, save)

    # Output subimage fits file
    if save_fits == True:
        output_to_fits(image, header, subwcs, units, filename, fitspath, nametag)