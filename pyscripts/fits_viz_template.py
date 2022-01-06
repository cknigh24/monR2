# Fits Visualization templates
#+
#-----------------------------------------------------------------------------
# These templates will look at an input fits file and output some quick statistics 
# and a histogram distribution to find an optimal scale to output an image  
# highlighting relevant emission within the FOV  
#-----------------------------------------------------------------------------
#-





# Define function that can output simple visualization of fits files
def fits_viz(filename,errorname= None, a1 = 0, b1 = 0, a2= 0, b2= 0, hist_display = False, 
            cumhist_display = False, im_display = True, show_stats = True,  
            scaling ="sqrt", lower_p =0.25, upper_p = 99.75, cmap='plasma', 
            grid = False, save = False,save_fits = False, plotpath = '', 
            nametag = "", fitspath = "", print_hdr = False, hdu = 0, 
            wcs_axes = False, convert_factor = None, pixel_size = None, 
            units = None, extract_units = False, sigma_cut = [3,False],
            fontsize = 18, figx = 10, figy = 10,tickspacing = [2.7777777e-4, 15],
            contourname1 = None , contour_levels1=None,
            contourname2 = None , contour_levels2=None):

    
    # Import libraries
    import matplotlib.pyplot as plt
    import numpy as np
    from astropy.visualization import astropy_mpl_style
    plt.style.use(astropy_mpl_style)
    np.set_printoptions(precision=4, suppress= True)
    from astropy.io import fits
    from astropy import units as u
    import matplotlib.ticker
  
  
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
        print([a1,a2,b1,b2])
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
                    
    # read in fits file for contours if desired (must be same size as image file)
    if contourname1 != None:
        header_data_unit_list_c1 = fits.open(contourname1)
        contours1 = header_data_unit_list_c1[hdu].data
        header_c1 = header_data_unit_list_c1[hdu].header
        #contours1 = contours1[b1:b2,a1:a2]
        if contour_levels1 == None:
            contour_levels1 = [0.25*i*np.nanmax(contours1) for i in range(1,5)]
            
    if contourname2 != None:
        header_data_unit_list_c2 = fits.open(contourname2)
        contours2 = header_data_unit_list_c2[hdu].data
        header_c2 = header_data_unit_list_c2[hdu].header
       # contours2 = contours2[b1:b2,a1:a2]
        if contour_levels2 == None:
            contour_levels2 = [0.25*i*np.nanmax(contours2) for i in range(1,5)]
                    
    
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
        image_vec=image.reshape(image.size) #reshape subimage into a 1D vector
        plt.rcParams['axes.facecolor'] = 'whitesmoke'
        fig1 = plt.figure(figsize=(12, 12))
        plt.hist(image_vec,bins = 50,range=(pcent[0],pcent[1] ),color='tab:blue')
        plt.xlabel(str(units))
        plt.ylabel('Counts')
        plt.grid(color='tab:grey')
        plt.show()
        print("Figure 1: Histogram \n")
        if save == True:
            fig1.savefig(plotpath + nametag + "_hist.png", dpi=300, bbox_inches='tight') # add more file output options
            print("Saving histogram titled: '" + plotpath + nametag + "_hist.png'" )
    
    
    #Make a cumulative histogram of emission distribution
    if cumhist_display == True:      
        plt.rcParams['axes.facecolor'] = 'whitesmoke'
        fig2 = plt.figure(figsize=(12, 12))
        plt.hist(image_vec,bins = 50,range=(pcent[0],pcent[1] ),color='tab:blue', cumulative = True, histtype= 'stepfilled' )
        plt.xlabel (str(units))
        plt.ylabel('Counts')
        plt.grid(color='tab:grey')
        plt.show()
        print("Figure 2: Cumulative Histogram \n")
        if save == True:
            fig2.savefig(plotpath + nametag + "_cumu_hist.png", dpi=300, bbox_inches='tight') # add more file output options
            print("Saving Cumulative histogram titled: '" + plotpath + nametag + "_cumu_hist.png'" )
    
     
    # Apply WCS to show sky coordinates on image axes?
    if wcs_axes == True or save == True: # 
        from astropy.wcs import WCS
        im_wcs = WCS(header)
        print(im_wcs)
        print("\n=================================================================\n")
        
        #Slice wcs to match latest subimage
        from astropy.wcs.wcsapi import SlicedLowLevelWCS
        slices=[slice(a1,a2), slice(b1,b2)]
        subwcs = SlicedLowLevelWCS(im_wcs,slices=slices )
        
        
    #plot image of fits array 
    
    
    # Can define different colormapping options, ranges, and colortables from function call 
    if im_display == True:
        import matplotlib.colors as colors
        fig3 =plt.figure(figsize=(figx,figy),tight_layout=True)
        if wcs_axes == False:
            ax =plt.gca()
            plt.xlabel(r'X (Pixels)',fontsize=fontsize)
            plt.ylabel(r'Y (Pixels)',fontsize=fontsize)
            ax.tick_params(labelsize =fontsize)
        else:
            ax = fig3.add_subplot(111,projection=subwcs)
            ra =ax.coords[1]
            dec= ax.coords[0]
            ra.set_axislabel(r'Declination',fontsize=fontsize)
            dec.set_axislabel(r'Right Ascension',fontsize=fontsize)
            ra.set_ticks(spacing=tickspacing[0] * u.hourangle, color='black', size = 8,
                         width =1, direction ='out')
            dec.set_ticks(spacing=tickspacing[1] * u.arcsec, color='black', size =8,
                          width =1, direction ='out')
            ra.set_ticklabel(exclude_overlapping=True, size = 8)
            dec.set_ticklabel(exclude_overlapping=True, size = 8)
            ra.set_ticks_position('bl')
            dec.set_ticks_position('bl')

            #ax.tick_params(labelsize =fontsize)
        
           # labels = ax.get_xticklabels()
            #ax.setp(labels, rotation=45, ha="right" rotation_mode="anchor")
       
           
            
        if scaling == "sqrt":
            norm = colors.PowerNorm(gamma=0.5,  vmax= pcent[1], vmin =pcent[0])
            im=ax.imshow(image, cmap=cmap,origin='lower',norm=norm) 
        elif scaling == "log":
            norm = colors.LogNorm(vmax= pcent[1], vmin =pcent[0])
            im=ax.imshow(image, cmap=cmap,origin='lower',norm=norm)
        elif scaling == "linear":
            norm = None
            im=ax.imshow(image, cmap=cmap,origin='lower',
                         vmax= pcent[1], vmin =pcent[0], extent = [a1,a2,b1,b2]) 
         
        
        im_ratio = image.shape[0]/image.shape[1] 
     
        v1 = np.floor(np.log10(np.nanmax(image)))
        


        cbar=fig3.colorbar(im, pad=0.02,fraction=0.046*im_ratio,
                            format=OOMFormatter(v1, mathText=False))
        
       
#        if scaling != "log":
#            cbar.formatter.set_powerlimits((np.nanmin(image),np.nanmax(image)))
#        
#        
#        
#        cbar.ax.set_yticklabels(['{:.1e}'.format(x) for x in v1],
#                                 fontsize=14, weight='bold')

#        for l in cbar.ax.yaxis.get_ticklabels():
#            l.set_weight("bold")
#            l.set_fontsize(16)
            
        ax.set_autoscale_on(False)
        if contourname1 != None:
            ax.contour(contours1, levels = contour_levels1, colors=('k',),
                       linestyles=('solid'),transform=ax.get_transform(WCS(header_c1)))
          
        if contourname2 != None:
            ax.contour(contours2, levels = contour_levels2, colors=('lightgrey'),
                       linestyles=('solid'),transform=ax.get_transform(WCS(header_c2)))
           
        #cbar.formatter.set_powerlimits(round(np.nanmin(image)),round(np.nanmax(image)))
       # cbar.ax.set_yticklabels(["{:4.2e}".format(i) for i in v1])
       
        
        if units != None:
           cbar.set_label(units)
        
        if grid == True:
            ax.grid(color='dimgrey', ls='dotted', lw = 2)
        plt.show()
        
        print("Figure 3: Image Display \n")
        if save == True:
            fig3.tight_layout()
            fig3.savefig(plotpath + nametag + "_im.png", dpi=300, bbox_inches='tight')
            print("Saving image titled: '" + plotpath + nametag + "_im.png'"  )

        
   
            
    '''
    Some colormap suggestions:
        inferno, magma, plasma, hot
    '''

    # Output subimage fits file
    if save_fits == True:
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
