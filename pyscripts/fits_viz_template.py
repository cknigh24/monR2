# Fits Visualization templates
#+
#-----------------------------------------------------------------------------
# These templates will look at an input fits file and output some quick statistics 
# and a histogram distribution to find an optimal scale to output an image  
# highlighting relevant emission within the FOV  
#-----------------------------------------------------------------------------
#-





# Define function that can output simple visualization of fits files
def fits_viz(filename, a1 = 0, b1 = 0, a2= 0, b2= 0, hist_display = False, 
            cumhist_display = False, im_display = True, show_stats = True,  
            scaling ="sqrt", lower_p =0.25, upper_p = 99.75, cmap='cividis', 
            grid = False, save = False, plotpath = '', nametag = "",
            fitspath = "", print_hdr = False, hdu = 0, wcs_axes = False):

    # Import libraries
    import matplotlib.pyplot as plt
    import numpy as np
    from astropy.visualization import astropy_mpl_style
    plt.style.use(astropy_mpl_style)
    np.set_printoptions(precision=4, suppress= True)
    from astropy.io import fits
    
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
    image = image[a1:a2,b1:b2]
     
    
    # Select colorbar range based on percentiles
    pcent = np.nanpercentile(image, [lower_p, upper_p]) 
    
    # Put call to unit converter here
    
    #Do some quick statistics excluding NaNs
    if show_stats == True:
        print('max = ',np.nanmax(image))
        print('min = ',np.nanmin(image))
        print('range =',np.nanmax(image) -np.nanmin(image))
        print('mean = ',np.nanmean(image) )
        print('median = ',np.nanmedian(image) )
        print('mean = ',np.nanmean(image) )
        print('std = ', np.nanstd(image))
        print("{}th percentile = {} & {}th percentile = {}".format(lower_p, pcent[0], upper_p, pcent[1]))
        print("\n=================================================================\n")
    
    
    
    #Make histogram of emission distribution
    
    plt.rcParams["axes.grid"] = False
    if hist_display == True:
        image_vec=image.reshape(image.size) #reshape subimage into a 1D vector
        plt.rcParams['axes.facecolor'] = 'whitesmoke'
        fig1 = plt.figure(figsize=(12, 12))
        plt.hist(image_vec,bins = 50,range=(pcent[0],pcent[1] ),color='tab:blue')
        plt.xlabel('Emission (DN)')
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
        plt.xlabel ('(DN)')
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
        fig3 =plt.figure(figsize=(12,12))
        if wcs_axes == False:
            ax =plt.gca()
            plt.xlabel(r'X (Pixels)')
            plt.ylabel(r'Y (Pixels)')
        else:
            ax = plt.subplot(projection=subwcs)
            plt.xlabel(r'Right Ascension')
            plt.ylabel(r'Declination')
         
        if scaling == "sqrt":
            norm = colors.PowerNorm(gamma=0.5,  vmax= pcent[1], vmin =pcent[0])
            im=ax.imshow(image, cmap=cmap,origin='lower',norm=norm) 
        elif scaling == "log":
            norm = colors.LogNorm(vmax= pcent[1], vmin =pcent[0])
            im=ax.imshow(image, cmap=cmap,origin='lower',norm=norm)
        elif scaling == "linear":
            norm = None
            im=ax.imshow(image, cmap=cmap,origin='lower', vmax= pcent[1], vmin =pcent[0]) 
        im_ratio = image.shape[0]/image.shape[1] 
        cbar=plt.colorbar(im,fraction=0.046*im_ratio, pad=0.05)
        label = "Emission (DN)"
        cbar.set_label(label)
        if grid == True:
            ax.grid(color='white', ls='dotted', lw = 2)
        plt.show()
        print("Figure 3: Image Display \n")
        if save == True:
            fig3.savefig(plotpath + nametag + "_im.png", dpi=300, bbox_inches='tight')
            print("Saving image titled: '" + plotpath + nametag + "_im.png'"  )
            
            
            
            
            
    '''
    Some colormap suggestions:
        inferno, magma, plasma, hot
    '''
    

    # Output subimage fits file
    if save == True:
        header_new=header.copy() # Need to update header
        subshape=subwcs.array_shape
        c_x= subshape[1]//2
        c_y= subshape[0]//2
        c_ra,c_dec=subwcs.pixel_to_world_values(c_x, c_y)
        c_ra,c_dec=float(c_ra),float(c_dec)
        header_new["CRVAL1"] = c_ra
        header_new["CRVAL2"] = c_dec
        header_new["CRPIX1"] = c_x
        header_new["CRPIX2"] = c_y
        header_new["NAXIS1"] = subshape[0]
        header_new["NAXIS2"] = subshape[1]
        header_data_unit_list[0].header = header_new
        header_data_unit_list[0].data = image
        header_data_unit_list.writeto(fitspath + nametag + "_cut.fits", overwrite=True)
        print("Saving fits titled: '" + fitspath + nametag + "_cut.fits'"  )

    print("\n=================================================================\n")