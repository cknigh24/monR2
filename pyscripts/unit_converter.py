# unit_converter.py
#+
#-----------------------------------------------------------------------------
#  Make use of astropy units module along with relevant conversion factors
#  to convert objects to appropriate unit types
#
# Currently only converts very specific unit types, will improve later
# 
# However, it may be better to refer to astropy.units API
# for conversions not considered here
#-----------------------------------------------------------------------------
#-




def convert(obj, cf = None, in_unit = "DN/pixel", out_unit = "MJy/sr", pixel_size = 1.375,
            verbose = False):
    
    from astropy import units as u
    
    # Apply a conversion factor to get to physical units
    # i.e. For Wise arrays, conversion factor in [Jy/DN]
    #   W1 cf = 1.9350e-6
    #   W2 cf = 2.7048e-6
    #   w3 cf = 2.9045e-6
    #   w4 cf = 5.2269e-5
    if cf != None:
        obj = obj*cf*(u.Jy) 
        if verbose == True:
            print("conversion factor applied:")
            print(obj)
        
    # Convert from per pixel to per sr
    # Assumes pixel_size is in arcsec
    if  pixel_size != None:
        pixel_size = pixel_size**2*(u.arcsec)**2
        size_sr = pixel_size.to(u.sr)
        if verbose == True:
            print("Pixels converted to sr:")
            print(size_sr)
        
    # Assumes units of Jy/pixel
    if in_unit == "DN/pixel" and out_unit =="MJy/sr":
        obj = obj.to(u.MJy)
        obj = obj/size_sr
        if verbose == True:
            print("Converted to MJy/sr")
            print(obj)
    
    # drop the unit and return the converted value
    return obj.value
    
        
