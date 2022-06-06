# unit_converter.py
#+
#-----------------------------------------------------------------------------
#  Make use of astropy units module along with relevant conversion factors
#  to convert objects to appropriate unit types
#
# Below includes multiple common unit conversions in astronomy 
# 
# However, it may be more useful to refer to astropy.units API
# for conversions not considered here
#-----------------------------------------------------------------------------
#-




def convert(obj, cf = None, in_unit = "DN/pixel", out_unit = "MJy/sr", pixel_size = 1.375,
            wave_0 = None, verbose = False):
    """
    NAME:
    convert
    
    DESCRIPTION:
    Converts from one set of predefined units to another.
    
    INPUT:
    obj = Object witn values of a given unit to be converted. 
    in_unit = String specifying units of input object.
    out_unit = String specifying desired units for object.
    pixel_size = Pixel size of input array in arcseconds. 
    wave_0 = Reference wavelength to use in specific unit conversions. 
    
    INPUT KEYWORD PARAMETERS:
    verbose = Boolean value to print output message of final units.
    
    
    """
    
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
    
    # Same as above but skips conversion from DN         
    if in_unit == "Jy/pixel" and out_unit =="MJy/sr" and cf == None:
        obj = obj.to(u.MJy)
        obj = obj/size_sr
        if verbose == True:
            print("Converted to MJy/sr")
            print(obj)
            
       
    if in_unit == "MJy/sr" and out_unit =="Jy/pixel" and cf == None:
        obj = obj.to(u.Jy)
        obj = obj*size_sr
        if verbose == True:
            print("Converted to Jy/pixel")
            print(obj)  
            
    if in_unit == "ergs cm^-2 s^-1 A^-1" and out_unit =="W m^-2 micron^-1" and cf == None:        
        obj = obj*10.0
        if verbose == True:
            print("Converted to W m^-2 micron^-1")
            print(obj)
            
    if in_unit == "W m^-2 micron^-1" and out_unit =="ergs cm^-2 s^-1 A^-1" and cf == None:        
        obj = obj/10.0
        if verbose == True:
            print("Converted to ergs cm^-2 s^-1 A^-1")
            print(obj)
            
    if in_unit == "W m^-2 micron^-1" and out_unit =="Jy" and wave_0 != None:         
        obj = obj*1e26*3.3e-15*wave_0**2
        if verbose == True:
            print("Converted to Jy")
            print(obj)
            
    if in_unit == "Jy" and out_unit =="W m^-2 micron^-1" and wave_0 != None:         
        obj = obj/1e26/3.3e-15/wave_0**2
        if verbose == True:
            print("Converted to W m^-2 micron^-1")
            print(obj)
            
    if in_unit == "W m^-2/pixel" and out_unit =="W m^-2 micron^-1/sr" and cf == None:         
        obj = obj/size_sr
        if verbose == True:
            print("Converted to W m^-2 micron^-1/sr")
            print(obj)
            
    if in_unit == "W m^-2sr" and out_unit =="W m^-2 micron^-1/pixel" and cf == None:         
        obj = obj*size_sr
        if verbose == True:
            print("Converted to W m^-2 micron^-1/pixel")
            print(obj)        
        
    # drop the unit and return the converted value
    return obj.value
    
        
