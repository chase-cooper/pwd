from astropy.io import fits
from astropy.table import Table
from astropy.modeling import functional_models
from astropy.convolution import convolve
from astroquery.mast import Observations
from scipy.interpolate import interp1d
from consts import *
import matplotlib.pyplot as plt
import numpy as np
import os,subprocess,time
import urllib.request

import corner,emcee
from multiprocessing import Pool


# Given a scenario whose abundances you wish to fit, and a spectrum to which the scenario is fitted
#   to, find best-fit abundances for a list of metal species.


################################################################################################

def chi2(y,ymodel,yerr):
    """
    Calculate the reduced chi-squared score of a model fit, assuming the model 
    is parameterized by two values ( T_eff, log(g) )
    """
    res = np.sum(np.nan_to_num((y-ymodel)**2/yerr**2,nan=0))
    return res/(len(y)-2)

##########################################################################################

### Functions to open observed and model spectra

def openModel(fnamem:str):
    """
    Read model spectrum
    """
    dat = np.loadtxt(fnamem)
    return dat[:,0],dat[:,1]

def openFile(fname:str):
    """
    Extract the spectrum from a single .fits file.

    I'm not sure that all files I'll be working with share the same format, so
    this function may be edited on the fly.
    """
    with fits.open(fname) as hdu:
        data_fuva   = hdu[1].data[1]
        data_fuvb   = hdu[1].data[0]
        wvln0       = np.hstack((data_fuva[3],data_fuvb[3]))    # Angstroms
        flux0       = np.hstack((data_fuva[4],data_fuvb[4]))    # erg/cm^2/s/A
        sigma0      = np.hstack((data_fuva[5],data_fuvb[5]))    # ^^^^^^^^^^^^

        # Remove zeros in spectrum
        wvln,flux,sigma = [],[],[]
        for i in range(len(wvln0)):
            if flux0[i] != 0:
                wvln.append(wvln0[i])
                flux.append(flux0[i])
                sigma.append(sigma0[i])
        wvln = np.array(wvln)
        flux = np.array(flux)
        sigma= np.array(sigma)
        hdu.close()
    return wvln,flux,sigma

##########################################################################################

### Functions to convolve model spectra with the proper LSF

def lsfParams(fname:str):
    """
    Get relevant parameters to determine which LSF file you need for a given
    FITS file. Borrowed from HST Notebook on LSF convolution
    """
    # Select the primary header
    fuvHeader0 = fits.getheader(fname, ext=0)
    # print(f"For the file {fname}, the relevant parameters are: ")

    # Make a dictionary to store what you find here
    param_dict = {}

    # We want data for the FUV detector, G130M grating at LP3 cenwave 1291
    keywords = ["DETECTOR", "OPT_ELEM", "LIFE_ADJ", "CENWAVE", "DISPTAB"]

    # Print out the relevant values:
    for hdrKeyword in keywords:
        # For DISPTAB
        try:
            # Save the key/value pairs to the dictionary
            value = fuvHeader0[hdrKeyword].split("$")[1]
            # DISPTAB needs the split here
            param_dict[hdrKeyword] = value
        # For other params
        except (IndexError, AttributeError):
            # Save the key/value pairs to the dictionary
            value = fuvHeader0[hdrKeyword]
            param_dict[hdrKeyword] = value

        # Print the key/value pairs
        # print(f"{hdrKeyword} = {value}")

    return param_dict

def fetch_files(det, grating, lpPos, cenwave, disptab):
    """
    Given all the inputs, this will download both
    the LSF and Disptab files to use in the convolution and return their paths.

    Input: 
    det (str): The detector used
    grating (str): Type of grating used
    lpPos (str): Lifetime position used
    cenwave (str): Central wavelength used
    disptab (str): DISPTAB used (will get the path in the function)

    Returns:
    LSF_file_name (str): filename of the new downloaded LSF file
    disptab_path (str): path to the new downloaded disptab file

    Borrowed from HST Notebook on LSF convolution
    """
    datadir = 'hstdata'

    # Link to where all the files live
    COS_site_rootname = (
        "https://www.stsci.edu/files/live/sites/www/files/"
        "home/hst/instrumentation/cos/"
        "performance/spectral-resolution/_documents/"
    )  
    # print(det)

    # Only one file for NUV
    if det == "NUV":
        LSF_file_name = "nuv_model_lsf.dat"

    # FUV files follow a naming pattern
    elif det == "FUV":
        LSF_file_name = f"aa_LSFTable_{grating}_{cenwave}_LP{lpPos}_cn.dat"

    # Where to find file online
    LSF_file_webpath = COS_site_rootname + LSF_file_name
    urllib.request.urlretrieve(
        LSF_file_webpath, f"{datadir}/{LSF_file_name}"
    )
    
    # Where to save file to locally
    # print(f"Downloaded LSF file to {f"{datadir}/{LSF_file_name}"}")

    # And we'll need to get the DISPTAB file as well
    disptab_path = f"{datadir}/{disptab}"
    urllib.request.urlretrieve(
        f"https://hst-crds.stsci.edu/unchecked_get/references/hst/{disptab}",
        disptab_path
    )
    
    # print(f"Downloaded DISPTAB file to {disptab_path}")

    return LSF_file_name, disptab_path

def read_lsf(filename):
    # This is the table of all the LSFs: called "lsf"
    # The first column is a list of the wavelengths corresponding to the line profile, so we set our header accordingly
    # If its an NUV file, header starts 1 line later
    # Borrowed from HST Notebook on LSF convolution
    if "nuv_" in filename:
        ftype = "nuv"
        # print(f"Detector used: {ftype}")
        hs = 1

    # Otherwise, assume its an FUV file
    else:
        ftype = "fuv"
        # print(f"Detector used: {ftype}")

    hs = 0
    lsf = Table.read(filename,
                     format="ascii",
                     header_start=hs)

    # This is the range of each LSF in pixels (for FUV from -160 to +160, inclusive)
    # The middle pixel of the lsf is considered zero; the center is relative zero

    # Integer division to yield whole pixels
    pix = np.arange(len(lsf)) - len(lsf) // 2

    # The column names returned as integers.
    lsf_wvlns = np.array([int(float(k)) for k in lsf.keys()])

    return lsf, pix, lsf_wvlns

def get_disp_params(disptab, cenwave, segment, x=[]):
    """
    Helper function to redefine_lsf(). Reads through a DISPTAB file and gives relevant
    dispersion relationship/wavelength solution over input pixels.

    Parameters:
    disptab (str): Path to your DISPTAB file.
    cenwave (str): Cenwave for calculation of dispersion relationship.
    segment (str): FUVA or FUVB?
    x (list): Range in pixels over which to calculate wvln with dispersion relationship (optional).
    
    Returns:
    disp_coeff (list): Coefficients of the relevant polynomial dispersion relationship
    wavelength (list; if applicable): Wavelengths corresponding to input x pixels 

    Borrowed from HST Notebook on LSF convolution
    """
    with fits.open(disptab) as d:
        wh_disp = np.where(
            (d[1].data["cenwave"] == cenwave)
            & (d[1].data["segment"] == segment)
            & (d[1].data["aperture"] == "PSA")
        )[0]
        # 0 is needed as this returns nested list [[arr]]
        disp_coeff = d[1].data[wh_disp]["COEFF"][0]
    
    # If given a pixel range, build up a polynomial wvln solution pix -> λ
    if len(x):
        wavelength = np.polyval(p=disp_coeff[::-1], x=np.arange(16384))
        return disp_coeff, wavelength
    
    # If x is empty:
    else:
        return disp_coeff

def redefine_lsf(lsf_file, cenwave, disptab, detector="FUV"):
    """
    Helper function to convolve_lsf(). Converts the LSF kernels in the LSF file from a fn(pixel) -> fn(λ)
    which can then be used by convolve_lsf() and re-bins the kernels.

    Parameters:
    lsf_file (str): path to your LSF file
    cenwave (str): Cenwave for calculation of dispersion relationship
    disptab (str): path to your DISPTAB file
    detector (str): FUV or NUV?

    Returns:
    new_lsf (numpy.ndarray): Remapped LSF kernels.
    new_w (numpy.ndarray): New LSF kernel's LSF wavelengths.
    step (float): first order coefficient of the FUVA dispersion relationship; proxy for Δλ/Δpixel.

    Borrowed from HST Notebook on LSF convolution
    """

    if detector == "FUV":
        xfull = np.arange(16384)

        # Read in the dispersion relationship here for the segments
        # FUVA is simple
        disp_coeff_a, wavelength_a = get_disp_params(disptab,
                                                     cenwave,
                                                     "FUVA",
                                                     x=xfull)
        
        # FUVB isn't taken for cenwave 1105, nor 800:
        if (cenwave != 1105) & (cenwave != 800):
            disp_coeff_b, wavelength_b = get_disp_params(disptab,
                                                         cenwave,
                                                         "FUVB",
                                                         x=xfull)
        elif cenwave == 1105:
            # 1105 doesn't have an FUVB so set it to something arbitrary and clearly not real:
            wavelength_b = [-99.0, 0.0]

        # Get the step size info from the FUVA 1st order dispersion coefficient
        step = disp_coeff_a[1]

        # Read in the lsf file
        lsf, pix, w = read_lsf(lsf_file)

        # Take median spacing between original LSF kernels
        deltaw = np.median(np.diff(w))

        # Resamples if the spacing of the original LSF wvlns is too narrow
        if (deltaw < len(pix) * step * 2):  
            # This is all a set up of the bins we want to use
            # The wvln difference between kernels of the new LSF should be about twice their width
            new_deltaw = round(len(pix) * step * 2.0)  

            # nw = Number of LSF wavelengths
            new_nw = (int(round((max(w) - min(w)) / new_deltaw)) + 1)  

            # New version of lsf_wvlns
            new_w = min(w) + np.arange(new_nw) * new_deltaw  

            # Populating the lsf with the proper bins:
            # Empty 2-D array to populate
            new_lsf = np.zeros((len(pix), new_nw)) 

            for i, current_w in enumerate(new_w):
                # Find closest original LSF wavelength to new LSF wavelength
                dist = abs(current_w - w)  
                lsf_index = np.argmin(dist)
                # Column name corresponding to closest orig LSF wvln
                orig_lsf_wvln_key = lsf.keys()[lsf_index]  
                # Assign new LSF wvln the kernel of the closest original lsf wvln
                new_lsf[:, i] = np.array(lsf[orig_lsf_wvln_key])  

        else:
            new_lsf = lsf
            new_w = w

        return new_lsf, new_w, step

    elif detector == "NUV":
        xfull = np.arange(1024)

        # Read in the dispersion relationship here for the segments
        disp_coeff_a, wavelength_a = get_disp_params(disptab,
                                                     cenwave,
                                                     "NUVA",
                                                     x=xfull)
        
        disp_coeff_b, wavelength_b = get_disp_params(disptab,
                                                     cenwave,
                                                     "NUVB",
                                                     x=xfull)
        
        disp_coeff_c, wavelength_c = get_disp_params(disptab,
                                                     cenwave,
                                                     "NUVC",
                                                     x=xfull)

        # Get the step size info from the NUVB 1st order dispersion coefficient
        step = disp_coeff_b[1]

        # Read in the lsf file
        lsf, pix, w = read_lsf(lsf_file)

        # Take median spacing between original LSF kernels
        deltaw = np.median(np.diff(w))

        # This section is a set up of the new bins we want to use:
        # The wvln difference between kernels of the new LSF should be about twice their width
        new_deltaw = round(len(pix) * step * 2.0) 

        # nw = Number of LSF wavelengths
        new_nw = (int(round((max(w) - min(w)) / new_deltaw)) + 1)  

        # New version of lsf_wvlns
        new_w = min(w) + np.arange(new_nw) * new_deltaw  

        # Populating the lsf with the proper bins:
        # Empty 2-D array to populate
        new_lsf = np.zeros((len(pix), new_nw))  

        for i, current_w in enumerate(new_w):
            # Find closest original LSF wavelength to new LSF wavelength
            dist = abs(current_w - w)  
            lsf_index = np.argmin(dist)

            # Column name corresponding to closest orig LSF wvln
            orig_lsf_wvln_key = lsf.keys()[lsf_index]  
            
            # Assign new LSF wvln the kernel of the closest original lsf wvln
            new_lsf[:, i] = np.array(lsf[orig_lsf_wvln_key])  
            
        return new_lsf, new_w, step

def convolve_lsf(wavelength, spec, cenwave, lsf_file, disptab, detector="FUV"):
    """
    Main function; Convolves an input spectrum - i.e. template or STIS spectrum - with the COS LSF.

    Parameters:
    wavelength (list or array): Wavelengths of the spectrum to convolve.
    spec (list or array): Fluxes or intensities of the spectrum to convolve.
    cenwave (str): Cenwave for calculation of dispersion relationship
    lsf_file (str): Path to your LSF file
    disptab (str): Path to your DISPTAB file
    detector (str) : Assumes an FUV detector, but you may specify 'NUV'.

    Returns:
    wave_cos (numpy.ndarray): Wavelengths of convolved spectrum.!Different length from input wvln
    final_spec (numpy.ndarray): New LSF kernel's LSF wavelengths.!Different length from input spec

    Borrowed from HST Notebook on LSF convolution.
    """
    # First calls redefine to get right format of LSF kernels
    new_lsf, new_w, step = redefine_lsf(lsf_file,
                                        cenwave,
                                        disptab,
                                        detector=detector)

    # Sets up new wavelength scale used in the convolution
    nstep = round((max(wavelength) - min(wavelength)) / step) - 1
    wave_cos = min(wavelength) + np.arange(nstep) * step

    # Resampling onto the input spectrum's wavelength scale:
    # Builds up interpolated function from input spectrum
    interp_func = interp1d(wavelength, spec)  

    # Builds interpolated initial spectrum at COS' wavelength scale for convolution
    spec_cos = interp_func(wave_cos)  

    # Initializes final spectrum to the interpolated input spectrum
    final_spec = interp_func(wave_cos)  

    # Loop through the redefined LSF kernels:
    for i, w in enumerate(new_w):  
        # First need to find the boundaries of each kernel's "jurisdiction": where it applies
        # The first and last elements need to be treated separately

        # First kernel:
        if i == 0:  
            diff_wave_left = 500
            diff_wave_right = (new_w[i + 1] - w) / 2.0

        # Last kernel
        elif i == len(new_w) - 1:  
            diff_wave_right = 500
            diff_wave_left = (w - new_w[i - 1]) / 2.0

        # All other kernels
        else:  
            diff_wave_left = (w - new_w[i - 1]) / 2.0
            diff_wave_right = (new_w[i + 1] - w) / 2.0

        # Splitting up the spectrum into slices around the redefined LSF kernel wvlns
        # Will apply the kernel corresponding to that chunk to that region of the spectrum - its "jurisdiction"
        chunk = np.where(
            (wave_cos < w + diff_wave_right) & (wave_cos >= w - diff_wave_left)
        )[0]
        if len(chunk) == 0:
            # Off the edge, go to the next chunk
            continue

        # Selects the current kernel
        current_lsf = new_lsf[:, i]  

        if len(chunk) >= len(
            current_lsf
        ):  # Makes sure that the kernel is smaller than the chunk
            final_spec[chunk] = convolve(
                spec_cos[chunk],
                # Applies the actual convolution
                current_lsf,  
                boundary="extend",
                normalize_kernel=True,
            )

    # Remember, not the same length as input spectrum data!
    return wave_cos, final_spec  

def convolveModelFromFile(fnamem,fname):
    """
    A function that simply puts all the previous steps together. Returns the
    convolved model spectrum
    """
    params  = lsfParams(fname)
    lsf_file,disp_file  = fetch_files(*params.values())
    lsf_file = 'hstdata/'+lsf_file

    wvln_m,spec_m = openModel(fnamem)
    wvln_m,spec_m = convolve_lsf(wvln_m,spec_m,cenwave=params['CENWAVE'],lsf_file=lsf_file,disptab=disp_file,detector=params['DETECTOR'])
    return wvln_m,spec_m

def convolveModel(wvln,flux,fname):
    """
    A function that simply puts all the previous steps together. Returns the
    convolved model spectrum
    """
    params  = lsfParams(fname)
    lsf_file,disp_file  = fetch_files(*params.values())
    lsf_file = 'hstdata/'+lsf_file

    wvln_m,spec_m = wvln,flux
    wvln_m,spec_m = convolve_lsf(wvln_m,spec_m,cenwave=params['CENWAVE'],lsf_file=lsf_file,disptab=disp_file,detector=params['DETECTOR'])
    return wvln_m,spec_m

##########################################################################################

def trilinear_interp(teff,logg,abn,spec_dict):
    """
    Given values for T_eff, log(g), and chondritic silicon abundance, interpolate over the
    model grid to get an approximate spectrum.
    """

    # From values, get adjacent grid points
    t_lo        = teff - (teff % 500)
    t_lo_str    = str(int(t_lo))
    t_hi        = t_lo + 500
    t_hi_str    = str(int(t_hi))
    diff_t      = teff - t_lo

    logg_lo     = logg - (logg % 0.5)
    logg_lo_str = np.format_float_positional(logg_lo,min_digits=1)
    logg_hi     = logg_lo + 0.5
    logg_hi_str = np.format_float_positional(logg_hi,min_digits=1)
    diff_logg   = logg - logg_lo

    abn_lo      = abn - (abn % 0.5)
    abn_lo_str  = np.format_float_positional(abn_lo,min_digits=1)
    abn_hi      = abn_lo + 0.5
    abn_hi_str  = np.format_float_positional(abn_hi,min_digits=1)
    diff_abn    = abn - abn_lo

    # print(t_lo_str,t_hi_str,logg_lo_str,logg_hi_str,abn_lo_str,abn_hi_str)

    # File names of adjacent model spectra
    # root        = 'tlusty/scenarios/'
    root    = 'tlusty/polluted white dwarfs grid/'

    # spec000     = np.genfromtxt(root + f"t{t_lo_str}_g{logg_lo_str}_si{abn_lo_str}"+"/synspec/fort.7")
    # spec001     = np.genfromtxt(root + f"t{t_lo_str}_g{logg_lo_str}_si{abn_hi_str}"+"/synspec/fort.7")
    # spec010     = np.genfromtxt(root + f"t{t_lo_str}_g{logg_hi_str}_si{abn_lo_str}"+"/synspec/fort.7")
    # spec011     = np.genfromtxt(root + f"t{t_lo_str}_g{logg_hi_str}_si{abn_hi_str}"+"/synspec/fort.7")
    # spec100     = np.genfromtxt(root + f"t{t_hi_str}_g{logg_lo_str}_si{abn_lo_str}"+"/synspec/fort.7")
    # spec101     = np.genfromtxt(root + f"t{t_hi_str}_g{logg_lo_str}_si{abn_hi_str}"+"/synspec/fort.7")
    # spec110     = np.genfromtxt(root + f"t{t_hi_str}_g{logg_hi_str}_si{abn_lo_str}"+"/synspec/fort.7")
    # spec111     = np.genfromtxt(root + f"t{t_hi_str}_g{logg_hi_str}_si{abn_hi_str}"+"/synspec/fort.7")
    spec000 = spec_dict[f"t{t_lo_str}_g{logg_lo_str}_si{abn_lo_str}"]
    spec001 = spec_dict[f"t{t_lo_str}_g{logg_lo_str}_si{abn_hi_str}"]
    spec010 = spec_dict[f"t{t_lo_str}_g{logg_hi_str}_si{abn_lo_str}"]
    spec011 = spec_dict[f"t{t_lo_str}_g{logg_hi_str}_si{abn_hi_str}"]
    spec100 = spec_dict[f"t{t_hi_str}_g{logg_lo_str}_si{abn_lo_str}"]
    spec101 = spec_dict[f"t{t_hi_str}_g{logg_lo_str}_si{abn_hi_str}"]
    spec110 = spec_dict[f"t{t_hi_str}_g{logg_hi_str}_si{abn_lo_str}"]
    spec111 = spec_dict[f"t{t_hi_str}_g{logg_hi_str}_si{abn_hi_str}"]

    ### Interpolate. Wavelength points are those of spec000
    wvln        = spec000['wvln']

    # Interpolate over abundances
    spec00x     = spec000['flux'] + np.interp(wvln,spec001['wvln'],spec001['flux']) * diff_abn/0.5
    spec01x     = np.interp(wvln,spec010['wvln'],spec010['flux'])*(1-diff_abn/0.5) + np.interp(wvln,spec011['wvln'],spec011['flux'])*diff_abn/0.5
    spec10x     = np.interp(wvln,spec100['wvln'],spec100['flux'])*(1-diff_abn/0.5) + np.interp(wvln,spec101['wvln'],spec101['flux'])*diff_abn/0.5
    spec11x     = np.interp(wvln,spec110['wvln'],spec110['flux'])*(1-diff_abn/0.5) + np.interp(wvln,spec111['wvln'],spec11['flux'])*diff_abn/0.5

    # Interpolate over log(g)
    spec0xx     = spec00x*(1-diff_logg/0.5) + spec01x*diff_logg/0.5
    spec1xx     = spec10x*(1-diff_logg/0.5) + spec11x*diff_logg/0.5

    # Interpolate over T_eff for final spectrum
    spec_final  = spec0xx*(1-diff_t/500) + spec1xx*diff_t/500

    # Plot result compared to extremes
    if False:
        fig,ax = plt.subplots()

        # ax.plot(spec000[:,0],spec000[:,1],c='salmon',label=t_lo_str+" K, "+logg_lo_str+" dex, "+abn_lo_str+" dex")
        # ax.plot(spec111[:,0],spec111[:,1],c='cornflowerblue',label=t_hi_str+" K, "+logg_hi_str+" dex, "+abn_hi_str+" dex")
        ax.plot(spec000[:,0],spec000[:,1],c='salmon')
        ax.plot(spec001[:,0],spec001[:,1],c='salmon')
        ax.plot(spec010[:,0],spec010[:,1],c='salmon')
        ax.plot(spec011[:,0],spec011[:,1],c='salmon')
        ax.plot(spec100[:,0],spec100[:,1],c='salmon')
        ax.plot(spec101[:,0],spec101[:,1],c='salmon')
        ax.plot(spec110[:,0],spec110[:,1],c='salmon')
        ax.plot(spec111[:,0],spec111[:,1],c='salmon')
        ax.plot(wvln,spec_final,c='black',label="Final")

        ax.set_yscale("log")
        ax.legend()

        plt.show()
        plt.close()

    return wvln,spec_final

def scaleUpModelSpectrum(wvln_model,flux_model,wvln_obs,flux_obs,sigma_obs):
    """
    Use chi2 fitting to find and apply the best scale factor to a modeled spectrum so it fits
    an observed spectrum
    """

    guess = np.average(flux_obs) / np.average(flux_model)
    log_guess = np.log10(guess)

    min_s = 1
    min_score = 1e100
    for s in np.logspace(log_guess-1,log_guess+1,100):
        scaled_flux = flux_model * s
        scaled_flux = np.interp(wvln_obs,wvln_model,scaled_flux)
        score = chi2(flux_obs,scaled_flux,sigma_obs)
        if score < min_score:
            min_score = score
            min_s = s

    return wvln_model,flux_model*min_s

################################################################################################

# Function that reruns SYNSPEC with an altered abundance of a specific abundance
def alter_synspec_abn(scen:str,anum:int,minabn,maxabn,ddex):
    """
    Using a provided scenario's base spectrum, produce a range of spectra with the abundance of
    a specified element 'species' ranging from 'minabn' to 'maxabn' in steps of 'ddex'.

    Inputs:
    - scen (string):    the scenario whose new spectra are being calculated
    - anum (int):       The atomic number of the species whose abundance is changed
    - minabn (int):     the minimum value of log(Sp/H)
    - maxabn (int):     the maximum value of log(Sp/H)
    - ddex (float):     the step in values of log(Sp/H)               
    """

    # Save base spectrum, if not already done
    scendir = SCEN+f"/{scen}/synspec/"
    
    if not 'base.7' in os.listdir(scendir):
        subprocess.run(['cp',f"{scendir}fort.7",f"{scendir}base.7"])

    # Update fort.55 to reflect a change in abundances
    file    = open(f"{scendir}fort.55",'r')
    lines   = file.readlines()
    file.close()

    l = lines[1]
    l = l.split()
    l[-1] = '1\n'
    lines[1] = '\t'.join(l)

    file    = open(f"{scendir}fort.55",'w')
    file.write('\t' + ''.join(lines))
    file.close()

    # Link files
    os.chdir(scendir)
    subprocess.run(['ln','-s','-f',PATH+'/data','data'])
    subprocess.run(['ln','-s','-f',FLAG,'cwd.flag'])

    # Rebuild line list
    subprocess.run(['cp',PATH+'/data/gfATO.dat','fort.19'])
    f = open('fort.19','r')
    lines = f.read().split('\n')
    f.close()

    f = open('fort.19','w')
    for line in lines[:-1]:
        elem = int(line[13:15])
        if int(elem)==anum:
            f.write(line+'\n')
    f.close()

    ### Iteratively run SYNSPEC with new abundances
    abn     = minabn
    while abn <= maxabn:
        
        val     = 10**abn
        valstr  = np.format_float_scientific(val,precision=3)
        abnstr  = np.format_float_positional(abn,precision=1,min_digits=1)
        if f"{anum}_{abnstr}.7" in os.listdir(scendir):
            abn += ddex
            continue

        # Write file
        file = open("fort.56",'w+')
        file.write('1\n')
        file.write(str(anum)+'\t'+valstr)
        file.close()

        # Run SYNSPEC
        # print(f"{anum}:  {valstr}")
        # print("Running SYNSPEC...")
        os.system(f"{SYNEXE} < fort.5 > fort.6")
        subprocess.run(['cp','fort.7',f"{anum}_{abnstr}.7"])
        
        # Modify abundance
        abn += ddex

    os.remove('data')
    os.remove('cwd.flag')
    os.remove('fort.19')
    os.chdir(HOME)

# Function that finds best fit to a line using chi2 minimization, by varying abundance of a species
#   and radial velocity?
def fit_line(spec:str,scen:str,anum:int,minabn,maxabn,ddex,linecen,width):

    # Open and trim observed flux
    wvln_all,flux_all,sigma_all = openFile(spec)
    # wvln    = wvln_all[(wvln_all > linecen-width) and (wvln_all < linecen+width)]
    wvln    = wvln_all[np.abs(wvln_all-linecen)<=width]
    flux    = flux_all[np.abs(wvln_all-linecen)<=width]
    sigma   = sigma_all[np.abs(wvln_all-linecen)<=width]

    vrad_range  = np.arange(-50,50,1)
    abn_range   = np.arange(minabn,maxabn+ddex/2,ddex)
    chi2_grid = np.zeros((len(abn_range),len(vrad_range)))

    abns,probs = [],[]

    for i in range(len(abn_range)):
        abn = abn_range[i]
        abnstr  = round(abn,1)

        # Open model spectrum
        fnamem  = SCEN+f"/{scen}/synspec/{str(anum)}_{abnstr}.7"
        wvln_m_all,flux_m_all   = openModel(fnamem)
        wvln_m_all,flux_m_all   = convolveModel(wvln_m_all,flux_m_all,spec)

        for j in range(len(vrad_range)):
            vrad = vrad_range[j]

            # Doppler shift spectrum
            beta = vrad/3e5
            wvln_m = wvln_m_all * np.sqrt( (1 + beta) / (1 - beta))
            wvln_m = wvln_m[np.abs(wvln_m_all-linecen)<=width]
            flux_m  = flux_m_all[np.abs(wvln_m_all-linecen)<=width]

            # Fit spectra segments
            wvln_m,flux_m = scaleUpModelSpectrum(wvln_m,flux_m,wvln,flux,sigma)
            flux_m = np.interp(wvln,wvln_m,flux_m)
            min_chi2 = chi2(flux,flux_m,sigma)
            chi2_grid[i,j] = min_chi2

    res = np.argmin(chi2_grid)
    min_chi2 = np.min(chi2_grid)
    min_abn_ind,min_vrad_ind = np.unravel_index(res,chi2_grid.shape)
    min_abn     = np.format_float_positional(abn_range[min_abn_ind],precision=1,min_digits=1)
    min_vrad    = vrad_range[min_vrad_ind]
    return min_abn,min_vrad,min_chi2

    # # Retrieve best fit spectrum
    # fname = SCEN+f"/{scen}/synspec/{str(anum)}_{min_abn}.7"
    # wvln_m,flux_m = openModel(fname)
    # wvln_m,flux_m = convolveModel(wvln_m,flux_m,spec)

    # # Doppler shift spectrum
    # beta            = min_vrad/3e5
    # wvln_m          = wvln_m * np.sqrt( (1 + beta) / (1 - beta))

    # wvln_m,flux_m = scaleUpModelSpectrum(wvln_m,flux_m,wvln,flux,sigma)

    # fig,ax = plt.subplots()
    # ax.errorbar(wvln,flux,sigma)
    # ax.plot(wvln_m,flux_m,zorder=100)
    # plt.show()
    # plt.close()

    # fig,ax = plt.subplots()
    # ax.imshow(chi2_grid)
    # plt.show()

################################################################################################

def log_prior(theta):
    """
    Calculate prior likelihood. Assuming flat priors for everything.
    """
    teff,logg,abn,vrad = theta

    # Constraints on parameters
    if (teff < 10000) or (teff > 25000)\
    or (logg < 6) or (logg > 9)\
    or (abn < -9) or (abn > -4)\
    or (vrad < -100) or (vrad > 100):
        return -np.inf

    else:
        return 0

def log_likelihood(theta,spec_dict,specfile,conv_args):
    """
    Calculate the log-likelihood of a forward model with parameters theta explaining the observed
    data in file specfile.
    """

    # Observed data
    x,y,y_err = openFile(specfile)

    # Parameters
    teff,logg,abn,vrad = theta

    # Produce forward model
    wvln_model,flux_model = trilinear_interp(teff,logg,abn,spec_dict)

    # Apply Doppler shift
    beta = vrad/3e5
    wvln_model *= np.sqrt( (1+beta) / (1-beta) )

    # Convolve model
    wvln_model,flux_model = convolve_lsf(wvln_model,flux_model,*conv_args)

    # Scale forward model to observed data
    wvln_model,flux_model = scaleUpModelSpectrum(wvln_model,flux_model,x,y,y_err)
    flux_model = np.interp(x,wvln_model,flux_model)

    score = chi2(y,flux_model,y_err)
    return -0.5*score

def log_prob(theta,spec_dict,specfile,conv_args):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    else:
        return lp + log_likelihood(theta,spec_dict,specfile,conv_args)

def specfit_MCMC():
    pass

################################################################################################

file_obs = 'fits files/gd56fuv.fits'
# x,y,yerr = openFile(file_obs)

# fig,ax = plt.subplots()
# ax.errorbar(x,y,yerr)
# plt.show()

# Carbon

# alter_synspec_abn('t15000_g8.0_si-6.5',
#                   6,
#                   -9,
#                   -4,
#                   0.2)

# fit_line(file_obs,
#          't15000_g8.0_si-6.5',
#          6,
#          -9,
#          -4,
#          0.2,
#          1335,
#          2)

# Silicon

   

# best fit abn = -5.2, best fit vrad = 11 km/s

# fit_line(file_obs,
#          't15000_g8.0_si-6.5',
#          14,
#          -9,
#          -4,
#          0.2,
#          1262.4,
#          5)

# best fit abn = -5.0, best fit vrad = 9 km/s


input()
input()

############################
### MCMC Quarantine zone ###
############################

# # Load spec dict
# start = time.time()
# spec_dict = {}

# for t in np.arange(10000,25001,500):
#     tstr = str(int(t))
#     for g in np.arange(6,9.1,0.5):
#         gstr = np.format_float_positional(g,min_digits=1)
#         for abn in np.arange(-9,-3.9,0.5):
#             abnstr = np.format_float_positional(abn,min_digits=1)
#             scen = f"t{tstr}_g{gstr}_si{abnstr}"
#             dat = np.genfromtxt('tlusty/scenarios/'+scen+'/synspec/fort.7')
#             spec_dict[scen] = {
#                 'wvln': dat[:,0],
#                 'flux': dat[:,1]
#             }
# print(f"Loaded spec dict in {time.time()-start}s")

# # Load necessary arguments for model convolution
# params      = lsfParams(file_obs)
# lsf_file,disp_file  = fetch_files(*params.values())
# lsf_file    = 'hstdata/'+lsf_file
# conv_args   = params['CENWAVE'],lsf_file,disp_file,params['DETECTOR']

# # variables for emcee
# ndim , nwalkers , nstep , nburn = 4 , 10 , 200 , 50

# # initial guess of parameters
# x0  = np.array([15000,8,-9,10])
# pos = [x0*(1 + 1e-2*np.random.randn((4))) for w in range(nwalkers)]

# # backup file
# fn = 'simple emcee.h5'
# backend = emcee.backends.HDFBackend( fn )
# backend.reset( nwalkers , ndim )

# # The big show
# sampler = emcee.EnsembleSampler(nwalkers,ndim,log_prob,args=(spec_dict,file_obs,conv_args),backend=backend)
# sampler.run_mcmc(pos,nstep,progress=True,skip_initial_state_check=True)

# # Get results
# reader  = emcee.backends.HDFBackend(fn)
# nthin   = 2
# samples = reader.get_chain(discard=nburn,thin=nthin,flat=True)

# # Get best fitting model
# lnprob  = reader.get_log_prob(discard=nburn,flat=True,thin=nthin)
# lnpmax  = np.amax(lnprob)
# xmax    = samples[np.where(lnprob==lnpmax)][0]
# print(xmax)

# # Plot results
# fig     = corner.corner(samples,quantiles=[0.16,0.5,0.84],show_titles=True,labels=["T",r"$\log(g)$","[Si/H]",r"$v_{\rm rad}$"],range=[[10000,25000],[6,9],[-9,-4],[-100,100]])
# plt.show()
