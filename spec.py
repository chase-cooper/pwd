from astropy.io import fits
from consts import *
import matplotlib.pyplot as plt
import numpy as np
import os,subprocess


# Given a scenario whose abundances you wish to fit, and a spectrum to which the scenario is fitted
#   to, find best-fit abundances for a list of metal species.

# For testing, I'm using a scenario, '944', with T = 17,000K, log(g) = 8, and a chondritic [Si/H]
#   of -5.5. Values are approximately those derived by Rogers+24a.
specf   = 'fits files/gd56fuv.fits'
scen    = 'gd56'

dat     = np.genfromtxt(SCEN+f"/{scen}/synspec/fort.7")
wvln    = dat[:,0]
flux    = dat[:,1]

################################################################################################

def chi2(y,ymodel,yerr):
    """
    Calculate the reduced chi-squared score of a model fit, assuming the model 
    is parameterized by two values ( T_eff, log(g) )
    """
    res = np.sum(np.nan_to_num((y-ymodel)**2/yerr**2,nan=0))
    return res/(len(y)-2)

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

    ### Iteratively run SYNSPEC with new abundances

    abn     = minabn
    while abn <= maxabn:
        
        val     = 10**abn
        valstr  = np.format_float_scientific(val,precision=3)
        if f"{anum}_{str(round(abn,1))}.7" in os.listdir(scendir):
            abn += ddex
            continue

        # Write file
        file = open(f"{scendir}fort.56",'w+')
        file.write('1\n')
        file.write(str(anum)+'\t'+valstr)
        file.close()

        # Run SYNSPEC
        os.chdir(scendir)
        print(f"{anum}:  {valstr}")
        print("Running SYNSPEC...")
        os.system(f"{SYNEXE} < fort.5 > fort.6")
        subprocess.run(['cp','fort.7',f"{anum}_{str(round(abn,1))}.7"])
        os.chdir(HOME)

        # Modify abundance
        abn += ddex

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
    print(f"For the file {fname}, the relevant parameters are: ")

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
        print(f"{hdrKeyword} = {value}")

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
    print(det)

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
    print(f"Downloaded LSF file to {f"{datadir}/{LSF_file_name}"}")

    # And we'll need to get the DISPTAB file as well
    disptab_path = f"{datadir}/{disptab}"
    urllib.request.urlretrieve(
        f"https://hst-crds.stsci.edu/unchecked_get/references/hst/{disptab}",
        disptab_path
    )
    
    print(f"Downloaded DISPTAB file to {disptab_path}")

    return LSF_file_name, disptab_path

def read_lsf(filename):
    # This is the table of all the LSFs: called "lsf"
    # The first column is a list of the wavelengths corresponding to the line profile, so we set our header accordingly
    # If its an NUV file, header starts 1 line later
    # Borrowed from HST Notebook on LSF convolution
    if "nuv_" in filename:
        ftype = "nuv"
        print(f"Detector used: {ftype}")
        hs = 1

    # Otherwise, assume its an FUV file
    else:
        ftype = "fuv"
        print(f"Detector used: {ftype}")

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

def convolveModel(fnamem,fname):
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

##########################################################################################

### Function to scale a modeled spectrum up to an observed spectrum, fitting for scale factor
###     and radial velocity
def scaleModel(wvln,flux,sigma,wvln_m,flux_m):
    """
    Fit an observed spectrum (or stacked spectrum) with a modeled one, 
    adjusting for a scale factor and Doppler shift.
    """
    
    # Get model data (usually my model data, but others work fine)
    model_wvln  = wvln_m
    model_flux  = np.interp(wvln,model_wvln,flux_m,left=0,right=0) # Interpolate at observed spectrum points

    # Rough scale fit to get approximate scale factor
    min_s   = 1
    min_score = 1e20
    for s in np.logspace(-30,0,num=5000):
        scaled_model_flux   = model_flux*s
        # score = np.sum(np.abs(flux-scaled_model_flux))
        score   = chi2(flux,scaled_model_flux,sigma)
        if score<min_score:
            min_s=s
            min_score=score
    
    # Proper chi2 grid; doppler shift fitting is pretty jank atm so don't trust it
    n = 50
    v_grid      = np.linspace(-100,100,n)
    s_grid      = np.logspace(np.log10(min_s)-1,np.log10(min_s)+1,n)
    chi2_grid   = np.zeros((n,n))
    for i in range(n):
        v = v_grid[i]
        beta = v/300_000
        for j in range(n):
            s = s_grid[j]

            # Scale model flux
            scaled_model_wvln   = wvln*np.sqrt((1-beta)/(1+beta))
            scaled_model_flux   = model_flux*s

            # Re-interpolate to observed spectrum points
            scaled_model_flux   = np.interp(wvln,scaled_model_wvln,scaled_model_flux)
            score = chi2(flux,scaled_model_flux,sigma)
            chi2_grid[i,j] = score
    
    # Fetch chi2 minimum
    ind_v,ind_s = np.unravel_index(np.argmin(chi2_grid),shape=(n,n))
    min_chi2    = chi2_grid[ind_v][ind_s]
    min_v       = v_grid[ind_v]
    min_s       = s_grid[ind_s]
    print(f"Chi-2 minimum:              {min_chi2}")
    print(f"Best-fit flux scale factor: {min_s}")
    print(f"Best-fit Doppler shift:     {min_v}km/s")

    # Rescale model flux for plotting
    beta        = min_v/300_000 # again, jank; trusting this leads to offset line features
    scaled_model_wvln   = wvln*np.sqrt((1-beta)/(1+beta))
    scaled_model_flux   = model_flux*min_s
    # Re-interpolate to observed spectrum points (not strictly necessary)
    # scaled_model_flux   = np.interp(wvln,scaled_model_wvln,scaled_model_flux)

    # # Plot that shit
    # fig,ax = plt.subplots()
    # ax.plot(wvln,flux)
    # ax.plot(scaled_model_wvln,scaled_model_flux)
    # plt.show()

    return min_chi2,min_v,min_s
  
# Function that finds best fit to a line using chi2 minimization, by varying abundance of a species
def fit_line(spec:str,scen:str,anum:int,minabn,maxabn,ddex,linecen,width):

    # Open and trim observed flux
    wvln_all,flux_all,sigma_all = openFile(spec)
    # wvln    = wvln_all[(wvln_all > linecen-width) and (wvln_all < linecen+width)]
    wvln    = wvln_all[np.abs(wvln_all-linecen)<=width]
    flux    = flux_all[np.abs(wvln_all-linecen)<=width]
    sigma   = sigma_all[np.abs(wvln_all-linecen)<=width]

    abns,probs = [],[]

    abn     = minabn
    while abn <= maxabn:
        abnstr  = round(abn,1)

        # Open and trim model spectrum
        fnamem  = SCEN+f"/{scen}/synspec/{str(anum)}_{abnstr}.7"
        wvln_m_all,flux_m_all   = openModel(fnamem)
        wvln_m  = wvln_m_all[np.abs(wvln_m_all-linecen)<=width]
        flux_m  = flux_m_all[np.abs(wvln_m_all-linecen)<=width]

        # Fit spectra segments
        min_chi2,min_v,min_s = scaleModel(wvln,flux,sigma,wvln_m,flux_m)
        prob    = np.exp(-0.5*min_chi2)
        abns.append(abn)
        probs.append(prob)

        abn += ddex

    ### Derive best-fit value and errorbars
    bestfit = abns[np.argmax(probs)]
    print(bestfit)
    
    pdf     = probs/np.sum(probs)
    cdf     = np.zeros_like(probs)
    for i in range(len(probs)):
        cdf[i] = np.sum(pdf[:i+1])

    fit = np.polyfit(abns,cdf,deg=5)
    fit = np.poly1d(fit)

    fig,ax = plt.subplots()
    ax.plot(abns,cdf)
    ax.plot(abns,fit(abns))
    plt.show()


    
alter_synspec_abn(scen=scen,anum=16,minabn=-7,maxabn=-5,ddex=0.2)

fit_line(specf,scen,anum=14,minabn=-6.5,maxabn=-4.5,ddex=0.1,linecen=1265,width=1)
