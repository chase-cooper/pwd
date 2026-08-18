### Catalogs
from astroquery.gaia import Gaia
from astroquery.ipac.irsa import Irsa
from astroquery.mast import Catalogs
from astroquery.mast import Mast
from astroquery.sdss import SDSS
from astroquery.simbad import Simbad
mast = Mast()

### Dust maps
from dustmaps.sfd import SFDQuery
from dustmaps.config import config
config['data_dir'] = '/Users/cooper/Documents/Research/pwd/dustmaps'
sfd     = SFDQuery()

### MCMC Fitting
import corner,emcee
from multiprocessing import Pool

### Other
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import LinearNDInterpolator
import os,requests,shutil,subprocess,sys,time
from consts import SCEN
from utils import newSpeciesDict

##################################################################################################################

# Constants
c       = 3e8
pc2m    = 3.086e16
pc2cm   = 3.086e18
grav    = 6.67e-11  # [m^3 kg^-1 s^-2]
msol    = 1.989e30  # [kg]
R_earth = 6.371e6   # [m]

##################################################################################################################

### Utility functions

def chi2(y,ymodel,yerr):
    """
    Calculate the reduced chi-squared score of a model fit, assuming the model 
    is parameterized by two values ( T_eff, log(g) )
    """
    res = (y - ymodel)**2 / yerr**2
    # print(res)
    return np.sum(res)/(len(y)-2)

def ABmag2flux(mag,wvln):
    """
    Calculate an object's flux given its AB magnitude. Flux can either be in units of
    Janskys (default) or W/m^2/nm. In the case of the latter, an effective wavelength
    is needed.
    """

    c = 3e8
    fluxJanskys = 10**(3.56 - 0.4*mag)                  # Jy = 10^-26 W/m^2/Hz
    # flux        = 1e-26*fluxJanskys*c/(wvln*wvln)   # W/m^2/nm
    return fluxJanskys
ABmag2flux = np.vectorize(ABmag2flux)

##################################################################################################################

### Catalog query functions

def query_simbad(name:str):

    ### Get designations from SIMBAD
    """
    Query simbad to get the Gaia DR3 designation for a given object.
    """

    name_gaia    = None

    # query SIMBAD
    res = Simbad.query_objectids(name)
    for des in res[:]:
        des = des[0]
        if 'Gaia DR3' in des:
            name_gaia = des

    # if no Gaia counterpart is found, quit
    if name_gaia==None:
        print("No Gaia ID found. Exiting...")
        exit(0)

    name_gaia = name_gaia.split()[-1]
    return name_gaia

def query_gaia(name:str,verbose:bool=False):

    ### Query Gaia main table for parallax, position, and photometry
    job = Gaia.launch_job("SELECT TOP 100 source_id, parallax, parallax_error, phot_g_mean_mag, phot_g_mean_flux, phot_g_mean_flux_error, phot_bp_mean_mag, phot_bp_mean_flux, phot_bp_mean_flux_error, phot_rp_mean_mag, phot_rp_mean_flux, phot_rp_mean_flux_error, b "
                                f"from gaiadr3.gaia_source where source_id = {name}",
                                dump_to_file=False)
    res = job.get_results()

    # Parallax and position
    plx         = res['parallax'][0]
    plx_err     = res['parallax_error'][0]
    b           = res['b'][0]   # Galactic latitude

    # Distance [pc]
    dist        = 1000/plx
    dist_err    = [dist-1000/(plx+plx_err),dist-1000/(plx-plx_err)]

    # Height above galactic plane
    z = dist*np.sin(np.deg2rad(b))
    z = np.abs(z)

    # Extinction fraction
    if dist < 100:
        ext_frac = 0
    else:
        ext_frac    = min( (dist-100) / (250/np.abs(np.sin(np.deg2rad(b))) - 100), 1)

    ### Query SDSS crossmatch table for SDSS neighbour
    job_sdss    = Gaia.launch_job_async(f"SELECT TOP 10 source_id, original_ext_source_id, angular_distance FROM gaiadr3.sdssdr13_best_neighbour WHERE source_id = {name}",
                                dump_to_file=False)
    sdss    = job_sdss.get_results()
    try:
        id_sdss = sdss['original_ext_source_id'][0]
    except IndexError:
        print("No SDSS DR13 neighbour found")
        id_sdss = None

    ### Query SkyMapper crossmatch table for SkyMapper neighbour
    job_skym    = Gaia.launch_job(f"SELECT TOP 10 source_id, original_ext_source_id, angular_distance FROM gaiadr3.skymapperdr2_best_neighbour WHERE source_id = {name}",
                                dump_to_file=False)
    skym    = job_skym.get_results()
    try:
        id_skym = skym['original_ext_source_id'][0]
    except IndexError:
        print("No SkyMapper DR2 neighbour found")
        id_skym = None

    ### Query PanSTARRS crossmatch table for PanSTARRS neighbour
    job_pans    = Gaia.launch_job(f"SELECT TOP 10 source_id, original_ext_source_id, angular_distance FROM gaiadr3.panstarrs1_best_neighbour WHERE source_id = {name}",
                                dump_to_file=False)
    pans    = job_pans.get_results()
    try:
        id_pans = pans['original_ext_source_id'][0]
    except IndexError:
        print("No PanSTARRS DR1 neighbour found")
        id_pans = None

    ### Query 2MASS crossmatch table for 2MASS neighbour
    job_tmas    = Gaia.launch_job(f"SELECT TOP 10 source_id, original_ext_source_id, angular_distance FROM gaiadr3.tmass_psc_xsc_best_neighbour WHERE source_id = {name}",
                                  dump_to_file=False)
    tmas    = job_tmas.get_results()
    try:
        id_2mas = tmas['original_ext_source_id'][0]
    except IndexError:
        print("No 2MASS neighbour found")
        id_2mas = None

    if verbose:
        print("*"*100)
        print("Parallax:        "+str(plx)+" +\\-"+str(plx_err)+" arcsec")
        print("Distance:        "+str(dist)+" +"+str(dist_err[0])+"\\-"+str(np.abs(dist_err[1]))+" pc")
        print("Gal. Latitude:   "+str(b)+" degrees")
        print("|z|:             "+str(z)+" pc")
        print("Ext. fraction:   "+str(ext_frac))
        print("\nNeighbours:")
        print(" - SDSS:         "+str(id_sdss))
        print(" - PanSTARRS:    "+str(id_pans))
        print(" - SkyMapper:    "+str(id_skym))
        print(" - 2MASS:        "+str(id_2mas))

    return plx,plx_err,ext_frac,id_sdss,id_pans,id_skym,id_2mas

def query_sdss(name:int,pos):

    # Query
    fields  = ['objID','u','g','r','i','z','err_u','err_g','err_r','err_i','err_z','extinction_u','extinction_g','extinction_r','extinction_i','extinction_z']
    xid     = SDSS.query_region(pos, radius='5 arcsec', spectro=False, fields=fields)
    
    # Get data
    for i in xid:
        # Check each row to see if it contains the SDSS designation from Gaia cross-matching
        if i['objID']==name:
            res = np.array(list(i))
    
            # Apply SDSS-to-AB mag corrections
            res[1] -= 0.04
            res[4] += 0.015
            res[5] += 0.03
            
            return res[1:6],res[6:11],res[11:]  # AB mags, errors, extinctions
    
    return None

def query_panstarrs(name:int,pos):
    try:
        catalog_data = Catalogs.query_region(pos,
                                        radius=0.01,
                                        catalog="Panstarrs",
                                        data_release="dr1",
                                        table="mean")
        
        # Find line whose object ID is the same one returned by Gaia cross-matching
        for line in catalog_data:
            if line['objID']==name:
                pans_mags   = np.array([line['gMeanApMag'],line['rMeanApMag'],line['iMeanApMag'],line['zMeanApMag'],line['yMeanApMag']])
                pans_errs   = np.array([line['gMeanApMagErr'],line['rMeanApMagErr'],line['iMeanApMagErr'],line['zMeanApMagErr'],line['yMeanApMagErr']])
                
                # Extinctions
                pans_rvals       = np.array([3.384,2.482,1.93,1.551,1.258])        # values of R, from Laura
                # pans_rvals       = np.array([3.172,2.271,1.682,1.322,1.087])       # values of R assuming R_V=3.1, Schlafly & Finkbeiner 2011

                return pans_mags,pans_errs,pans_rvals
    
    except:
        print("\nFailed to query PanSTARRS for object "+str(name_pans))
        return np.array([]),np.array([])
    
def query_skymapper(name,pos):

    # Query parameters
    ra,dec = pos.to_string(style='decimal').split()
    params = {
        'RA':   ra,
        'DEC':  dec,
        'SR':   '0.1',
        'VERB': '3',
        'RESPONSEFORMAT':   'CSV'
    }

    # make query
    url = 'https://skymapper.anu.edu.au/sm-cone/public/query?'
    for key in params.keys():
        url  += key + '=' + params[key] + '&'
    url = url[:-1]
    req     = requests.get(url)

    # save data
    txt     = req.text
    txt     = txt.split('\n')
    for line in txt:
        line = line.split(',')

        # Find line with correct name
        if line[1] == str(name):
            
            # Important indices:
            #   1  --> object id
            #   43 --> u mag
            #   57 --> g mag
            #   64 --> r mag
            #   71 --> i mag
            #   78 --> z mag
            #   95 --> closest Gaia neighbor

            skym_mags   = np.array([line[43],line[57],line[64],line[71],line[78]],dtype=np.float32)
            skym_errs   = np.array([line[44],line[58],line[65],line[72],line[79]],dtype=np.float32)

            # R-values
            skym_rvals  = np.array([5,3.4,2.5,1.8,1.4])     # from Laura
            
            return skym_mags,skym_errs,skym_rvals
    
    # If object isn't found
    return np.array([]),np.array([]),np.array([])

def query_2mas(name,pos):
    try:
        catalog_data = Irsa.query_region(coordinates=pos,
                                         catalog='fp_psc')
        for line in catalog_data:
            if line['designation']==name:
                tmas_mags = np.array([line['j_m'],line['h_m'],line['k_m']])
                tmas_errs = np.array([line['j_cmsig'],line['h_cmsig'],line['k_cmsig']])
                tmas_rvals = np.array([0.72,0.46,0.306])    # from Yuan, Liu, and Xiang (2013)

        return tmas_mags,tmas_errs,tmas_rvals
    except:
        print("\nFailed to query 2MASS for object "+name)
        return np.array([]),np.array([]),np.array([])

##################################################################################################################

### Retrieve all extinction-corrected photometric fluxes

def photometry(name:str):

    ########################
    ###  Object queries  ###
    ########################

    ### Get object position (RA, DEC) and get line-of-sight E(B-V)
    pos     = mast.resolve_object(name,resolver=None)
    ebv     = sfd(pos)
    print("E(B-V):      "+str(ebv))


    ### Query SIMBAD for object's designation in Gaia DR3. If no such designation exists, exit.
    print("Querying SIMBAD...",end='')
    name_gaia   = query_simbad(name)
    print("Done!")
    print('\n' + '*'*100 + '\n')


    ### Query Gaia for parallax and galactic latitude to derive distances and height above/below
    ###     the galactic plane |z|; use these to determine fraction of extinction to apply 
    ###     following procedure in Harris (2006).
    ###
    ###     Additionally, query Gaia's crossmatching tables for object designations in SDSS DR12,
    ###     PanSTARRS DR1, and SkyMapper DR2 surveys.
    print("Querying Gaia DR3...")
    name_gaia   =   name.split()[-1]
    plx,plx_err,ext_frac,name_sdss,name_pans,name_skymapper,name_2mass = query_gaia(name_gaia,verbose=True)
    print("\nFinished querying Gaia DR3")
    print("\n" + '*'*100 + '\n')


    ### Query SDSS DR19 if an SDSS object ID is returned from Gaia query. Returns object mags,
    ###     errors, and extinctions for each ugriz passband.
    if name_sdss:
        print("Querying SDSS...")
        sdss_res    = query_sdss(name_sdss,pos)
        if sdss_res!=None:
            sdss_mags,sdss_errs,sdss_exts   = sdss_res
            sdss_mags   = np.array(sdss_mags)
            sdss_errs   = np.array(sdss_errs)
            sdss_exts   = np.array(sdss_exts)

            # Correct SDSS mags for reddening
            sdss_mags   -= sdss_exts*ext_frac

            print(sdss_mags)
            print(sdss_errs)
            print(sdss_exts)
        
        else:
            # If no SDSS object ID is found, create empty arrays for SDSS quantities
            sdss_mags = np.array([])
            sdss_errs = np.array([])
            sdss_exts = np.array([])

        print("Done!")
        print('\n' + '*'*100 + '\n')
    
    else:
        # If no SDSS object ID is found, create empty arrays for SDSS quantities
        sdss_mags = np.array([])
        sdss_errs = np.array([])
        sdss_exts = np.array([])


    ### Query PanSTARRS DR1 if a PanSTARRS object ID is returned from Gaia query. Returns object
    ###     magnitudes, errors, and R-values for each grizy passband. R-values are provided in the
    ###     function call.
    if name_pans:
        print("Querying PanSTARRS...")
        pans_mags,pans_errs,pans_rvals = query_panstarrs(name_pans,pos) # already numpy arrays

        # Extinction correction
        pans_exts   = pans_rvals * ebv
        pans_mags   -= pans_exts * ext_frac

        print(pans_mags)
        print(pans_errs)
        print(pans_exts)

        print("Done!")
        print('\n' + '*'*100 + '\n')

    else:
        # If no PanSTARRS object ID is found, create empty arrays for PanSTARRS quantities
        pans_mags = np.array([])
        pans_errs = np.array([])
        pans_exts = np.array([])
    

    ### Query SkyMapper DR2 if a SkyMapper object ID is returned from the Gaia query. Returns
    ###     object magnitudes, errors, and R-values for each ugriz passband. R-values are provided
    ###     in the function call.
    if name_skymapper:
        print("Querying SkyMapper DR2...")

        skym_mags,skym_errs,skym_rvals  = query_skymapper(name_skymapper,pos)
        skym_exts   = skym_rvals * ebv
        skym_mags   -= skym_exts * ext_frac

        print(skym_mags)
        print(skym_errs)
        print(skym_exts)

        print("Done!")
        print('\n' + '*'*100 + '\n')
    
    else:
        # If no SkyMapper object ID is found, create empty arrays for SkyMapper quantities
        skym_mags = np.array([])
        skym_errs = np.array([])
        skym_exts = np.array([])

    ### Query 2MASS if a 2MASS object ID is returned from the Gaia query. Returns object
    ###     magnitudes, errors, and R-values for each jhk passband. R-values are provided in the
    ###     function call.
    if name_2mass:
        print("Querying 2MASS...")

        tmas_mags,tmas_errs,tmas_rvals = query_2mas(name_2mass,pos)
        tmas_exts = tmas_rvals*ebv
        tmas_mags -= tmas_exts

        print(tmas_mags)
        print(tmas_errs)
        print(tmas_exts)
        print("Done!")
        print('\n' + '*'*100 + '\n')

    else:
        tmas_mags = np.array([])
        tmas_errs = np.array([])
        tmas_exts = np.array([])

    ##################################################################################################################

    ### Formatting of data and conversion of magnitudes to fluxes in Jy
    ###
    ### This kinda sucks but whatever it works

    # Filter effective wavelengths [nm] from SVO
    sdss_eff_wvlns  = np.array([360.8,467.18,614.11,745.79,892.28]) * 1e-9
    pans_eff_wvlns  = np.array([481.02,615.55,750.3,866.84,961.36]) * 1e-9
    skym_eff_wvlns  = np.array([350.02,501.6,607.69,773.28,912.03]) * 1e-9
    tmas_eff_wvlns  = np.array([1235,1662,2159]) * 1e-9

    # Mags to fluxes and their errors
    all_wvlns       = np.array([])
    all_fluxes      = np.array([])
    all_fluxerrs_hi = np.array([])
    all_fluxerrs_lo = np.array([])
    
    if np.any(sdss_mags):   
        sdss_fluxes     = ABmag2flux(sdss_mags,sdss_eff_wvlns)
        sdss_fluxes_hi  = ABmag2flux(sdss_mags - sdss_errs,sdss_eff_wvlns)
        sdss_fluxes_lo  = ABmag2flux(sdss_mags + sdss_errs,sdss_eff_wvlns)
        sdss_fluxerr_hi = sdss_fluxes_hi - sdss_fluxes
        sdss_fluxerr_lo = sdss_fluxes - sdss_fluxes_lo
    else:
        sdss_fluxes     = np.array([0,0,0,0,0])
        sdss_fluxerr_hi = np.array([1e10,1e10,1e10,1e10,1e10])

    if np.any(pans_mags):
        pans_fluxes     = ABmag2flux(pans_mags,pans_eff_wvlns)
        pans_fluxes_hi  = ABmag2flux(pans_mags - pans_errs,pans_eff_wvlns)
        pans_fluxes_lo  = ABmag2flux(pans_mags + pans_errs,pans_eff_wvlns)
        pans_fluxerr_hi = pans_fluxes_hi - pans_fluxes
        pans_fluxerr_lo = pans_fluxes - pans_fluxes_lo
    else:
        pans_fluxes     = np.array([0,0,0,0,0])
        pans_fluxerr_hi = np.array([1e10,1e10,1e10,1e10,1e10])

    if np.any(skym_mags):
        skym_fluxes     = ABmag2flux(skym_mags,skym_eff_wvlns)
        skym_fluxes_hi  = ABmag2flux(skym_mags - skym_errs,skym_eff_wvlns)
        skym_fluxes_lo  = ABmag2flux(skym_mags + skym_errs,skym_eff_wvlns)
        skym_fluxerr_hi = skym_fluxes_hi - skym_fluxes
        skym_fluxerr_lo = skym_fluxes - skym_fluxes_lo
    else:
        skym_fluxes     = np.array([0,0,0,0,0])
        skym_fluxerr_hi = np.array([1e10,1e10,1e10,1e10,1e10])

    if np.any(tmas_mags):   
        # NOTE: 2MASS mags are not absolute (AB) magnitudes; fluxes are calculated using zero-
        #   point fluxes in Cohen+03
        tmas_fluxes     = np.array([1594,1024,666.7]) * 10**(-0.4*tmas_mags)
        tmas_fluxes_lo  = np.array([1594,1024,666.7]) * 10**(-0.4 * (tmas_mags+tmas_errs) )
        tmas_fluxes_hi  = np.array([1594,1024,666.7]) * 10**(-0.4 * (tmas_mags-tmas_errs) )
        # tmas_fluxes     = ABmag2flux(tmas_mags,tmas_eff_wvlns)
        # tmas_fluxes_hi  = ABmag2flux(tmas_mags - tmas_errs,tmas_eff_wvlns)
        # tmas_fluxes_lo  = ABmag2flux(tmas_mags + tmas_errs,tmas_eff_wvlns)
        tmas_fluxerr_hi = tmas_fluxes_hi - tmas_fluxes
        tmas_fluxerr_lo = tmas_fluxes - tmas_fluxes_lo
    else:
        tmas_fluxes     = np.array([0,0,0,0,0])
        tmas_fluxerr_hi = np.array([1e10,1e10,1e10,1e10,1e10])

    all_wvlns       = np.concat([sdss_eff_wvlns,pans_eff_wvlns,skym_eff_wvlns,tmas_eff_wvlns])
    all_fluxes      = np.concat([sdss_fluxes,pans_fluxes,skym_fluxes,tmas_fluxes])
    all_fluxerrs_hi = np.concat([sdss_fluxerr_hi,pans_fluxerr_hi,skym_fluxerr_hi,tmas_fluxerr_hi])

    return all_wvlns,all_fluxes,all_fluxerrs_hi,plx,plx_err

##################################################################################################################

### Convolve a model spectrum with a filter response function

def convolveModelWithBandpass(mwvln,mflux,band_fname):
    """
    Convolve a TLUSTY/SYNSPEC spectrum (at location 'model_fname') with a filter response function
    (at location 'band_fname').
    """

    # Speed of light in AA/s
    cAA = 3e18

    # Datatype change
    mwvln       = np.array(mwvln)
    mflux       = np.array(mflux)

    model_wvln0 = mwvln[1:]
    model_flux0 = mflux[1:]
    model_dlam0 = mwvln[1:] - mwvln[:-1]

    # TLUSTY/SYNSPEC spectra are Eddington fluxes, so require a factor of 4*pi to get fluxes
    model_flux0 *= 4*np.pi

    filt    = np.genfromtxt(band_fname)
    filt_wvlns  = filt[:,0]
    filt_resp   = filt[:,1]

    # Interpolate response to model wavelengths
    model_response  = np.interp(model_wvln0,filt_wvlns,filt_resp,left=0,right=0)

    # Resulting flux in [erg/s/cm^2/Hz]
    conv_flux   = (1/cAA) * np.sum(model_response*model_flux0*model_wvln0*model_dlam0) / np.sum(model_response * model_dlam0 / model_wvln0)
    conv_flux   *= 1e23     # Convert to Jy
    return conv_flux

##################################################################################################################

### Forward model functions

def loadSpectra():
    # Spectrum dictionary
    start = time.time()
    print("Loading spectra...",end='')

    spec_dict = {}
    for t in np.arange(10000,25001,500):
        tstr = str(t)
        for g in np.arange(6,9.1,0.5):
            gstr = np.format_float_positional(g,precision=1,min_digits=1)

            scen = f"t{tstr}_g{gstr}"
            spec_dict[scen] = {}

            dat = np.genfromtxt(f'tlusty/hydrogengrid/{scen}.spec')
            fwvln = dat[:,0]
            fflux = dat[:,1]
            spec_dict[scen] = {
                'wvln':fwvln,
                'flux':fflux
            }

    end = time.time()
    print(f"Loaded spectra in {round(end-start,2)}s")
    return spec_dict

def interp_spectrum_dict(t,logg,spec_dict):
    """
    Given a temperature (anything) and log(g) (to 2 decimal places), linearly interpolate between
    adjacent spectra to produce a new spectrum. Instead of reading in data from files, the data
    are accessed through a dictionary where all the data has been preloaded.
    """
    t_lo        = t - (t % 500)
    t_lo_str    = str(int(t_lo))
    t_hi        = t_lo + 500
    t_hi_str    = str(int(t_hi))
    diff_t      = t - t_lo

    logg_lo     = logg - (logg % 0.5)
    logg_lo_str = np.format_float_positional(logg_lo,min_digits=1)
    logg_hi     = logg_lo + 0.5
    logg_hi_str = np.format_float_positional(logg_hi,min_digits=1)
    diff_logg   = logg - logg_lo

    # Bilinear interpolation over model grid
    spec00  = spec_dict[f"t{t_lo_str}_g{logg_lo_str}"]
    spec01  = spec_dict[f"t{t_lo_str}_g{logg_hi_str}"]
    spec10  = spec_dict[f"t{t_hi_str}_g{logg_lo_str}"]
    spec11  = spec_dict[f"t{t_hi_str}_g{logg_hi_str}"]

    # Interpolate fluxes to lower model wavelength spectrum
    wvln    = spec00['wvln']

    spec0x  = spec00['flux'] + np.interp(wvln,spec01['wvln'],spec01['flux']) * diff_logg/0.5
    spec1x  = np.interp(wvln,spec10['wvln'],spec01['flux'])*(1-diff_logg/0.5) + np.interp(wvln,spec11['wvln'],spec11['flux'])*diff_logg/0.5

    spec    = spec0x*(1-diff_t/500) + spec1x*(diff_t/500)
    return wvln,spec

def model_flux_to_phot(mwvln,mflux):
    """
    
    """
    sdss_model_flux = np.array([])
    pans_model_flux = np.array([])
    skym_model_flux = np.array([])
    tmas_model_flux = np.array([])

    for f in ['u','g','r','i','z']:
        pass_fname  = f"filters/sdss/{f}.dat"
        model_flux  = convolveModelWithBandpass(mwvln,mflux,pass_fname)
        sdss_model_flux  = np.concat([sdss_model_flux,np.array([model_flux])])
        
    for f in ['g','r','i','z','y']:
        pass_fname  = f"filters/panstarrs/{f}.dat"
        model_flux  = convolveModelWithBandpass(mwvln,mflux,pass_fname)
        pans_model_flux = np.concat([pans_model_flux,np.array([model_flux])])
        
    for f in ['u','g','r','i','z']:
        pass_fname  = f"filters/skymapper/{f}.dat"
        model_flux  = convolveModelWithBandpass(mwvln,mflux,pass_fname)
        skym_model_flux = np.concat([skym_model_flux,np.array([model_flux])])

    for f in ['j','h','k']:
        pass_fname = f"filters/2mass/{f}.dat"
        model_flux = convolveModelWithBandpass(mwvln,mflux,pass_fname)
        tmas_model_flux = np.concat([tmas_model_flux,np.array([model_flux])])

    return np.concat([sdss_model_flux,pans_model_flux,skym_model_flux,tmas_model_flux])

##################################################################################################################

### MCMC

def log_likelihood(theta,x,y,y_err,spec_dict):
    """
    Log likelihood function
    """

    teff,logg,r,plx   = theta
    d = (1000 / plx) * pc2m
    rm = r*R_earth

    # Forward model
    wvln,flux       = interp_spectrum_dict(teff,logg,spec_dict)
    phot_model_emi  = model_flux_to_phot(wvln,flux)
    phot_model_obs  = phot_model_emi*rm*rm/d/d

    # Chi2
    score       = chi2(y,phot_model_obs,y_err)
    if np.nan_to_num(score,nan=-1)==-1:
        print(y)
        print(phot_model_obs)
        print(y_err)
        print("FUCKFUCKFUCKFUCKFUCKFUCK")

    return  -0.5*score

def log_prior(theta):           # contains constants dependent on the model grid bounds
    """
    Log prior function
    """

    teff,logg,r,plx   = theta

    # distance prior from gaia data

    # log prior -- containts hard-coded constants set by model grid
    if (teff < 10000)   or (teff > 25000)\
    or (logg < 6.0)    or (logg > 9.0)\
    or (r < 0.2)        or (r > 2):
        return -np.inf

    else:
        return -0.5 * ( (plx - P0) / sig_P0 )**2

def log_prob(theta,x,y,y_err,spec_dict):
    """
    Log probability
    """
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    else:
        return lp + log_likelihood(theta,x,y,y_err,spec_dict)

##################################################################################################################

### MISC

def shrinkScenarios():
    """
    Shrink scenario folder sizes by removing large files that are not used.
    """
    for scen in os.listdir(SCEN):
        if scen=='.DS_Store' or scen=='.gitignore': continue
        if 'fort.81' in os.listdir(f'{SCEN}/{scen}/lte'):
            os.remove(f'{SCEN}/{scen}/lte/fort.81')
        if 'fort.82' in os.listdir(f'{SCEN}/{scen}/lte'):
            os.remove(f'{SCEN}/{scen}/lte/fort.82')
        if 'fort.83' in os.listdir(f'{SCEN}/{scen}/lte'):
            os.remove(f'{SCEN}/{scen}/lte/fort.83')

##################################################################################################################

spec_dict = loadSpectra()

name    = 'Gaia DR3 3251748915515143296'

# Get photometry and parallax
wvln_obs,flux_obs,err_obs,plx,plx_err = photometry(name)
P0      = plx
sig_P0  = plx_err

# variables for emcee
ndim , nwalkers , nstep , nburn = 4 , 400 , 1000 , 200

# initial guess of parameters
x0  = np.array([15000,8,1,10])
pos = [x0*(1 + 1e-2*np.random.randn((4))) for w in range(nwalkers)]

# backup file
fn = 'simple emcee.h5'
backend = emcee.backends.HDFBackend( fn )
backend.reset( nwalkers , ndim )

# The big show
sampler = emcee.EnsembleSampler(nwalkers,ndim,log_prob,args=(wvln_obs,flux_obs,err_obs,spec_dict),backend=backend)
sampler.run_mcmc(pos,nstep,progress=True,skip_initial_state_check=True)

# Get results
reader  = emcee.backends.HDFBackend(fn)
nthin   = 2
samples = reader.get_chain(discard=nburn,thin=nthin,flat=True)

# Get best fitting model
lnprob  = reader.get_log_prob(discard=nburn,flat=True,thin=nthin)
lnpmax  = np.amax(lnprob)
xmax    = samples[np.where(lnprob==lnpmax)][0]
print(xmax)

# Plot results
fig     = corner.corner(samples,quantiles=[0.16,0.5,0.84],show_titles=True,labels=["T",r"$\log(g)$","R",r"$\pi$"])#,range=[[13000,25000],[7.75,8.25],[3e6,1.5e7],[1,100]])
plt.show()

input()


