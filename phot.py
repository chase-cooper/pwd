from astropy import coordinates as coords
from astroquery.gaia import Gaia
from astroquery.mast import Catalogs
from astroquery.mast import Mast
from astroquery.sdss import SDSS
from astroquery.simbad import Simbad
from dustmaps.sfd import SFDQuery
from dustmaps.config import config
# config['data_dir'] = "/Users/cooper/Documents/Research/polluted_white_dwarfs/dustmaps"
# import dustmaps.sfd
# dustmaps.sfd.fetch()
# input()
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import LinearNDInterpolator
import os,requests,shutil,subprocess,sys,time
from consts import SCEN
from utils import newSpeciesDict

import corner,emcee
from multiprocessing import Pool

# initializing MAST queries
mast = Mast()

# Extinction
sfd     = SFDQuery()

# Constants
c       = 3e8
pc2m    = 3.086e16
pc2cm   = 3.086e18
grav    = 6.67e-11  # [m^3 kg^-1 s^-2]
msol    = 1.989e30  # [kg]
R_earth = 6.371e6   # [m]

##################################################################################################################

def chi2(y,ymodel,yerr):
    """
    Calculate the reduced chi-squared score of a model fit, assuming the model 
    is parameterized by two values ( T_eff, log(g) )
    """
    res = (y - ymodel)**2 / yerr**2
    # print(res)
    return np.sum(res)/(len(y)-2)

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

    return plx,plx_err,ext_frac,id_sdss,id_pans,id_skym

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

### Magnitude-to-flux conversion

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

### Model spectrum convolution

def convolveModelWithPassband(model_fname,band_fname):
    """
    Convolve a TLUSTY/SYNSPEC spectrum (at location 'model_fname') with a filter response function
    (at location 'band_fname').
    """

    # Speed of light in AA/s
    cAA = 3e18

    model       = np.genfromtxt(model_fname)
    model_wvlns = model[1:,0]
    model_flux  = model[1:,1]
    model_dlam  = model[1:,0] - model[:-1,0]

    # TLUSTY/SYNSPEC spectra are Eddington fluxes, so require a factor of 4*pi to get fluxes
    model_flux  *= 4*np.pi

    filt    = np.genfromtxt(band_fname)
    filt_wvlns  = filt[:,0]
    filt_resp   = filt[:,1]

    # Interpolate response to model wavelengths
    model_response  = np.interp(model_wvlns,filt_wvlns,filt_resp,left=0,right=0)

    # Resulting flux in [erg/s/cm^2/Hz]
    conv_flux   = (1/cAA) * np.sum(model_response*model_flux*model_wvlns*model_dlam) / np.sum(model_response * model_dlam / model_wvlns)
    conv_flux   *= 1e23     # Convert to Jy
    return conv_flux

def convolve2(mwvln,mflux,band_fname):
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
    # clear

    # Interpolate response to model wavelengths
    model_response  = np.interp(model_wvln0,filt_wvlns,filt_resp,left=0,right=0)
    # clear

    # Resulting flux in [erg/s/cm^2/Hz]
    conv_flux   = (1/cAA) * np.sum(model_response*model_flux0*model_wvln0*model_dlam0) / np.sum(model_response * model_dlam0 / model_wvln0)
    conv_flux   *= 1e23     # Convert to Jy
    return conv_flux

##################################################################################################################

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
    plx,plx_err,ext_frac,name_sdss,name_pans,name_skymapper = query_gaia(name_gaia,verbose=True)
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

    ##################################################################################################################

    ### Formatting of data and conversion of magnitudes to fluxes in Jy
    ###
    ### This kinda sucks but whatever it works

    # Filter effective wavelengths [nm] from SVO
    sdss_eff_wvlns  = np.array([360.8,467.18,614.11,745.79,892.28]) * 1e-9
    pans_eff_wvlns  = np.array([481.02,615.55,750.3,866.84,961.36]) * 1e-9
    skym_eff_wvlns  = np.array([350.02,501.6,607.69,773.28,912.03]) * 1e-9

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

    all_wvlns       = np.concat([sdss_eff_wvlns,pans_eff_wvlns,skym_eff_wvlns])
    all_fluxes      = np.concat([sdss_fluxes,pans_fluxes,skym_fluxes])
    all_fluxerrs_hi = np.concat([sdss_fluxerr_hi,pans_fluxerr_hi,skym_fluxerr_hi])

    return all_wvlns,all_fluxes,all_fluxerrs_hi,plx,plx_err

    ##################################################################################################################

    #########################
    ###   Model fitting   ###
    #########################

    def iterative_fit(sdss_flux =np.array([]),
                      pans_flux =np.array([]),
                      skym_flux =np.array([]),
                      sdss_errs =np.array([]),
                      pans_errs =np.array([]),
                      skym_errs =np.array([]),
                      num_iter  = 5):

        ### Compile fluxes and errors
        all_fluxes  = np.concat([sdss_flux,pans_flux,skym_flux])
        all_fluxerrs= np.concat([sdss_errs,pans_errs,skym_errs])

        ### Establish a function that, using derived temperature and log(g), interpolates over the
        ###     Bédard (2010) evolutionary models to get white dwarf mass. With mass and radius (from
        ###     the scale factor), one can get a new value of log(g).
        x,y,z = [],[],[]
        for file in os.listdir('mass_radius/AllSequences/'):
            model_mass    = int(file[4:7])/100
            f       = np.loadtxt(f"mass_radius/AllSequences/{file}",skiprows=5,usecols=(0,1,2))
            t       = f[::3,1]
            logg    = f[::3,2]
            for i in range(len(t)):
                x.append(t[i])
                y.append(logg[i])
                z.append(model_mass)
        X, Y = np.meshgrid(x,y)  # 2D grid for interpolation

        # Mass interpolation is called by supplying temperature and log(g) in that order
        mass_interp = LinearNDInterpolator(list(zip(x, y)), z)

        # Start with an assumed log(g) of 8
        log_g = 8

        ### Iteratively calculate the objects temperature and radius, and use that to get a new,
        ###     hopefully more accurate, value of log(g).

        for niter in range(num_iter):

            # Variables for finding chi2 minimum
            min_t               = 9999
            min_chi2            = 1e100
            min_scale           = 1e100
            best_model_fluxes   = np.array([])

            ### Iterate over models with different temperatures
            num_temp    = 51
            num_scale   = 500
            chi2_grid   = np.zeros((num_temp,num_scale))

            temps       = np.linspace(20000,25000,num_temp)

            for i in range(len(temps)):
                t   = temps[i]
                tstr = str(int(t))
                gstr = str(int(log_g*100))

                # If TLUSTY model doesn't already exist, run it. Model name follows Ivan Hubeny's 
                #   naming convention. This may be obsolete once we have a grid.
                model_name  = f"t{tstr}g{gstr}"
                model_exists= True
                num_tries   = 6
                if model_name not in os.listdir('tlusty/scenarios/'):
                    for j in range(num_tries):  # Try running the model up to 3 times

                        # If not on the first iteration, something must have gone wrong, so remove the
                        #   folder 'model_name'
                        if model_name in os.listdir('tlusty/scenarios/'):
                            shutil.rmtree(f"tlusty/scenarios/{model_name}")

                        # Run the model
                        subprocess.run(['python','tlusty.py',model_name,tstr,str(log_g)])

                        # Check model outputs for errors
                        try:
                            dat     = np.genfromtxt(f"tlusty/scenarios/{model_name}/lte/fort.14")
                            flux    = np.nan_to_num(dat[:,1],nan=-1)
                            if -1 not in flux:
                                model_exists = True
                                break
                            else:
                                model_exists = False
                                print(f"{model_name}: NaN's in spectrum. {num_tries-j-1} more attempts")
                                continue
                        except FileNotFoundError:
                            model_exists = False
                            print(f"{model_name}: fort.14 not found. {num_tries-j-1} more attempts")
                            continue
                    
                    # If after all attempts the model is still invalid, exit.
                    if model_exists == False:
                        print(f"{model_name}: Model failed to produce a valid spectrum. Exiting")
                        shutil.rmtree(f"tlusty/scenarios/{model_name}")
                        exit(0)

                # The file path of the resulting model flux in [erg/s/cm^2/AA] is something that needs
                #   to be hard-coded in, unfortunately. 
                # For now, taking the emergent spectra from TLUSTY
                spec_fname  = f"tlusty/scenarios/{model_name}/lte/fort.14"

                # Model convolved with passbands
                sdss_model_flux = np.array([])
                pans_model_flux = np.array([])
                skym_model_flux = np.array([])

            # if np.any(sdss_mags):
                for f in ['u','g','r','i','z']:
                    pass_fname  = f"filters/sdss/{f}.dat"
                    model_flux  = convolveModelWithPassband(spec_fname,pass_fname)
                    sdss_model_flux  = np.concat([sdss_model_flux,np.array([model_flux])])
                    
                for f in ['g','r','i','z','y']:
                    pass_fname  = f"filters/panstarrs/{f}.dat"
                    model_flux  = convolveModelWithPassband(spec_fname,pass_fname)
                    pans_model_flux = np.concat([pans_model_flux,np.array([model_flux])])
                    
                for f in ['u','g','r','i','z']:
                    pass_fname  = f"filters/skymapper/{f}.dat"
                    model_flux  = convolveModelWithPassband(spec_fname,pass_fname)
                    skym_model_flux = np.concat([skym_model_flux,np.array([model_flux])])


                
                all_model_fluxes    = np.concat([sdss_model_flux,pans_model_flux,skym_model_flux])
                print(all_model_fluxes)
                input()

                ### Iterate over scale factors
                approx_scale    = all_fluxes[0] / all_model_fluxes[0]   # <-- div by 0?
                approx_log      = np.log10(approx_scale)
                scales          = np.logspace(approx_log-2,approx_log+2,500)

                for k in range(len(scales)):
                    scale = scales[k]
                    scaled_model_fluxes = all_model_fluxes * scale
                    print(all_model_fluxes)
                    print(scale)
                    input()
                    score   = chi2(all_fluxes,scaled_model_fluxes,all_fluxerrs)
                    chi2_grid[i][k] = score
                    if score < min_chi2:
                        min_chi2            = score
                        min_t               = t
                        min_scale           = scale
                        best_model_fluxes   = scaled_model_fluxes
                    # print(t,scale,score)

            ### With best-fit scale factor, get best-fit radius. Use radius and temperature to get
            ###     a value of mass using evolutionary tables from Bédard+2010, then use that to get
            ###     new log(g).
            print('*'*100)
            print("Iteration "+str(niter+1))

            print(name+" Best Fit\n"+'-'*60)
            print("Temperature:     "+str(min_t))
            print("Scale factor:    "+str(min_scale))
            print(r"$\chi^2$:       ",end='')
            print(str(min_chi2))
            
            dist_m  = dist*pc2m
            print("Distance in meters:  "+str(dist_m))
            radius  = np.sqrt(min_scale*dist_m*dist_m)
            print("Radius in meters:    "+str(radius))

            mass_sol= mass_interp(min_t,log_g)  # Mass in solar masses
            mass_kg = mass_sol * msol           # kg
            print("Mass, solar masses:  "+str(mass_sol))

            new_g   = 100*grav*mass_kg/radius/radius    # surface gravity [cm/s^2]
            new_logg= np.log10(new_g)
            print("New log(g):          "+str(new_logg))

            print('*'*100)

            ### Use chi2 grid to get errorbars on temperature and scale factor
            likelihood  = np.exp(-0.5*chi2_grid)

            # Temperature
            dscale      = np.zeros_like(likelihood[0])
            dscale[1:]  = scales[1:] - scales[:-1]
            dscale[0]   = dscale[1]
            t_marg      = np.sum(likelihood*dscale,axis=1)
            t_pdf       = t_marg / np.sum(t_marg)
            t_cdf       = np.zeros_like(t_pdf)
            for i in range(len(t_cdf)):
                t_cdf[i] = np.sum(t_pdf[:i+1])

            # Scale
            sc_marg = np.sum(likelihood,axis=0)
            sc_pdf  = sc_marg / np.sum(sc_marg)
            sc_cdf  = np.zeros_like(sc_pdf)
            for i in range(len(sc_cdf)):
                sc_cdf[i] = np.sum(sc_pdf[:i+1])

            fig,ax = plt.subplots()
            ax.plot(t_pdf)
            ax.plot(t_cdf)
            plt.show()
            plt.close()

            ### Set log(g) to new value
            log_g   = new_logg

        return best_model_fluxes

    def lev_marq(sdss_flux  = np.array([]), # incomplete
                 pans_flux  = np.array([]),
                 skym_flux  = np.array([]),
                 sdss_errs  = np.array([]),
                 pans_errs  = np.array([]),
                 skym_errs  = np.array([]),
                 logg       = 8):

        ### Compile fluxes and errors
        all_fluxes  = np.concat([sdss_flux,pans_flux,skym_flux])
        all_fluxerrs= np.concat([sdss_errs,pans_errs,skym_errs])

    iterative_fit(sdss_fluxes,pans_fluxes,skym_fluxes,
                  sdss_fluxerr_hi,pans_fluxerr_hi,skym_fluxerr_hi)

################################################################################################

### MCMC

def loadSpectra():
    # Spectrum dictionary
    start = time.time()
    print("Loading spectra...",end='')

    global spec_dict
    spec_dict = {}
    for f in os.listdir(SCEN):
        if f in ['.DS_Store','.gitignore','metals','old grid','supermetals']: continue
        spec_dict[f] = {}
        dat = np.genfromtxt(SCEN+f'/{f}/synspec/fort.7')
        fwvln = dat[:,0]
        fflux = dat[:,1]
        spec_dict[f] = {
            'wvln':fwvln,
            'flux':fflux
        }

    end = time.time()
    print(f"Loaded spectra in {round(end-start,2)}s")

# Forward model functions, modified to assume log(g) = 8

def interp_spectrum(t,logg=8):
    """
    Given a temperature (anything) and log(g) (to 2 decimal places), linearly interpolate between
    adjacent spectra to produce a new spectrum
    """
    t_grid_res  = 500
    g_grid_res  = 0.1
    tstr_lo     = str(int(t - (t % t_grid_res) + 0.001))
    tstr_hi     = str(int(t+t_grid_res - (t%t_grid_res) + 0.001))
    gstr_lo     = str(int(100 * (logg - (logg % g_grid_res) + 0.001)))
    gstr_hi     = str(int(100 * (logg + g_grid_res - (logg*g_grid_res) + 0.001)))
    print(logg,gstr_lo,gstr_hi)

    ### Interpolate
    
    # Open lower model spectrum
    fname   = f'{SCEN}/t{tstr_lo}g{gstr_lo}/lte/fort.14'
    dat_lo  = np.genfromtxt(fname)
    wvln    = dat_lo[:,0]
    flux_lo = dat_lo[:,1]

    # Open upper model spectrum
    fname   = f'{SCEN}/t{tstr_hi}g{gstr_hi}/lte/fort.14'
    dat_hi  = np.genfromtxt(fname)
    wvln_hi = dat_hi[:,0]
    flux_hi = dat_hi[:,1]

    # Interpolate upper model flux to lower model wavelength spectrum
    flux_hi = np.interp(wvln,wvln_hi,flux_hi)

    # Interpolate between fluxes
    flux_diff   = flux_hi - flux_lo
    diff_t      = (t - int(tstr_lo))/t_grid_res
    diff_g      = (logg - int(gstr_lo)/100)/g_grid_res
    diff        = np.sqrt(diff_t**2 + diff_g**2)
    flux        = flux_lo + diff*flux_diff

    # write to temp file
    file    = open('temp.dat','w')
    for i in range(len(flux)):
        file.write(str(wvln[i])+'\t'+str(flux[i])+'\n')
    file.close()
    return wvln,flux

def interp_spectrum_dict(t,logg=8,spec_dict=newSpeciesDict()):
    """
    Given a temperature (anything) and log(g) (to 2 decimal places), linearly interpolate between
    adjacent spectra to produce a new spectrum. Instead of reading in data from files, the data
    are accessed through a dictionary where all the data has been preloaded.
    """
    t_grid_res  = 500
    g_grid_res  = 0.1
    tstr_lo     = str(int(t - (t % t_grid_res) + 0.001))
    tstr_hi     = str(int(t+t_grid_res - (t%t_grid_res) + 0.001))
    gstr_lo     = str(int(100 * (logg - (logg % g_grid_res) + 0.001)))
    gstr_hi     = str(int(100 * (logg + g_grid_res - (logg%g_grid_res) + 0.001)))
    # print(logg,gstr_lo,gstr_hi)

    # Open lower model spectrum
    mname   = f"t{tstr_lo}g{gstr_lo}"
    wvln    = spec_dict[mname]['wvln']
    flux_lo = spec_dict[mname]['flux']

    # Open upper model spectrum
    mname   = f"t{tstr_hi}g{gstr_hi}"
    wvln_hi = spec_dict[mname]['wvln']
    flux_hi = spec_dict[mname]['flux']

    # Interpolate upper model flux to lower model wavelength spectrum
    flux_hi = np.interp(wvln,wvln_hi,flux_hi)

    # Interpolate between fluxes
    flux_diff   = flux_hi - flux_lo
    diff_t      = (t - int(tstr_lo))/t_grid_res
    diff_g      = (logg - int(gstr_lo)/100)/g_grid_res
    diff        = np.sqrt(diff_t**2 + diff_g**2)
    flux        = flux_lo + diff*flux_diff

    # write to temp file
    # file    = open('temp.dat','w')
    # for i in range(len(flux)):
    #     file.write(str(wvln[i])+'\t'+str(flux[i])+'\n')
    # file.close()
    return wvln,flux

def model_flux_to_phot(mwvln,mflux):
    """
    
    """
    sdss_model_flux = np.array([])
    pans_model_flux = np.array([])
    skym_model_flux = np.array([])

    for f in ['u','g','r','i','z']:
        pass_fname  = f"filters/sdss/{f}.dat"
        model_flux  = convolve2(mwvln,mflux,pass_fname)
        sdss_model_flux  = np.concat([sdss_model_flux,np.array([model_flux])])
        
    for f in ['g','r','i','z','y']:
        pass_fname  = f"filters/panstarrs/{f}.dat"
        model_flux  = convolve2(mwvln,mflux,pass_fname)
        pans_model_flux = np.concat([pans_model_flux,np.array([model_flux])])
        
    for f in ['u','g','r','i','z']:
        pass_fname  = f"filters/skymapper/{f}.dat"
        model_flux  = convolve2(mwvln,mflux,pass_fname)
        skym_model_flux = np.concat([skym_model_flux,np.array([model_flux])])

    return np.concat([sdss_model_flux,pans_model_flux,skym_model_flux])

def model_flux_to_phot_fname(fname):
    """
    
    """
    sdss_model_flux = np.array([])
    pans_model_flux = np.array([])
    skym_model_flux = np.array([])

    for f in ['u','g','r','i','z']:
        pass_fname  = f"filters/sdss/{f}.dat"
        model_flux  = convolveModelWithPassband(fname,pass_fname)
        sdss_model_flux  = np.concat([sdss_model_flux,np.array([model_flux])])
        
    for f in ['g','r','i','z','y']:
        pass_fname  = f"filters/panstarrs/{f}.dat"
        model_flux  = convolveModelWithPassband(fname,pass_fname)
        pans_model_flux = np.concat([pans_model_flux,np.array([model_flux])])
        
    for f in ['u','g','r','i','z']:
        pass_fname  = f"filters/skymapper/{f}.dat"
        model_flux  = convolveModelWithPassband(fname,pass_fname)
        skym_model_flux = np.concat([skym_model_flux,np.array([model_flux])])

    return np.concat([sdss_model_flux,pans_model_flux,skym_model_flux])

def log_likelihood(theta,x,y,y_err):
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

# contains constants dependent on the model grid bounds
def log_prior(theta):
    """
    Log prior function
    """

    teff,logg,r,plx   = theta

    # distance prior from gaia data

    # log prior -- containts hard-coded constants set by model grid
    if (teff < 13000)   or (teff > 24500)\
    or (logg < 7.5)    or (logg > 8.5)\
    or (r < 0.2)        or (r > 2):
        return -np.inf

    else:
        return -0.5 * ( (plx - P0) / sig_P0 )**2

def log_prob(theta,x,y,y_err):
    """
    Log probability
    """
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    else:
        return lp + log_likelihood(theta,x,y,y_err)

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

# fname = SCEN+'/t16000g805/lte/fort.14'

# dat     = np.genfromtxt(fname)
# wvln    = dat[:,0].tolist()
# flux    = dat[:,1].tolist()

# print(wvln)
# print(flux)
# res     = model_flux_to_phot(wvln,flux)
# print(res)

# # res = model_flux_to_phot_fname(fname)
# # print(res)
# input()


loadSpectra()

name    = 'Gaia DR3 3251748915515143296'

# Get photometry and parallax
wvln_obs,flux_obs,err_obs,plx,plx_err = photometry(name)
P0      = plx
sig_P0  = plx_err

# variables for emcee
ndim , nwalkers , nstep , nburn = 4 , 200 , 500 , 200

# initial guess of parameters
x0  = np.array([15000,8,1,10])
pos = [x0*(1 + 1e-2*np.random.randn((4))) for w in range(nwalkers)]

# backup file
fn = 'simple emcee.h5'
backend = emcee.backends.HDFBackend( fn )
backend.reset( nwalkers , ndim )

# The big show
sampler = emcee.EnsembleSampler(nwalkers,ndim,log_prob,args=(wvln_obs,flux_obs,err_obs),backend=backend)
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


