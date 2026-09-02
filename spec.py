from consts import *
import matplotlib.pyplot as plt
import numpy as np
import os,subprocess,time

### MCMC
import corner,emcee
from multiprocessing import Pool

### Other
from utils import *

##########################################################################################

def bilinear_interp(teff,logg,spec_dict):

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

    # Retrieve spectra
    spec00      = spec_dict[f"t{t_lo_str}_g{logg_lo_str}"]
    spec01      = spec_dict[f"t{t_lo_str}_g{logg_hi_str}"]
    spec10      = spec_dict[f"t{t_hi_str}_g{logg_lo_str}"]
    spec11      = spec_dict[f"t{t_hi_str}_g{logg_hi_str}"]
    wvln        = spec00['wvln']

    # Interpolate
    spec0x      = spec00['flux']*(1-diff_logg/0.5) + np.interp(wvln,spec01['wvln'],spec01['flux'])*diff_logg/0.5
    spec1x      = spec10['flux']*(1-diff_logg/0.5) + np.interp(wvln,spec11['wvln'],spec11['flux'])*diff_logg/0.5
    spec_final  = spec0x*(1-diff_t/500) + spec1x*diff_t/500

    return wvln,spec_final
    
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
    spec11x     = np.interp(wvln,spec110['wvln'],spec110['flux'])*(1-diff_abn/0.5) + np.interp(wvln,spec111['wvln'],spec110['flux'])*diff_abn/0.5

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
def fit_line(spec:str,scen:str,anum:int,minabn,maxabn,ddex,linecen,width,plotit:bool=False):

    # Open and trim observed flux
    wvln_all,flux_all,sigma_all = openCOSfile(spec)
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
        wvln_m_all,flux_m_all   = convolveModelWithCOS(wvln_m_all,flux_m_all,spec)

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

        # fig,ax = plt.subplots()
        # ax.plot(wvln,flux_m,zorder=100)
        # ax.plot(wvln,flux)
        # ax.set_title(abnstr)
        # plt.show()

    res = np.argmin(chi2_grid)
    min_chi2 = np.min(chi2_grid)
    min_abn_ind,min_vrad_ind = np.unravel_index(res,chi2_grid.shape)
    min_abn     = np.format_float_positional(abn_range[min_abn_ind],precision=1,min_digits=1)
    min_vrad    = vrad_range[min_vrad_ind]

    prob_grid   = np.exp(-0.5*chi2_grid)
    vrad_marg   = np.sum(prob_grid,axis=0)/np.sum(prob_grid)
    abn_marg    = np.sum(prob_grid,axis=1)/np.sum(prob_grid)

    if plotit:
        

        fig,ax = plt.subplots(ncols=2)
        ax[0].plot(abn_range,abn_marg)
        ax[0].set_title("Abundance")

        ax[1].plot(vrad_range,vrad_marg)
        ax[1].set_title("Radial velocity")

        plt.show()
        plt.close()

    # return min_abn,min_vrad,min_chi2
    return abn_range,abn_marg

    # # Retrieve best fit spectrum
    # fname = SCEN+f"/{scen}/synspec/{str(anum)}_{min_abn}.7"
    # wvln_m,flux_m = openModel(fname)
    # wvln_m,flux_m = convolveModelWithCOS(wvln_m,flux_m,spec)

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

# Forward model
def log_likelihood(theta,wvln,flux,sigma,width,spec_dict,fit_halpha:bool=False,plotit:bool=False):

    # unpack model parameters
    teff,logg,vrad = theta

    # use bilinear interpolation to get synthetic spectrum
    wvln_model,flux_model = bilinear_interp(teff,logg,spec_dict)

    # model spectrum uses vacuum wavelengths -- switch to air (Morton 1991) (not sure how much this matters)
    wvln_model = wvln_model / (1 + 2.735e-4 + 131.4182/(wvln_model**2) + 2.7625e8/(wvln_model**4) ) 

    # convolve model spectrum with a gaussian
    wvln_model,flux_model = convolveModelWithGaussian(wvln_model,flux_model,wvln,width)

    # Radial velocity shift
    beta = vrad/3e5
    wvln_model *= np.sqrt((1+beta)/(1-beta))
    flux_model = np.interp(wvln,wvln_model,flux_model)

    ### Fit first 5 Balmer lines
    line_centers    = [6562.8,4861.4,4340.5,4101.7,3970.1,3889.1][(not fit_halpha):]
    widths          = [150,120,80,50,30,20][(not fit_halpha):]
    labels          = [r"H$\alpha$",r"H$\beta$",r"H$\gamma$",r"H$\delta$",r"H$\epsilon$",r"H$\zeta$"][(not fit_halpha):]
    sum_chi2        = 0
    if plotit: fig,ax = plt.subplots()

    # For each of the first 4/5 Balmer lines:
    for i in range(len(line_centers)):
        # print(i)

        # normalize observed and modelled spectrua to average spectra at a fixed distance from the line center
        blue_ave        = np.average(flux[np.abs(wvln-(line_centers[i]-widths[i])) <= 2])
        red_ave         = np.average(flux[np.abs(wvln-(line_centers[i]+widths[i])) <= 2]) #
        cut_wvln        = wvln[np.abs(wvln-line_centers[i]) <= widths[i]]
        cut_flux_obs    = flux[np.abs(wvln-line_centers[i]) <= widths[i]]
        
        norm_obs        = np.interp(cut_wvln,
                                    [cut_wvln[0],cut_wvln[-1]],
                                    [blue_ave,red_ave])
                                    # [cut_flux_obs[0],cut_flux_obs[-1]])

        cut_flux_model  = flux_model[np.abs(wvln-line_centers[i]) <= widths[i]]
        norm_model      = np.interp(cut_wvln,
                                    [cut_wvln[0],cut_wvln[-1]],
                                    [cut_flux_model[0],cut_flux_model[-1]])

        # calculate chi^2 of the fit and add it to the total chi2
        # to get new sigmas for normalized flux, divide by original flux to get percentage differences
        cut_sigma       = sigma[np.abs(wvln-line_centers[i]) <= widths[i]]/cut_flux_obs

        score = chi2(cut_flux_obs/norm_obs,cut_flux_model/norm_model,cut_sigma)
        sum_chi2 += score

        if plotit:
            ax.plot(cut_wvln-line_centers[i],cut_flux_obs/norm_obs + 0.3*i,c='black')
            ax.plot(cut_wvln-line_centers[i],cut_flux_model/norm_model + 0.3*i,c='red',zorder=100)
            ax.text(-widths[i]-20,1+0.3*i,labels[i])

    if plotit:
        ax.set_xlabel(r"$\Delta\lambda\,\, [\AA]$",size='large')
        ax.set_ylabel("Normalized flux",size='large')
        ax.set_xlim(-175,175)
        plt.show()
        plt.close()

    ### Return log-likelihood using the sum of chi^2 values from each balmer line
    return -0.5*sum_chi2

def log_prior(theta):

    teff,logg,vrad = theta
    if (teff < 10000) or (teff > 25000)\
    or (logg < 6) or (logg >= 9)\
    or (vrad < -200) or (vrad > 200):
        return -np.inf
    else:
        return 0

def log_prob(theta,wvln,flux,sigma,width,spec_dict):

    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    else:
        return lp + log_likelihood(theta,wvln,flux,sigma,width,spec_dict)

################################################################################################

### Open observed spectrum

# X-Shooter
file_uvb = 'fits files/gd56_xsh_uvb.fits'
wvln,flux,sigma,res_uvb = openXshooter(file_uvb)
# file_vis = 'fits files/g29-38_uves_vis.fits'
# wvln2,flux2,sigma2,res_vis = openUVES(file_vis)

# wvln    = np.concat([wvln,wvln2])
# flux    = np.concat([flux,flux2])
# sigma   = np.concat([sigma,sigma2])

# fig,ax = plt.subplots()
# ax.plot(wvln,flux)
# plt.show()

# Convolve model spectrum with a Gaussian with width given by the lowest resolution
# min_res = min(res_uvb,res_vis)
min_res = res_uvb
width = 5000/min_res

# Load spectra dictionary
spec_dict = loadSpectra()

### TESTING

res = log_likelihood([15000,8.0,0],wvln,flux,sigma,width,spec_dict,fit_halpha=False,plotit=True)
input()
input()
### MCMC

# variables for emcee
ndim , nwalkers , nstep , nburn = 3 , 200 , 500 , 200

# initial guess of parameters
x0  = np.array([15000,8,10])
pos = [x0*(1 + 1e-3*np.random.randn((3))) for w in range(nwalkers)]

# backup file
fn = 'simple emcee.h5'
backend = emcee.backends.HDFBackend( fn )
backend.reset( nwalkers , ndim )

# The big show
sampler = emcee.EnsembleSampler(nwalkers,ndim,log_prob,args=(wvln,flux,sigma,width,spec_dict),backend=backend)
sampler.run_mcmc(pos,nstep,progress=True,skip_initial_state_check=True)

# Get results
reader  = emcee.backends.HDFBackend(fn)
nthin   = 2
samples = reader.get_chain(discard=nburn,thin=nthin,flat=True)

# Get best-fitting individual model values
lnprob  = reader.get_log_prob(discard=nburn,flat=True,thin=nthin)
lnpmax  = np.amax(lnprob)
xbest   = samples[np.where(lnprob==lnpmax)][0]
print(lnpmax,xbest)
# xmax    = [11960,8.23,1.3,57.02]

# Get mean values
xmean = []
for i in range(ndim):
    mean = np.percentile(samples[:,i],50)
    print(mean)
    xmean.append(mean)
lnpmean = log_likelihood(xmean,wvln,flux,sigma,width,spec_dict,plotit=False)
print(lnpmean,xmean)

### Plot results
# MCMC
fig     = corner.corner(samples,quantiles=[0.16,0.5,0.84],show_titles=True,labels=["T",r"$\log(g)$",r"$v_{\rm rad}$"])#,range=[[13000,25000],[7.75,8.25],[3e6,1.5e7],[1,100]])
plt.show()
plt.close()

### Plot fits to Balmer lines

# unpack model parameters
teff,logg,vrad = xbest

# use bilinear interpolation to get synthetic spectrum
wvln_model,flux_model = bilinear_interp(teff,logg,spec_dict)

# model spectrum uses vacuum wavelengths -- switch to air (Morton 1991) (not sure how much this matters)
wvln_model = wvln_model / (1 + 2.735e-4 + 131.4182/(wvln_model**2) + 2.7625e8/(wvln_model**4) ) 

# convolve model spectrum with a gaussian
wvln_model,flux_model = convolveModelWithGaussian(wvln_model,flux_model,wvln,width)

# Radial velocity shift
beta = vrad/3e5
wvln_model *= np.sqrt((1+beta)/(1-beta))
flux_model = np.interp(wvln,wvln_model,flux_model)

### Fit first 5 Balmer lines
fit_halpha = False
line_centers    = [6562.8,4861.4,4340.5,4101.7,3970.1,3889.1][(not fit_halpha):]
widths          = [150,120,80,50,30,20][(not fit_halpha):]
labels          = [r"H$\alpha$",r"H$\beta$",r"H$\gamma$",r"H$\delta$",r"H$\epsilon$",r"H$\zeta$"][(not fit_halpha):]
sum_chi2        = 0

fig,ax = plt.subplots()

# For each of the first 4/5 Balmer lines:
for i in range(len(line_centers)):
    # print(i)

    # normalize observed and modelled spectrua to average spectra at a fixed distance from the line center
    blue_ave        = np.average(flux[np.abs(wvln-(line_centers[i]-widths[i])) <= 2])
    red_ave         = np.average(flux[np.abs(wvln-(line_centers[i]+widths[i])) <= 2]) #
    cut_wvln        = wvln[np.abs(wvln-line_centers[i]) <= widths[i]]
    cut_flux_obs    = flux[np.abs(wvln-line_centers[i]) <= widths[i]]
    
    norm_obs        = np.interp(cut_wvln,
                                [cut_wvln[0],cut_wvln[-1]],
                                [blue_ave,red_ave])
                                # [cut_flux_obs[0],cut_flux_obs[-1]])

    cut_flux_model  = flux_model[np.abs(wvln-line_centers[i]) <= widths[i]]
    norm_model      = np.interp(cut_wvln,
                                [cut_wvln[0],cut_wvln[-1]],
                                [cut_flux_model[0],cut_flux_model[-1]])

    # calculate chi^2 of the fit and add it to the total chi2
    # to get new sigmas for normalized flux, divide by original flux to get percentage differences
    cut_sigma       = sigma[np.abs(wvln-line_centers[i]) <= widths[i]]/cut_flux_obs

    score = chi2(cut_flux_obs/norm_obs,cut_flux_model/norm_model,cut_sigma)
    sum_chi2 += score

    ax.plot(cut_wvln-line_centers[i],cut_flux_obs/norm_obs + 0.3*i,c='black')
    ax.plot(cut_wvln-line_centers[i],cut_flux_model/norm_model + 0.3*i,c='red',zorder=100)
    ax.text(-widths[i]-20,1+0.3*i,labels[i])

# unpack model parameters
teff,logg,vrad = xmean

# use bilinear interpolation to get synthetic spectrum
wvln_model,flux_model = bilinear_interp(teff,logg,spec_dict)

# model spectrum uses vacuum wavelengths -- switch to air (Morton 1991) (not sure how much this matters)
wvln_model = wvln_model / (1 + 2.735e-4 + 131.4182/(wvln_model**2) + 2.7625e8/(wvln_model**4) ) 

# convolve model spectrum with a gaussian
wvln_model,flux_model = convolveModelWithGaussian(wvln_model,flux_model,wvln,width)

# Radial velocity shift
beta = vrad/3e5
wvln_model *= np.sqrt((1+beta)/(1-beta))
flux_model = np.interp(wvln,wvln_model,flux_model)

### Fit first 5 Balmer lines
fit_halpha = False
line_centers    = [6562.8,4861.4,4340.5,4101.7,3970.1,3889.1][(not fit_halpha):]
widths          = [150,120,80,50,30,20][(not fit_halpha):]
labels          = [r"H$\alpha$",r"H$\beta$",r"H$\gamma$",r"H$\delta$",r"H$\epsilon$",r"H$\zeta$"][(not fit_halpha):]
sum_chi2        = 0

# fig,ax = plt.subplots()

# For each of the first 4/5 Balmer lines:
for i in range(len(line_centers)):
    # print(i)

    # normalize observed and modelled spectrua to average spectra at a fixed distance from the line center
    blue_ave        = np.average(flux[np.abs(wvln-(line_centers[i]-widths[i])) <= 2])
    red_ave         = np.average(flux[np.abs(wvln-(line_centers[i]+widths[i])) <= 2]) #
    cut_wvln        = wvln[np.abs(wvln-line_centers[i]) <= widths[i]]
    cut_flux_obs    = flux[np.abs(wvln-line_centers[i]) <= widths[i]]
    
    norm_obs        = np.interp(cut_wvln,
                                [cut_wvln[0],cut_wvln[-1]],
                                [blue_ave,red_ave])
                                # [cut_flux_obs[0],cut_flux_obs[-1]])

    cut_flux_model  = flux_model[np.abs(wvln-line_centers[i]) <= widths[i]]
    norm_model      = np.interp(cut_wvln,
                                [cut_wvln[0],cut_wvln[-1]],
                                [cut_flux_model[0],cut_flux_model[-1]])
    

    # calculate chi^2 of the fit and add it to the total chi2
        # to get new sigmas for normalized flux, divide by original flux to get percentage differences
    cut_sigma       = sigma[np.abs(wvln-line_centers[i]) <= widths[i]]/cut_flux_obs

    score = chi2(cut_flux_obs/norm_obs,cut_flux_model/norm_model,cut_sigma)
    sum_chi2 += score

    ax.plot(cut_wvln-line_centers[i],cut_flux_obs/norm_obs + 0.3*i,c='black')
    ax.plot(cut_wvln-line_centers[i],cut_flux_model/norm_model + 0.3*i,c='forestgreen',zorder=100)
    ax.text(-widths[i]-20,1+0.3*i,labels[i])

ax.plot([],[],c='black',label='Observed')
ax.plot([],[],c='red',label=r"Mean values, $\chi_\nu^2=$"+str(round(-2*lnpmean,2)))
ax.plot([],[],c='forestgreen',label=r"Best fit model, $\chi_\nu^2=$"+str(round(-2*lnpmax)))
ax.set_xlabel(r"$\Delta\lambda\,\, [\AA]$",size='large')
ax.set_ylabel("Normalized flux",size='large')
ax.set_xlim(-175,175)
ax.set_ylim(bottom=0)
ax.legend()
plt.show()
plt.close()

log_likelihood(xbest,wvln,flux,sigma,width,spec_dict,plotit=True)


