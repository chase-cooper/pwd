# NOT LSF PURETLY SCRATCH FILE FOR FUCKING AROUND
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import os,subprocess
from scipy.interpolate import LinearNDInterpolator
from consts import *


def doubleModelComp():

    fig,ax = plt.subplots(nrows=3)

    dat1    = np.loadtxt('tlusty/scenarios/test/synspec/fort.7')
    # dat1    = np.loadtxt('tlusty/scenarios/t20000g800/synspec/fort.7')
    dat2    = np.loadtxt('tlusty/grid/hubeny/t20000g800n.spec')
    dat3    = np.loadtxt('tlusty/grid/hubeny2/t200g800n.spec')
    # dat3    = np.loadtxt('tlusty/grid/koester/da20000_800.dk.dat.txt')
    # dat3[:,1] = dat3[:,1]/(4*np.pi)

    ax[2].plot(dat1[:,0],dat1[:,1],label='My run')
    ax[2].plot(dat2[:,0],dat2[:,1],label='Cool metal-rich grid')
    ax[2].plot(dat3[:,0],dat3[:,1],label='Hot pure-H grid')
    ax[2].set_xlim(900,5000)
    ax[2].set_yscale('log')
    ax[2].legend()

    flux1 = np.interp(dat2[:,0],dat1[:,0],dat1[:,1])
    flux2 = np.interp(dat3[:,0],dat1[:,0],dat1[:,1])

    # Axis 1: Lyman series
    ax[0].plot(dat2[:,0],np.abs(1-flux1/dat2[:,1]),label='Model / metal grid')
    ax[0].plot(dat3[:,0],np.abs(1-flux2/dat3[:,1]),label='Model / pure-H grid')

    ax[0].set_yscale('log')
    ax[0].set_xlim(900,1500)
    ax[0].set_ylim(top=1)
    ax[0].hlines(0.01,xmin=900,xmax=1500,colors='black',linestyles='dotted',zorder=-100)
    ax[0].text(1450,0.015,'1%')
    ax[0].set_xlabel(r"Wavelength $\AA$")
    ax[0].set_ylabel("|Model / grid|")
    ax[0].legend()

    # # Axis 2: Balmer series
    ax[1].plot(dat2[:,0],np.abs(1-flux1/dat2[:,1]),label='Model / metal grid')
    ax[1].plot(dat3[:,0],np.abs(1-flux2/dat3[:,1]),label='Model / pure-H grid')

    ax[1].set_yscale('log')
    ax[1].set_xlim(3500,5000)
    ax[1].set_ylim(top=1)
    ax[1].hlines(0.01,xmin=3500,xmax=6000,colors='black',linestyles='dotted',zorder=-100)
    ax[1].text(4800,0.015,'1%')
    ax[1].set_xlabel(r"Wavelength $\AA$")
    ax[1].set_ylabel("|Model / grid|")
    ax[1].legend()

    fig.set_figheight(7)
    plt.tight_layout()
    # plt.savefig('figs/gridcomp_pureh_metal_t20000g800')
    plt.show()

def other():

    # dat1    = np.loadtxt('tlusty/metals/h3.14')
    dat1    = np.loadtxt('tlusty/scenarios/coolda/metals/fort.14')
    # dat1    = np.loadtxt('tlusty/scenarios/t20000g800_test/lte/fort.14')
    # dat2    = np.loadtxt('fort.7')
    dat3    = np.loadtxt('tlusty/metals 2/h2.14')
    # dat3    = np.loadtxt('coolda/t20000g800l.14')

    fig,ax = plt.subplots(nrows=2)

    ax[0].plot(dat1[:,0],dat1[:,1],label='My run')
    # ax[0].plot(dat2[:,0],dat2[:,1],label='Ivan\'s run')
    ax[0].plot(dat3[:,0],dat3[:,1],label='Ivan\'s Run')
    # ax[0].plot(dat3[:,0],dat3[:,1],label='Hubeny grid')
    ax[0].set_xscale('log')
    ax[0].set_yscale('log')
    ax[0].set_xlim(left=900,right=10000)

    flux51 = np.interp(dat3[:,0],dat1[:,0],dat1[:,1])
    # flux54 = np.interp(dat3[:,0],dat2[:,0],dat2[:,1])
    ax[1].plot(dat3[:,0],flux51/dat3[:,1],label="Me / Ivan")
    # ax[1].plot(dat3[:,0],flux54/dat3[:,1],label="SYNSPEC54 / Ivan")
    # ax[1].plot(dat1[:,0],dat1[:,1]/dat3[:,1],label='TLUSTY301 / TLUSTY205')
    ax[1].set_xscale('log')
    # ax[1].set_yscale('log')
    ax[1].set_ylim(0.8,1.2)
    ax[1].set_xlim(left=900,right=10000)

    ax[0].legend()
    ax[1].legend()
    plt.show()
    plt.close()

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

fig,ax = plt.subplots()

if True:
    tstr = '15000'
    for g in np.arange(6,9.1,0.5):
        gstr = np.format_float_positional(g,precision=1,min_digits=1)

        scen = f"t{tstr}_g{gstr}"

        dat = np.genfromtxt(f'tlusty/hydrogengrid/{scen}.spec')
        ax.plot(dat[:,0],dat[:,1],label=gstr)

ax.legend()
ax.set_yscale('log')
plt.show()