from matplotlib import pyplot as plt
import numpy as np
import os,sys
import xml.etree.ElementTree as ET

if __name__ == '__main__':
    """
    Plot output spectrum from a converged SYNSPEC run
    """
    numargs     = len(sys.argv)
    scenario    = sys.argv[1]
    if '--lte' in sys.argv:
        path        = os.getcwd()+f'/tlusty/scenarios/{scenario}/synspec_lte/'
    else:
        path        = os.getcwd()+f'/tlusty/scenarios/{scenario}/synspec/'
    try:
        spec_file = np.loadtxt(path+'fort.7',dtype=str)
    except FileNotFoundError:
        print("[ERROR]      File \'fort.7\' does not exist. Exiting.")
        sys.exit(0)
    lam = np.array(spec_file[:,0],dtype=np.float32)
    flux = []
    for i in range(len(spec_file)):
        try:
            val = np.float32(spec_file[i,1])
        except ValueError:
            val = np.nan
        flux.append(val)
    
    compfile = None
    if numargs>2:
        compfile = sys.argv[2:]

    if compfile:
        fig,axes = plt.subplots(ncols=2)
        ax = axes[0]
        for f in compfile:
            if f=='--lte': continue

            ### Parse files from Koester grid (just some leading lines, plus flux to Eddington flux conversion)
            if 'koester' in f:
                dat = np.loadtxt(f,skiprows=6)
                grid_lam    = dat[:,0]
                grid_flux   = dat[:,1]/(4*np.pi)    # Flux --> Eddington Flux

            ### txt files from previous model runs, Bohlin+20 grid or the Hubeny grid Siyi uploaded
            else:
                dat = np.loadtxt(f,dtype=str)
                grid_lam = np.array(dat[:,0],dtype=np.float32)
                grid_flux = []
                for k in range(len(dat)):       # 1000nm - 3000nm
                    try:
                        val = np.float32(dat[k,1])
                    except ValueError:
                        val = np.nan
                    grid_flux.append(val)
                grid_flux = np.nan_to_num(np.array(grid_flux),nan=0)
                # grid_lam    = dat[:,0]
                # grid_flux   = dat[:,1]
        
            grid_lam = np.array(grid_lam,dtype=np.float64)
            grid_flux = np.array(grid_flux,dtype=np.float64)

            ax.plot(grid_lam,grid_flux,label=f)

            # Plot the ratio of the spectrum of interest, extrapolated to grid wavelengths,
            #   over the grid spectrum
            x = grid_lam
            y = np.interp(x,lam,flux)
            axes[1].plot(x,y/grid_flux,label=f)
            axes[1].legend()
        
        axes[1].set_xlabel("Wavelength [Å]",size='large')
        axes[1].set_xlim(900,10000)
        axes[1].set_ylim(0,2)
    else:
        fig,ax = plt.subplots()
    ax.plot(lam,flux,c='black',linewidth=1,label="Model",zorder=1000)
    ax.set_xlim(900,10000)
    ax.set_xlabel("Wavelength [Å]",size='large')
    ax.set_ylabel(r"Flux [erg/cm$^2$/s/Å]",size='large')
    ax.set_yscale('log')
    # ax.set_xlim(left=900,right=2000)
    # ax.set_ylim(1e5)
    ax.legend()

    fig.set_figwidth(10)
    plt.tight_layout()
    plt.savefig(path+'spec.png',dpi=200)
    plt.show()
