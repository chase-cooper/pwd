from matplotlib import pyplot as plt
import numpy as np
import os,sys

if __name__ == '__main__':
    scenario = sys.argv[1]
    step = sys.argv[2]

    ### Establish variables for relevant output files. Exit if any are missing
    unit6,unit7,unit9,unit69    = None,None,None,None
    for f in os.listdir(os.getcwd()+f"/tlusty/scenarios/{scenario}/{step}"):
        if f[-2:] == '.6':  unit6 = os.getcwd()+f"/tlusty/scenarios/{scenario}/{step}/{f}"
        if f[-2:] == '.7':  unit7 = os.getcwd()+f"/tlusty/scenarios/{scenario}/{step}/{f}"
        if f[-2:] == '.9':  unit9 = os.getcwd()+f"/tlusty/scenarios/{scenario}/{step}/{f}"
        if f[-3:] == '.69': unit69= os.getcwd()+f"/tlusty/scenarios/{scenario}/{step}/{f}"

    print(unit6,unit7,unit9,unit69)

    ### Begin plots, gather quantities
    fig,ax = plt.subplots(nrows=2,ncols=3)

    ### Unit 7: # depth points, # parameters, layer depths, layer values
    f7_raw  = open(unit7,'r')
    f7      = f7_raw.read().split()
    f7_raw.close()

    ndepth,nparams = int(f7[0]),int(f7[1])
    depths  = np.zeros((ndepth))
    temps   = np.zeros((ndepth))
    for i in range(ndepth):
        val = f7[2+i]
        val = val.replace('D','E')
        depths[i] = np.float32(val)

        tval    = f7[2+ndepth+nparams*(i)+1]
        tval    = tval.replace('D','E')
        temps[i] = np.float32(tval)

    ### Unit 9: convergence criteria
    f9 = np.loadtxt(unit9,dtype=str,skiprows=3)
    iteration   = np.int8(f9[:,0])
    layer       = np.int16(f9[:,1])
    delta_t     = np.float32([float(val.replace('D','E')) for val in f9[:,2]])
    delta_ne    = np.float32([float(val.replace('D','E')) for val in f9[:,3]])
    delta_pop   = np.float32([float(val.replace('D','E')) for val in f9[:,4]])
    delta_rad   = np.float32([float(val.replace('D','E')) for val in f9[:,5]])
    delta_max   = np.float32([float(val.replace('D','E')) for val in f9[:,6]])
    lvl_id      = np.int8(f9[:,7])
    freq_id     = np.int8(f9[:,8])


    ### Plot relative change in temperature with log depth
    for i in range(max(iteration)):
        x = depths
        y = delta_t[ndepth*(i):ndepth*(i+1)][::-1]
        ax[0,0].plot(x,y,label=f"Iter{i+1}")
    ax[0,0].set_xscale('log')
    ax[0,0].set_xlabel("Log depth (mass)")
    ax[0,0].set_ylabel(r"$\Delta$T/T")
    # ax[0,0].legend(draggable=True)

    ### Plot relative change in temperature with log depth IN LOG SPACE
    for i in range(max(iteration)):
        x = depths
        y = delta_t[ndepth*(i):ndepth*(i+1)][::-1]
        ax[0,1].plot(x,abs(y))
    ax[0,1].set_xscale('log')
    ax[0,1].set_yscale('log')
    ax[0,1].set_xlabel("Log depth (mass)")
    ax[0,1].set_ylabel(r"$\Delta$T/T")

    ### Plot max temperature changes per iteration
    x = range(1,np.max(iteration)+1)
    y = []
    for i in range(max(iteration)):
        y.append(np.max(np.abs(delta_t[70*i:70*(i+1)])))
    ax[0,2].plot(x,y,marker='s',c='k')
    ax[0,2].set_yscale('log')
    ax[0,2].set_xlabel("Iteration")
    ax[0,2].set_ylabel(r"Max $\Delta$T/T")
        
    ### Plot max relative change in state vector
    for i in range(max(iteration)):
        x = depths
        y = delta_max[ndepth*(i):ndepth*(i+1)][::-1]
        ax[1,0].plot(x,y)
    ax[1,0].set_xscale('log')
    ax[1,0].set_xlabel("Log depth (mass)")
    ax[1,0].set_ylabel(r"Max $\Delta\Phi$/$\Phi$")
    
    ### Plot max relative change in state vector IN LOG SPACE
    for i in range(max(iteration)):
        x = depths
        y = delta_max[ndepth*(i):ndepth*(i+1)][::-1]
        ax[1,1].plot(x,abs(y))
    ax[1,1].set_xscale('log')
    ax[1,1].set_yscale('log')
    ax[1,1].set_xlabel("Log depth (mass)")
    ax[1,1].set_ylabel(r"Max $\Delta\Phi$/$\Phi$")

    ### Plot max change in state vector per iteration
    x = range(1,np.max(iteration)+1)
    y = []
    for i in range(max(iteration)):
        y.append(np.max(np.abs(delta_max[70*i:70*(i+1)])))
    ax[1,2].plot(x,y,marker='s',c='k')
    ax[1,2].set_yscale('log')
    ax[1,2].set_xlabel("Iteration")
    ax[1,2].set_ylabel(r"Max $\Delta\Phi$/$\Phi$")

    # plt.tight_layout()
    plt.show()
    

