import numpy as np
import sys
import time
import shutil
import phydrus as ps
import pandas as pd
import contextlib


sys.path.append("../Phydrus_MRS_modules-py")
#path to the working directory where python modules have to be
from MRS_load_interp_mod import makeZarray_sim_extra_layers_Phydrus, read_interp_ker_csv, estimate_MRS_sig_from_sim_v4

def PHYDRUS_MRS_wrapper_numexp(X, Parameters):
    '''
    Function to forward simulate MRS signals for hydrologic params in X and fixed parameters in the dictionnary Parameters.
    To be used in the UQpyLab framework

    Parameters
    ----------
    X : list (Size Nsample x Nparams) [theta_r, theta_s, alpha, n, log10(Ks)
                                       ...
                                       same for each layer of the sim zone
                                       ...
                                       theta_u_percent_base, S_lim, theta_u_percent_max
                                       ...
                                       same for each layer of the sim zone
                                       ...
                                       theta_MRS_bedrock,
                                       dhpiezo]

    Parameters dict 

    Returns
    -------
    Y_E1 : list of numpy arrays (size : nsamples x n_print_times)
        List of estimated MRS signals for each parameter set contained in X for q1
    Y_E2 : list of numpy arrays (size : nsamples x n_print_times)
        List of estimated MRS signals for each parameter set contained in X for q2
    Y_E3 : list of numpy arrays (size : nsamples x n_print_times)
        List of estimated MRS signals for each parameter set contained in X for q3
    Y_E4 : list of numpy arrays (size : nsamples x n_print_times)
        List of estimated MRS signals for each parameter set contained in X for q4
    Y_E5 : list of numpy arrays (size : nsamples x n_print_times)
        List of estimated MRS signals for each parameter set contained in X for q5 
    '''

    #1. Retrieve the static configuration parameters
    workdir = Parameters['workdir'] # String, Directory for Hydrus files to be stored
    exe = Parameters['exe'] # String, Path to and name of the Hydrus executable file
    desc = Parameters['desc'] # String, Description of the case
    tprint_fn = Parameters['tprint_fn'] # String, Path to and name of the file containing print times
    atm_fn = Parameters['atm_fn'] # String, Path to and name of the .in file containing BC in Hydrus format
    ker_fn = Parameters['ker_fn'] # String, name of the .csv file containing MRS kernel
    savefn_num = Parameters['savefn_num'] # String, number X in output file name : "GSA5_X_E1.csv"
    zmax_MRS = Parameters['zmax_MRS']# [cm] maximum depth of the non-simulated zone
    zmax = Parameters['zmax'] # [cm] Maximum depth of Phydrus profile
    ncells = Parameters['ncells'] #Number of Phydrus simulation cells
    layer_bound = [Parameters['layer_bound']] # list of float, boundary between Phydrus soil layers
    tinit = Parameters['tinit'] # Float, [h] Simulation initial time
    tmax = Parameters['tmax'] # Float, [h] Simulation final time
    dtinit = Parameters['dtinit'] # Float, [h] Initial time step value 
    dtmin = Parameters['dtmin'] # Float, [h] Minimal time step value
    dtmax = Parameters['dtmax'] # Float, [h] Maximal time step value
    maxit = Parameters['maxit'] # Float, [-] Maximal number of iteration 
    tolth = Parameters['tolth'] # Float, [-] Tolerance on theta variable
    tolh = Parameters['tolh'] # Float, [cm] Tolerance on h variable
    top_bc = Parameters['top_bc'] # Int, Top boundary condition type (2 : Atmospheric Boundary Condition with Surface Layer.)
    bot_bc = Parameters['bot_bc'] # Int, Bottom boundary condition type (2 : variable pressure head)
    nlayers_phydrus = Parameters['nlayers_phydrus'] # Int, Number of VG params layers 
    ncells_bedrock = Parameters['ncells_bedrock'] # Number of cells attritbuted to the bedrock zone
    nvalues_theta_MRS_bed = Parameters['nvalues_theta_MRS_bed'] #Number of values of theta_MRS in the bedrock zone
    ncells_per_value_thetaMRSbed = [Parameters['ncells_per_value_thetaMRSbed']] # Number of cells attributed to each values in theta_MRS_bedrock
    nCells_MRS_interp = Parameters['nLayers_MRS_interp'] # Number of cells in the always-saturated zone

    del_outfiles = Parameters['del_outfiles'] # Boolean, if True, delete folder containing full simulation results after computing MRS signals
    nSave = Parameters['nSave'] # int, print outputs in a file every nSave simulations


    #make MRS experiment data
    tprint = np.loadtxt(tprint_fn, delimiter=';')

    Nq = 100

    qmin = 60
    qmax = 2260
    
    qrange = np.logspace(np.log10(qmin), np.log10(qmax), Nq)

    iq = [0, 24, 49, 74, 99]
    qrange_numexp = []

    for idx in iq : 
        qrange_numexp.append(qrange[idx])
    
    nMes = tprint.shape[0]


    kernel = read_interp_ker_csv(ker_fn)

    #2. Calculate the model response on each sample
    npar_layer = 5 
    Y_E1 = []; Y_E2 = []; Y_E3 = []; Y_E4 = []; Y_E5 = []

    savefilename = 'GSA5_5_'+ savefn_num
    np.savetxt(savefilename+'_X.csv', np.array(X), delimiter = ',')

    for k, theta in enumerate(X) :
        # Initialize Phydrus simulation 

        start_time = time.time()

        ml = ps.Model(exe_name=exe, ws_name=workdir, name="model", description=desc, mass_units="mmol",
              time_unit="hours", length_unit="cm", print_screen=False)

        ml.add_time_info(tinit= tinit, tmax= tmax,
                 dt = dtinit, dtmin = dtmin, dtmax = dtmax,
                 print_times=True,
                 print_array=tprint)
        
        ml.add_waterflow(maxit = maxit, tolth = tolth, tolh = tolh,
                 top_bc = top_bc, bot_bc = bot_bc)

        #define VGparams : [theta_r theta_s alpha n Ksat lvg]
        VG_params = []
        m = ml.get_empty_material_df(n=nlayers_phydrus)

        for i in range(nlayers_phydrus):
            VG_params.append([theta[i*npar_layer], theta[i*npar_layer+1], theta[i*npar_layer+2], theta[i*npar_layer+3], 10**theta[i*npar_layer+4], 0.5])

        m.loc[1:nlayers_phydrus] = VG_params
        ml.add_material(m)

        # Load atmospheric data and add delta_h to bottom BC
        dh_piezo = theta[-1]
        atm = np.loadtxt(atm_fn, skiprows= 9, max_rows=tmax)
        atm[:,6] = atm[:,6] + dh_piezo
        #convert to dataframe
        columns = ['tAtm','Prec','rSoil','rRoot','hCritA','rB','hB','hT']
        atm = pd.DataFrame(atm, columns=columns)
        ml.add_atmospheric_bc(atm)

        # Add profile (arbitrary, will be replaced by copied profile)
        profile = ps.create_profile(bot=-zmax, dx=1, h=-70.5614)
        ml.add_profile(profile)

        # Write input files 
        ml.write_input()

        # Replace Profile file 
        shutil.copyfile('Data/PROFILE_Hydrus_hetero.DAT', workdir+'/PROFILE.DAT')
        # Modify IC to take into account dhpiezo
        
        
        prof = np.loadtxt(workdir+'/PROFILE.DAT', skiprows=5, max_rows=ncells)
        prof[:,2] = prof[:,2] + dh_piezo

        PROFILE_DAT = open(workdir+'/PROFILE.DAT', 'w')

        with contextlib.redirect_stdout(PROFILE_DAT):
            print('Pcp_File_Version=4')
            print('2')
            print('1  0.000000e+000  1.000000e+000  1.000000e+000')
            print('2 '+str(-zmax)+'  1.000000e+000  1.000000e+000')
            print(str(ncells)+' 0 0 0         x        h  Mat  Lay  Beta  Axz  Bxz  Dxz  Temp Conc SConc')

            for i in range(ncells):
                print(*[int(prof[i,0]), prof[i,1], prof[i,2], int(prof[i,3])
                    , int(prof[i,4]), int(prof[i,5])
                    , prof[i,6], prof[i,7], prof[i,8]])
            print('0')

        PROFILE_DAT.close()

        ml.simulate()

        # Load outputs 

        Out_fn = workdir + "/NOD_INF.OUT"

        n_print = tprint.shape[0] + 1  
        N_rows_per_t = 701

        N_rows_header = 12
        N_rows_between_out = 9

        PHYDRUS_outputs = []

        for i in range(n_print+1):
            PHYDRUS_outputs.append(np.loadtxt(Out_fn, skiprows=N_rows_header+i*(N_rows_per_t+N_rows_between_out), max_rows = N_rows_per_t))

        
        PHYDRUS_outputs = PHYDRUS_outputs[1:] #No MRS signal at t=0h
        PHYDRUS_outputs.pop(-1) # no MRS signal at simulation end time (here 3700h)
        
        # Catch non-convergent simulations and put NaN as signal value

        if(PHYDRUS_outputs[-1].shape[0] == 0):
                
                print('Simulation '+str(k)+'has crashed')
                nanarray = np.empty(nMes)
                nanarray[:] = np.nan

                Y_E1.append(nanarray)
                Y_E2.append(nanarray)
                Y_E3.append(nanarray)
                Y_E4.append(nanarray)
                Y_E5.append(nanarray)
                continue
        
        depth = makeZarray_sim_extra_layers_Phydrus(PHYDRUS_outputs, nCells_MRS_interp, zmax/100, zmax_MRS/100)
        
        layer_bound_idx = [0]
        layer_bound_m = np.array(layer_bound) / 100. # convert depth to m.

        for z in layer_bound_m:
            layer_bound_idx.append(np.where(depth > z)[0][0])

        layer_bound_idx.append(np.where(depth > zmax/100)[0][0]+nCells_MRS_interp-ncells_bedrock)
        theta_MRS_bedrock = []
        
        for i in range(nvalues_theta_MRS_bed):
            theta_MRS_bedrock.append(theta[nlayers_phydrus*npar_layer + nlayers_phydrus * 3]+i)


        theta_s_sapro = theta[npar_layer+1]

        theta_s = []; theta_r = []
        theta_i_percent_base = []; S_lim = []; theta_i_percent_max = []


        # theta_i_percent_base, S_lim, theta_i_percent_max

        for i in range(nlayers_phydrus):
            theta_r.append(theta[i*npar_layer])
            theta_s.append(theta[i*npar_layer+1])
            theta_i_percent_base.append(theta[3*i + 2*npar_layer])
            S_lim.append(theta[3*i + 2*npar_layer +1])
            theta_i_percent_max.append(theta[3*i + 2*npar_layer +2])
        
        theta_PHYDRUS, theta_MRS_tot, est_sig = estimate_MRS_sig_from_sim_v4(PHYDRUS_outputs, depth, kernel, ncells_bedrock, ncells_per_value_thetaMRSbed, 
                                                                         theta_MRS_bedrock, theta_s_sapro, theta_s, theta_r, layer_bound_idx,
                                                                           theta_i_percent_base, S_lim, theta_i_percent_max)

        
        # Make signal array corresponding to injected pulse at each t_print
        # One signal array for each selected pulse
        # q1
        est_sig_in_pulse = np.zeros(nMes)

        for i in range(len(est_sig)):
            est_sig_in_pulse[i] = est_sig[i][np.where(qrange==qrange_numexp[0])][0]
        
        Y_E1.append(est_sig_in_pulse)

        # q2
        est_sig_in_pulse = np.zeros(nMes)

        for i in range(len(est_sig)):
            est_sig_in_pulse[i] = est_sig[i][np.where(qrange==qrange_numexp[1])][0]
        
        Y_E2.append(est_sig_in_pulse)

        # q3
        est_sig_in_pulse = np.zeros(nMes)

        for i in range(len(est_sig)):
            est_sig_in_pulse[i] = est_sig[i][np.where(qrange==qrange_numexp[2])][0]
        
        Y_E3.append(est_sig_in_pulse)

        # q4
        est_sig_in_pulse = np.zeros(nMes)

        for i in range(len(est_sig)):
            est_sig_in_pulse[i] = est_sig[i][np.where(qrange==qrange_numexp[3])][0]
        
        Y_E4.append(est_sig_in_pulse)

        # q5
        est_sig_in_pulse = np.zeros(nMes)

        for i in range(len(est_sig)):
            est_sig_in_pulse[i] = est_sig[i][np.where(qrange==qrange_numexp[4])][0]
        
        Y_E5.append(est_sig_in_pulse)

        if del_outfiles:
            # delete simulation files to avoid disk space saturation
            shutil.rmtree(workdir)


        end_time = time.time()
        ctime = end_time-start_time
        print('Sim '+str(k+1)+'/'+str(len(X))+ ' done. Computation time: '+str(ctime))

        if(k%nSave == 0):
            np.savetxt(savefilename+'_E1.csv', np.array(Y_E1), delimiter = ',')
            np.savetxt(savefilename+'_E2.csv', np.array(Y_E2), delimiter = ',')
            np.savetxt(savefilename+'_E3.csv', np.array(Y_E3), delimiter = ',')
            np.savetxt(savefilename+'_E4.csv', np.array(Y_E4), delimiter = ',')
            np.savetxt(savefilename+'_E5.csv', np.array(Y_E5), delimiter = ',')


    return(Y_E1, Y_E2, Y_E3, Y_E4, Y_E5)
