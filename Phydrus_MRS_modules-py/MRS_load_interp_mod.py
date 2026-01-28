"""
This python library includes the MRS forward model and the coupling between the hydrological model and the MRS forward model.

@author: Guillaume Gru
"""

import numpy as np
import matplotlib.pyplot as plt 

def loadMRSkernel(filename):
    '''
    Load MRS kernel from .mrm file.

    Parameters
    ----------
    filename : str
        name of the Samovar kernel data file in the same folder as python file.

    Returns
    -------
    info : list of float64
        Samovar parameters for this kernel
         -> incl [deg]     : local geomagnetic field inclination 
         -> antt [-]       : type of the loop setup (1 : cirlce, 2: square, 3: cicrle 8 , 4: square 8)
         -> length [m]     : side length
         -> frequence [Hz] : Emission frequency
    model : list of float64
        resistivity model: depth (m) / resistivity (Ohm.m).
    depth : ndarray of float64
        array containing the values of depth for the kernel [m].
    pulse : ndarray of float64
        array containing the values of injected pulse amplitude for the kernel [A.ms].
    ker : ndarray of float64
        Samovar kernel [nV/m].

    '''
    
    row = 0
    
    info = np.loadtxt(filename+'.mrm', skiprows=row, max_rows = 1)
    row +=1
    
    model = np.loadtxt(filename+'.mrm', skiprows=row, max_rows = 7)
    row += 7
    
    data = np.loadtxt(filename+'.mrm', skiprows=row, max_rows = 100)
    depth = data[:,1]
    pulse = data[:,0]
    
    row+=100
    
    ker_0 = np.loadtxt(filename+'.mrm', skiprows=row)
    ker = np.array(ker_0[:,0] + 1j*ker_0[:,1], dtype='complex')
    ker = ker/10 # nV /10000 * 1000 (normalized by 1e5 and for 1 mm thick)
    ker = np.transpose(ker.reshape((100,100)))
    
    return info, model, depth, pulse, ker

def plotMRSkernel(pulse, depth, ker, Qmin, Qmax, title, titlepng = "MRS_kernel", savefig = False, fontsize = 20, xlabel = 'Pulse [A.ms]', ylabel = 'Depth [m]', cmap = 'jet'):
    '''
    Plot MRS Kernel kappa(q,z)

    Parameters
    ----------
    pulse : ndarray of float64
        pulse distribution.
    depth : ndarray of float64
        depth distribution.
    ker : ndarray of float64
        kernel to print.
    Qmin : float
        pulse min value to print.
    Qmax : float
        pulse max value to print.
    title : string
        figure title.
    titlepng : string
        png filename (if savefig = True).
    savefig : boolean, optional
        save the figure. 
        The default is False.
    fontsize : float, optional
        font size for figure. 
        The default is 20.
    xlabel : string, optional
        label for xaxis. 
        The default is 'Pulse [A.ms]'
    ylabel : string, optional
        label for yaxis.
        The default is 'Depth [m]'
    cmap : string, optional
        pyplot color map.
        The default is 'jet'.

    Returns
    -------
    None.

    '''
    
    fig = plt.figure(figsize=(10,5))
    ax = fig.add_subplot(111)
    im = ax.pcolor(pulse, depth, np.transpose(np.absolute(ker)), cmap = cmap)
    plt.ylim(max(depth), min(depth))
    plt.xlim(Qmin,Qmax)
    plt.xscale('log')
    cbar = fig.colorbar(im, ticks = [50,100,150,200,250,300])
    ax.tick_params(labelsize = fontsize-3)
    ax.set_title(title, fontsize = fontsize)
    ax.set_xlabel(xlabel, fontsize = fontsize)
    ax.set_ylabel(ylabel, fontsize = fontsize)
    cbar.ax.tick_params(labelsize=fontsize)

    if(savefig):
        plt.savefig(titlepng+'.png', format = 'png', transparent = True, dpi=300, bbox_inches='tight')
    
    
def read_interp_ker_csv(filename):
    '''
    Read interpolated kernel from csv file.

    Parameters
    ----------
    filename : string
        Name of the interpolated kernel csv file.
    
    Returns
    -------
    ker : numpy array
        Interpolated SAMOVAR Kernel
    '''
    
    ker_0 = np.loadtxt(filename, delimiter = ';', skiprows=1, usecols = (0,1))
    ker = np.array(ker_0[:,0] + 1j*ker_0[:,1], dtype='complex')
    
    info = np.loadtxt(filename, delimiter = ';', skiprows=1, max_rows = 2, usecols = (2))
    
    ker = ker.reshape((int(info[0]), int(info[1])))
    
    return ker


def makeZarray_sim_extra_layers_Phydrus(Phydrus_outputs, nlayers_MRS, extraLbegin, extraLend):
    '''
    Make depth array from Phydrus simulation domain and extra layers
    
    Parameters
    ----------
    Phydrus_outputs : list of ndarrays
        Phydrus simulation results, used to extract depth values.
    nlayers_MRS : int
        number of extra layers for MRS signal estimation.
    extraLbegin : float64
        depth of the end of Phydrus simulation domain (beginning of extra layers).
    extraLend : float646
        depth of the end of the extra layers.

    Returns
    -------
    depth : ndarray of float64.
        depth array containing Phydrus simulation domain and extra layers
    '''
    
    depth_Phydrus_m = (-Phydrus_outputs[0][1:,1])/100 # Use Phydrus output file to get grid depth values and convert it to m

    depth_size = depth_Phydrus_m.shape[0] + nlayers_MRS 
    depth = np.zeros(depth_size)

    depth[:depth_Phydrus_m.shape[0]] = depth_Phydrus_m

    extra_l_thick = (extraLend - extraLbegin)/ nlayers_MRS # thickness of extra layers

    depth[depth_Phydrus_m.shape[0]:depth_Phydrus_m.shape[0]+nlayers_MRS] = np.linspace(extraLbegin + extra_l_thick,extraLend,nlayers_MRS)
    
    return(depth)


def compute_theta_u_percent(theta, theta_s, theta_r, layer_bound_idx, theta_u_percent_base, S_lim, theta_u_percent_max):
    '''
    Compute MRS undetectable water content fraction. 
    Version 2 : takes into account the effects of saturation on MRS undetectable water content fraction.
    Ref: Boucher et al., 2011 - The detectability of water by NMR considering the instrumental dead time – a laboratory analysis of unconsolidated materials
        doi : https://doi.org/10.3997/1873-0604.2010056

    Parameters
    ----------
    theta : ndarray of float64 (Size: Nz + nlayers_MRS)
        Phydrus simulation outputs : water content distributions.
    theta_s : list of float64 (Size: nLayers_simulated_zone)
        List of saturated water content parameters for the simulated zone.
    theta_r : list of float64 (Size: nLayers_simulated_zone)
        List of residual water content parameters for the simulated zone.
    layer_bound_idx : list of int (nLayers_simulated_zone + 2)
        List of indices of boundaries between layers.
    theta_u_percent_base : list of float64 in [0,1]
        Base value of theta_u_percent: when saturation is above S_lim.
    S_lim : list of float64 in [0,1]
        Threshold saturation below which theta_u_percent increases.
    theta_u_percent_max : list of float64 in [0,1]
        Maximum fraction of MRS undetectable water.

    Returns
    -------
    theta_u_percent : ndarray of float64 (Size : Nz + NlayersMRS)
        Fraction of MRS undetectable water.
    '''

    theta_u_percent = np.zeros(theta.shape[0])

    # Compute eff_sat array : (theta-theta_r)/(theta_s-theta_r)
    eff_sat = np.zeros(theta.shape[0])
    for i in range(len(layer_bound_idx)-1):
        eff_sat[layer_bound_idx[i]:layer_bound_idx[i+1]] = (theta[layer_bound_idx[i]:layer_bound_idx[i+1]]-theta_r[i]) / (theta_s[i]-theta_r[i])

    k=1 #Iterative int for layer number 
    
    for i, s in enumerate(eff_sat):
        if(i>layer_bound_idx[k]): 
            k += 1
        
        if(s > S_lim[k-1]): theta_u_percent[i] = theta_u_percent_base[k-1]
        else :
            # linear from theta_u_percent(0) = theta_u_percent_max, to theta_u_percent(Slim) = theta_u_percent_base 
            theta_u_percent[i] = (theta_u_percent_base[k-1]-theta_u_percent_max[k-1])/S_lim[k-1] * s + theta_u_percent_max[k-1]

    return eff_sat, theta_u_percent

def estimateSignal_v4(kernel, theta_z, depth, theta_u_percent, ncells_bedrock, ncells_per_value_thetaMRSbed,
                      theta_MRS_bedrock):
    '''
    Estimate MRS signal from water content distribution and Samovar (interpolated) kernel
    Version 4 : takes into account an undetectable water content, and its dependence on saturation  
        theta_MRS = theta * (1 - theta_u_percent)
        
    No theta_u_percent value for bedrock, the parameter theta_MRS_bedrock already is equal to theta_s - theta_u_percent_base in the bedrock

    Parameters
    ----------
    kernel : Numpy array of complex128
        Samovar or interpolated kernel. ([Nv/m], size = Nq x Nz)
    theta_z : Numpy array of float64
        Water content distribution above bedrock. ([-], size = Nz)
    depth : Numpy array of float64
        Depth array ([m] size = Nz).
    theta_u_percent : list of float64 (size = Nz - ncells_bedrock)
        Fraction of MRS undetectable water (in cells above the bedrock)
    ncells_bedrock : int
        Number of cells attributed to bedrock (treated separately because always saturated)
    theta_MRS_bedrock : list of float in [0,1] (size < ncells_bedrock)
        List of theta_MRS_values in the bedrock zone 
    ncells_per_value_thetaMRSbed : list of int
        Number of cells for each theta_MRS value in theta_MRS_bedrock.
    Returns
    -------
    signal : Numpy array of float64
        Estimated signal (for all pulses considered in kernel). ([Nv], size = Nq)
 
    '''

    if(ncells_bedrock!=sum(ncells_per_value_thetaMRSbed)):
        print('Error : the sum of the values in ncells_per_value_thetaMRSbed should be equal to ncells_bedrock')
        return(-1)
    
    # construction of the theta_MRS vector
    ncells_tot = depth.shape[0]
    theta_u_percent_fullarr = np.ones(ncells_tot)
    theta_u_percent_fullarr[:theta_u_percent.shape[0]] = theta_u_percent

    theta_bedrock_fullarr = np.zeros(ncells_tot)  

    theta_z_fullarr = np.zeros(ncells_tot)  
    theta_z_fullarr[:theta_z.shape[0]] = theta_z

    k = 0 #iteraror : index in bedrock zone of theta_z array
    for i, th in enumerate(theta_MRS_bedrock):
        for j in range(ncells_per_value_thetaMRSbed[i]):
            theta_bedrock_fullarr[-ncells_bedrock+k] = th
            k+=1

    theta_MRS = theta_z_fullarr * (1-theta_u_percent_fullarr) + theta_bedrock_fullarr
    
    thickness = np.diff(np.concatenate((np.array([0]),depth)))
    signal = np.absolute(kernel) @ (thickness*theta_MRS)
    
    return signal, theta_MRS

def estimate_MRS_sig_from_sim_v4(Phydrus_outputs, depth, kernel, ncells_bedrock, ncells_per_value_thetaMRSbed, theta_MRS_bedrock,
                                 theta_s_sapro, theta_s, theta_r, layer_bound_idx, theta_u_percent_base, S_lim, theta_u_percent_max):
    '''
    Use estimateSignal_v4 function to compute MRS signals from Phydrus simulation and extra layers.
    Version 4 : v2 implementation of MRS undetectable water content 
    Parameters
    ----------
    Phydrus_outputs : list of ndarrays
        Phydrus simulation results
    depth : Numpy array of float64
        Depth array ([m] size = Nz).
    kernel : Numpy array of complex128
        SAMOVAR or interpolated kernel. ([Nv/m], size = Nq x Nz)
    ncells_bedrock : int
        Number of cells attributed to bedrock (treated separately because always saturated)
    ncells_per_value_thetaMRSbed : list of int
        Number of cells for each theta_MRS value in theta_MRS_bedrock.
    theta_MRS_bedrock : list of float in [0,1] (size < ncells_bedrock)
        List of theta_MRS_values in the bedrock zone 
    theta_s_sapro : float64
        Saturated water content in the saprolite layer.
    theta_s : list of float64 (Size: nLayers_simulated_zone)
        List of saturated water content parameters for the simulated zone.
    theta_r : list of float64 (Size: nLayers_simulated_zone)
        List of residual water content parameters for the simulated zone.
    layer_bound_idx : list of int (nLayers_simulated_zone + 2)
        List of indices of boundaries between layers.
    theta_u_percent_base : list of float64 (Size: nLayers_simulated_zone)
        Base value of theta_u_percent: when saturation is above S_lim.
    S_lim : list of float64 (Size: nLayers_simulated_zone)
        Threshold saturation below which theta_u_percent increases.
    theta_u_percent_max : list of float64 (Size: nLayers_simulated_zone)
        Maximum fraction of MRS undetectable water.
    
    Returns
    -------
    theta_Phydrus : list of ndarray (Nq)
        water content distribution resulting from Hydrus simulation.
    theta_MRS_tot : list of ndarray (Nq)
        MRS-detectable water content distribution with additionnal layers.
    est_sig : list of ndarray  (Nq)
        estimated MRS signal from theta.
    '''
    theta_Phydrus = []
    theta_MRS_tot = []
    est_sig = []

    ncells_tot = depth.shape[0]
    
    for k in range(len(Phydrus_outputs)):
    
        theta_Phydrus.append(Phydrus_outputs[k][1:,3])
        
        theta_above_bed = np.zeros(ncells_tot-ncells_bedrock)
        theta_above_bed[:theta_Phydrus[k].shape[0]] = theta_Phydrus[k]
        theta_above_bed[theta_Phydrus[k].shape[0]:(ncells_tot-ncells_bedrock)] = theta_s_sapro

        eff_sat, theta_u_percent = compute_theta_u_percent(theta_above_bed, theta_s, theta_r, layer_bound_idx, theta_u_percent_base,
                                                            S_lim, theta_u_percent_max)
        
        # Estimate MRS signal
        
        signal, theta_MRS = estimateSignal_v4(kernel, theta_above_bed, depth, theta_u_percent, ncells_bedrock, ncells_per_value_thetaMRSbed,
                      theta_MRS_bedrock)

        est_sig.append(signal)
        
        theta_MRS_tot.append(theta_MRS)
    
    return theta_Phydrus, theta_MRS_tot, est_sig