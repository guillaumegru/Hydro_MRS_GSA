"""
Miscelaneous useful functions

@author: Guillaume Gru
"""

import matplotlib.pyplot as plt

def plot_theta_z_simpler(depth, theta_soilsapro, theta_bedrock, nlayers_bedrock, n_theta_MRS_values_bedrock,
                 leg_pos = 'lower right', savefig = False, colors = ['#1f77b4', 'g', '#9467bd', '#7f7f7f', 'r'], 
                 labels = ['simu zone'], xlabel = r'$\theta$ [-]', linestyles = ['-', '--'], pngtitle = 'theta_z.png'):
    '''
    Plot theta(z) profile
    Considering extra layers for MRS signal computing
    Parameters
    ----------
    depth : ndarray of float64
        depth array.
    theta_soilsapro : ndarray of float64
        water content distribution above the bedrock.
    theta_bedrock : ndarray of float64
        water content distribution int the bedrock.
    nlayers_bedrock : int
        number of layers in the bedrock 
    n_theta_MRS_values_bedrock : int
        number of theta_MRS_bedrock values
    leg_pos : string, optional
        legend position.
        The default is 'lower right'
    savefig : boolean, optional
        save the figure.
        The default is False
    colors : list of strings, optional
        plot lines colors.
        The default is ['#1f77b4', 'g', '#9467bd', '#7f7f7f', 'r']
    labels : list of strings, optional
        line labels.
        The default is ['simu_zone']
    xlabel : string, optional
        xlabel.
        The default is r'$\theta$ [-]'
    linestyles : list of strings, optional
        linestyles options for pyplot plot.
        The default is ['-', '--']
    pngtitle : string, optional
        title of the png file to save (if savefig==True).
        The default is 'theta_z.png'

    Returns
    -------
    None.

    '''

    fontsize = 15

    k=0
    for i, theta in enumerate(theta_soilsapro):
        plt.plot(theta, depth[:theta.shape[0]], linewidth = 2, label= labels[i], color = colors[0], linestyle = linestyles[i])
        k+=1
    plt.ylim(max(depth), min(depth))
    plt.grid()
    
    for j, theta_b in enumerate(theta_bedrock):
        depth_MRS_toplot = [depth[-(nlayers_bedrock+1)], depth[-(nlayers_bedrock+1)]+0.001]
        theta_MRS_toplot = [theta_soilsapro[j][-1], theta_b[0]]

        for i in range(n_theta_MRS_values_bedrock-1):
            depth_MRS_toplot.append(depth[-(nlayers_bedrock+1)+i+1])
            depth_MRS_toplot.append(depth[-(nlayers_bedrock+1)+i+1]+0.001)

            theta_MRS_toplot.append(theta_b[i])
            theta_MRS_toplot.append(theta_b[i+1])
        

        depth_MRS_toplot.append(depth[-1])
        theta_MRS_toplot.append(theta_b[-1])
        plt.plot(theta_MRS_toplot, depth_MRS_toplot, color = 'black', label = labels[len(theta_bedrock)+j], linewidth = 2, linestyle = linestyles[j])
    

    plt.legend(loc = leg_pos, fontsize = fontsize+2)
    plt.xlabel(xlabel, fontsize = fontsize)
    plt.ylabel('depth [m]',fontsize = fontsize)
    plt.tick_params(labelsize =fontsize)


    if(savefig):
        plt.savefig(pngtitle, format = 'png', dpi=300, bbox_inches='tight')
        
        









