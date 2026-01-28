# Code used to run locally the forward model computations 

import numpy as np
from GSA5_func_forUQ_numexp import PHYDRUS_MRS_wrapper_numexp
from numpy.random import uniform


Parameters = {}

Parameters['workdir'] = str(uniform())[2:6] # String, Directory for Hydrus files to be stored
Parameters['exe'] = 'hydrus.exe' # String, Path to and name of the Hydrus executable file
Parameters['desc'] = 'test' # String, Description of the case
Parameters['tprint_fn'] = 'tprint.csv' # String, Path to and name of the file containing print times
Parameters['atm_fn'] = 'Data/ATMOSPH_fromhydrus.IN' # String, Path to and name of the .in file containing BC in Hydrus format
Parameters['ker_fn'] = 'GSA5_ker.csv' # String, name of the .csv file containing MRS kernel
Parameters['savefn_num'] = 'test' # String, number X in output file name : "GSA5_X_E1.csv"
Parameters['zmax'] = 700 # [cm] maximum depth of the Phydrus profile 
Parameters['zmax_MRS'] = 5000 # [cm] maximum depth of the non-simulated zone
Parameters['ncells'] = 701 # int, number of Phydrus simulation cells 
Parameters['layer_bound'] = [210] # list of float, boundary between Phydrus soil layers 
Parameters['tinit'] = 0 # Float, [h] Simulation initial time
Parameters['tmax'] = 3700 # Float, [h] Simulation final time
Parameters['dtinit'] = 0.005 # Float, [h] Initial time step value 
Parameters['dtmin'] = 0.001 # Float, [h] Minimal time step value
Parameters['dtmax'] = 0.01 # Float, [h] Maximal time step value
Parameters['maxit'] = 50 # Float, [-] Maximal number of iteration 
Parameters['tolth'] = 0.01 # Float, [-] Tolerance on theta variable
Parameters['tolh'] = 1 # Float, [cm] Tolerance on h variable
Parameters['top_bc'] = 2 # Int, Top boundary condition type (2 : Atmospheric Boundary Condition with Surface Layer.)
Parameters['bot_bc'] = 2 # Int, Bottom boundary condition type (2 : variable pressure head)
Parameters['nlayers_phydrus'] = 2 # Int, Number of VG params layers 
Parameters['nLayers_MRS_interp'] = 86 #Int,  Number of cells in the always-saturated zone
Parameters['nvalues_theta_MRS_bed'] = 1 #Int, Number of values of theta_MRS in the bedrock zone
Parameters['ncells_bedrock'] = 74 # Number of cells attritbuted to the bedrock zone
Parameters['ncells_per_value_thetaMRSbed'] = 74 # Number of cells attributed to each values in theta_MRS_bedrock
Parameters['nSave'] = 1 # int, print outputs in a file every nSave simulations
Parameters['del_outfiles'] = False # Boolean, if True, delete folder containing full simulation results after computing MRS signals

X = np.loadtxt('GSA5_5_'+Parameters['savefn_num']+'_X.csv', delimiter=',')

Y_E1, Y_E2, Y_E3, Y_E4, Y_E5 = PHYDRUS_MRS_wrapper_numexp(X, Parameters)
