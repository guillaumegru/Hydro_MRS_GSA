from GSA5_func_forUQ_numexp import PHYDRUS_MRS_wrapper_numexp
from uqpylab import sessions
from numpy.random import uniform

#STEP 1 : start UQ session

#Uncomment and fill with your UQ[py]Lab credentials
#myToken =  # The user's token to access the UQCloud API
#UQCloud_instance =  # The UQCloud instance to use

# Start the session
mySession = sessions.cloud(host=UQCloud_instance, token=myToken)
# (Optional) Get a convenient handle to the command line interface
uq = mySession.cli
# Reset the session
mySession.reset()

#STEP2 : computational model 


Model1Opts = {
    'Type': 'Model', 
    'ModelFun':'GSA5_func_forUQ_numexp.PHYDRUS_MRS_wrapper_numexp',
    'Parameters' : {
        'workdir' : str(uniform())[2:6], 
        'exe' : 'hydrus.exe',
        'desc' :  'test',
        'tprint_fn' : 'tprint.csv',
        'atm_fn' : 'Data/ATMOSPH_fromhydrus.IN',
        'ker_fn' : 'GSA5_ker.csv',
        'savefn_num' : '4',
        'zmax' : 700 , # [cm]
        'zmax_MRS' : 5000, #[cm]
        'ncells' : 701, 
        'layer_bound' : [210],
        'tinit' : 0,
        'tmax' : 3700,
        'dtinit' : 0.005,
        'dtmin' : 0.001,
        'dtmax' : 0.01,
        'maxit' : 50,
        'tolth' : 0.01,
        'tolh' : 1,
        'top_bc' : 2, # Int, Top boundary condition type (2 : Atmospheric Boundary Condition with Surface Layer.)
        'bot_bc' : 2, # Int, Bottom boundary condition type (2 : variable pressure head)
        'nlayers_phydrus' : 2, # Number of layers in the simulated zone
        'nLayers_MRS_interp' : 86, # Int, Number of layers in the zone which is always staturated 
        'nvalues_theta_MRS_bed' : 1, #Int, Number of values for theta_MRS in the non-simulated zone
        'ncells_bedrock' : 70,#Int, Number of layers of the always saturated zone associated with the bedrock zone 
        'ncells_per_value_thetaMRSbed' : 70, #Int, Number of layers in the interpolated kernel, associated with one theta_MRS value 

        'nSave' : 1, # int, print outputs in a file every nSave simulations

        'del_outfiles' : True # Boolean, if True, delete folder conaining full simulation results after computing MRS signals

    }
}

#create a model object
myModel = uq.createModel(Model1Opts)

# STEP 3 : probabilistic model input 

InputOpts = {
    "Marginals": [
        # Soil Layer 
        {"Type": "Uniform",
         "Parameters": [0.00, 0.10] # tr1
        },
        {"Type": "Uniform",
         "Parameters": [0.30, 0.50] #ts1
        },
        {"Type": "Uniform",
         "Parameters": [1e-2, 0.15] # alpha1
        },
        {"Type": "Uniform",
         "Parameters": [1.1, 3] # n1
        },
        {"Type": "Uniform",
         "Parameters": [-1., 2.18] # log10Ks1
        },
        # Saprolite Layer
        {"Type": "Uniform",
         "Parameters": [0.00, 0.04] # tr2
        },
        {"Type": "Uniform",
         "Parameters": [0.08, 0.20] #ts2
        },
        {"Type": "Uniform",
         "Parameters": [1e-2, 0.15] # alpha2
        },
        {"Type": "Uniform",
         "Parameters": [1.1, 3] # n2
        },
        {"Type": "Uniform",
         "Parameters": [-1, 0.18] # log10Ks2
        },
        # MRS undetectable water parameters
        # Soil layer
        {"Type": "Uniform",
         "Parameters": [0.20, 0.45] # theta_u_percent_base
        },
        {"Type": "Uniform",
         "Parameters": [0.40, 0.50] # S_lim
        },
        {"Type": "Uniform",
         "Parameters": [0.50, 0.70] # theta_u_percent_max
        },
        # Saprolite layer
        {"Type": "Uniform",
         "Parameters": [0.20, 0.45] # theta_u_percent_base
        },
        {"Type": "Uniform",
         "Parameters": [0.40, 0.50] # S_lim
        },
        {"Type": "Uniform",
         "Parameters": [0.50, 0.70] # theta_u_percent_max
        },

        # Non simulated zone parameters
        {"Type": "Uniform",
         "Parameters": [0.01, 0.05] # tMRS_Bedrock
        },
        # delta_h piezo
        {"Type": "Uniform",
         "Parameters": [-20.0, 20.0] 
        },
    ]
}

myInput = uq.createInput(InputOpts)

#STEP 4 : sample parameter distribution and make simulations 

X = uq.getSample(myInput, 1000)
Y_E1, Y_E2, Y_E3, Y_E4, Y_E5 = uq.evalModel(myModel,X) #, Y_th1, Y_th2, Y_th3,...

mySession.quit()


