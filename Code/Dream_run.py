from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import numpy as np
import pandas as pd
import spotpy
import os
from scipy import stats
import sys
import subprocess
import spotpy


#Function to write to DHSVM config file simulated inputs from DREAM.
def change_setting(config_file, setting_name, new_value, occurrence_loc='g'):
    sed_cmd = "sed -i 's:{setting_name} = .*:{setting_name} = {new_value}:{occurrence_loc}' {config_file}"
    sed_cmd = sed_cmd.format(setting_name = setting_name, new_value = new_value
                             , config_file = config_file
                             , occurrence_loc = occurrence_loc)
    return subprocess.call(sed_cmd, shell=True)

#Files and their input paths

config_file = 'Input.sauk.wrf.dynG.global.2006flood'
streamflow_only = 'output/DHSVMoutput_wrf_global_2006flood/Streamflow.Only'
validation_csv = 'validation.csv'
dhsvm_cmd = 'DHSVM3.1.3 ' + config_file 


class Dream_run_setup(object):
    def __init__(self):

        self.params = [spotpy.parameter.Uniform('Precipitation_Lapse_Rate',low=0.0001 , high=0.001,  optguess=0.0006),
                       spotpy.parameter.Uniform('Temperature_Lapse_Rate',low=-0.0025 , high=-0.008,  optguess=-0.0048)
                       ]
        evaluation = pd.read_csv(validation_csv)
        self.evals = evaluation['value'].values
        print(len(self.evals))
        
    def parameters(self):
        return spotpy.parameter.generate(self.params)
    
    #setting up simulation for location:12189500 with predefined params and writing to config file 
    def simulation(self,x):
        #write DREAM parameter input to config file.
        change_setting(config_file, "Precipitation Lapse Rate", str(round(x[0],5)))
        change_setting(config_file, "Temperature Lapse Rate", str(round(x[1],5)))
#         change_setting(config_file, "Overstory Monthly LAI    9"," ".join([str(round(x[3],3))]*12))
        #run DHSVM with modified parameters in config file
        retcode = subprocess.call(dhsvm_cmd, shell=True)
        print("Ran DHSVM: ", str(retcode))
        
        simulations=[]
        #read streamflow data from DHSVM output file
        with open(streamflow_only, 'r') as file_output:
            header_name = file_output.readlines()[0].split(' ')

        with open(streamflow_only) as inf:
            next(inf)
            date_q = []
            q_12189500 = []
            for line in inf:
                parts = line.split() 
                if len(parts) > 1:
                    date_q.append(parts[0])
                    q_12189500.append(float(parts[2])/(3600*3))

        Simulation_streamflow = pd.DataFrame({'x[0]':date_q, 'x[2]':q_12189500})
        Simulation_streamflow.columns = [header_name[0], header_name[2]]
        Simulation_streamflow = Simulation_streamflow[:-14]
        print(Simulation_streamflow.shape)
        simulations = Simulation_streamflow['12189500'].values 
        return simulations
    
    def evaluation(self):
        return self.evals.tolist()
    
    def objectivefunction(self,simulation,evaluation, params=None):
        model_fit = spotpy.objectivefunctions.nashsutcliffe(evaluation,simulation)
        print('Nashsutcliffe: ', model_fit)
        return model_fit

# Initialize the Dream Class
dream_run=Dream_run_setup()

# Create the Dream sampler of spotpy, al_objfun is set to None to force SPOTPY
# to jump into the def objectivefunction in the Dream_run_setup class with 
# nashsutcliffe as objectivefunction.
dream_sampler=spotpy.algorithms.dream(dream_run, dbname='DREAM_dhsvm', dbformat='csv',
                                alt_objfun=None)

#Select number of maximum repetitions
rep=4

# Select five chains and set the Gelman-Rubin convergence limit
nChains                = 4
convergence_limit      = 1.2
runs_after_convergence = 1

r_hat = dream_sampler.sample(rep,nChains=nChains,convergence_limit=convergence_limit, 
                       runs_after_convergence=runs_after_convergence, acceptance_test_option = 6)
r_hat = np.array(r_hat)

np.save(outfile, r_hat)
