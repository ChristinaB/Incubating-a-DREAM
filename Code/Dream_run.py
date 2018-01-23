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

#Function to write to DHSVM config file simulated inputs from DREAM.
def change_setting(config_file, setting_name, new_value, occurrence_loc='g'):
    sed_cmd = "sed -i 's:{setting_name} = .*:{setting_name} = {new_value}::{occurrence_loc}' {config_file}"
    sed_cmd = sed_cmd.format(setting_name = setting_name, new_value = new_value
                             , config_file = config_file
                             , occurrence_loc = occurrence_loc)
    return subprocess.call(sed_cmd, shell=True)

#Files and their input paths
config_file = 'Input.sauk.Livneh2013.dynG.rbm'
streamflow_only = 'output/Streamflow.Only'
validation_csv = 'validation.csv'
dhsvm_cmd = 'DHSVM3.1.3 ' + config_file

class Dream_run_setup(object):
    def __init__(self):

        self.params = [spotpy.parameter.Uniform('Lateral_Conductivity',low=0.001000 , high=0.003000,  optguess=0.001500),
                       spotpy.parameter.Uniform('Depth_Threshold',low=0.0200000000 , high=0.0900000000,  optguess=0.0500000000),
                       spotpy.parameter.Uniform('Maximum_Infiltration',low=0.010 , high=0.020,  optguess=0.017)
                       ]
        evaluation = pd.read_csv(validation_csv)
        self.evals = evaluation['value'].values                   
        print(len(self.evals))
        
    def parameters(self):
        return spotpy.parameter.generate(self.params)

