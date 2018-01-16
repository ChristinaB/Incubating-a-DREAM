#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 10:40:32 2017

@author: barnhark
"""

# herbie driver

import sys
from subprocess import call
from yaml import load
import numpy as np

# set files and directories used to set input templates.
# Files and directories.
input_file = 'inputs.txt'
input_template = 'inputs_template.txt'

# Use `dprepro` (from $DAKOTA_DIR/bin) to substitute parameter
# values from Dakota into the SWASH input template, creating a new
# inputs.txt file.
call(['dprepro', sys.argv[1], input_template, input_file])

# now prepare to run model.
# load the params file to get the correct file names
with open(input_file, 'r+') as f:
    params = load(f)

# get parameter values from params dict
X1 = params['X1']
X2 = params['X2']

# evaluate herbie
herbie = -1 * (np.exp(-1.0 * (X1 - 1)**2) + np.exp(-0.8*(X1 + 1)**2) - 0.05*np.sin(8.*(X1 + 0.1))  +
           np.exp(-1.0 * (X2 - 1)**2) + np.exp(-0.8*(X2 + 1)**2) - 0.05*np.sin(8.*(X2 + 0.1)))

## calculate the residual # this could also be done by dakota by providing a data file
#value_out = herbie + 2.1230702895623965

value_out = herbie

# write out residual to the sys argument where the file dakota will look for is located
with open(sys.argv[2], 'w') as fp:
    fp.write(str(value_out)+'\n')
