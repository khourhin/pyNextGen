#! /usr/bin/env python3

import os
import sys
import subprocess

# Identify the ipynb and store arguments
ipynb = sys.argv[1]
args = sys.argv[2:]
os.environ['NB_ARGS'] = " ".join(args)

# Execute the notebook with the modified environment
cmd = ['jupyter', 'nbconvert', '--to', 'pdf', ipynb, '--execute']
subprocess.run(cmd)
