#!/bin/bash

# Before running this script, please edit the phd_config.json file with the 
# correct information regarding paths to DMD binaries, parameters, and 
# TURBOMOLE.
# 
# If you are using a queuing system like SLURM or UGE, please fill in the 
# queuing system so that PHD3 can submit jobs to the system (automates). 
# Though, you will have the option to manually or automatically resubmit jobs.
#
# This script will create a directory in your home directory labeled .phd3.
# Inside will be where this software stores some options, configurations, and 
# logger information. It will not be necessary to do much with this directory
# unless you move your TURBOMOLE or DMD directories. 
#
# Additionally, if you will be using the built-in auto resubmitter, then you 
# will need to place a submit.j2 file in the .phd3 directory that will act 
# as a template for submission. Albeit, there are several places where you need 
# to supply variables for PHD3 to fill in. These are denoted by {{ option }}. 
# Feel free to view the sample submit.j2 in phd3/templates directory (UGE 
# system).
#
# PYTHON INSTALLS THE EXECUTABLES TO YOUR ~/.local/bin DIRECTORY SO PLEASE
# INCLUDE THIS IN YOUR PATH
# 
#==============================================================================

echo "
>>>> Installing PHD3 >>>>"
mkdir ~/.phd3

echo ">>>> Copying files to ${HOME}/.phd3"
cp phd_config.json ~/.phd3/
cp phd3/resources/logger_config.json ~/.phd3/

pip install --user ./ || exit 1
echo "export PATH=$HOME/.local/bin:\$PATH" >> ~/.bashrc
echo ">>>> Successfully installed PHD3
"

