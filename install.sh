#!/bin/sh

# Can be installed on any system
# POSIX compliant

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
[ -d "${XDG_CONFIG_HOME:-$HOME/.config}" ] && PHDDIR="${XDG_CONFIG_HOME:-$HOME/.config}/phd3" || PHDDIR="${HOME}/.phd3"

# Make the phddirectory
[ -d  "$PHDDIR" ] || mkdir $PHDDIR

echo ">>>> Copying files to ${PHDDIR}"
cp phd_config.json ${PHDDIR}/
cp phd3/resources/logger_config.json ${PHDDIR}/

[ -f "${HOME}/.turbopy/submit.j2" ] && cp ${HOME}/.turbopy/submit.j2 ${PHDDIR} && echo "Submit file installed" || echo "Not submit file installed"

# Uninstall conflicting modules
pip list --format=columns | grep turbopy &> /dev/null && pip uninstall turbopy && echo "Remove all turbopy scripts from ~/.local/bin"
pip list --format=columns | grep dmdpy &> /dev/null && pip uninstall dmdpy && echo "Remove all turbopy scripts from ~/.local/bin"

# Actually load
pip list --format=columns | grep phd3 &> /dev/null || pip install --user ./ || echo "Failed to install using pip" && exit 1

echo "
Please 'export PATH=\$HOME/.local/bin:\$PATH' in your profile so 
that you have access to the binaries
" 

echo ">>>> Successfully installed PHD3
"

