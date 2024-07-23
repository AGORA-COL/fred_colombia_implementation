#!/bin/bash

### JOB NUMBER
PROCESS_NUMBER=$1
PROCESS_NUMBER=$((PROCESS_NUMBER))

### PATHS
#export HOME=/zine/HPC02S1/ex-dveloza/AGORA/apps

export HOME=/zine/HPC02S1/ex-dveloza
### SET ENV PATH
export PATH=/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin
export PATH=$HOME/mambaforge/bin:$PATH

### PATHS
export PERL5LIB=/zine/HPC02S1/ex-dveloza/mambaforge/lib/perl5/vendor_perl
export LD_LIBRARY_PATH=/zine/HPC02S1/ex-dveloza/mambaforge/lib/:/usr/lib/:/usr/lib64/
### FRED SETUP
export condor_local=$(realpath ./)
export FRED_HOME=$condor_local/FRED
echo "FRED_HOME is set to $FRED_HOME"

export PATH=${PATH}:${FRED_HOME}/bin
echo "PATH is set to $PATH"

echo "####################################################################"
echo "################### collect_calibration.R #########################"
echo "####################################################################"
Rscript collect_calibration.R; 