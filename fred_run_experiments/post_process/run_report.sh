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


# Overriding defaults based on the number of provided arguments
calibration_label="${1}"
ss="${2}"
outdir_simul="${3}"
outdir_md="${4}"
asymp_inf="${5}"
fm_ef="${6}"
ksus="${7}"
stage="${8}"
var_in="${9}"
vax_in="${10}"
report_label="${11}"

echo "####################################################################"
echo "################### get_fit_icu_duration.R #########################"
echo "####################################################################"
Rscript get_fit_icu_duration.R $calibration_label $ss $outdir_simul $asymp_inf $fm_ef $ksus $stage $var_in $vax_in $report_label; 

echo "##########################################################"
echo "################### get_report.R #########################"
echo "##########################################################"
Rscript get_report.R $calibration_label $ss $outdir_simul $asymp_inf $fm_ef $ksus $stage $var_in $vax_in $report_label; 

echo "######################################################"
echo "################### get_rt.R #########################"
echo "######################################################"
Rscript get_rt.R $calibration_label $ss $outdir_simul $asymp_inf $fm_ef $ksus $stage $var_in $vax_in; 
