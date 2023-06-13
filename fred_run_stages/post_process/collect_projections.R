#!/usr/bin/env Rscript
##=========================================#
## Author: Guido Espana
## Collect data from param sweep
## Year: 2019
## 
## requires:
##          densimtools library
##=========================================#
## User's input----------------
##=========================================#
library(tidyverse)
library(fredtools)

####################################
## FRED enviromental variables
####################################
 
Sys.setenv(FRED_HOME='/shared/home/Azure-ASIS/FRED')
Sys.setenv(FRED_RESULTS='/shared/home/Azure-ASIS/FRED_Implementation/FRED_results')
Sys.setenv(PATH='/bin:/shared/home/Azure-ASIS/anaconda3/condabin:/usr/local/cuda/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/usr/local/bin:/usr/local/sbin:/shared/home/Azure-ASIS/.local/bin:/shared/home/Azure-ASIS/bin:/shared/home/Azure-ASIS/FRED/bin:/shared/home/Azure-ASIS/.local/bin:/shared/home/Azure-ASIS/bin')
Sys.setenv(scratch_dir='/shared/home/Azure-ASIS/FRED_Implementation/scratch')

fred_home = Sys.getenv('FRED_HOME')
fred_defaults = sprintf("%s/input_files/defaults", fred_home)

##========================================#
## Fix files -------------
##========================================#
check_finished_jobs <- function(params_in){

    job_id = params_in$job_id[1] ## This should be only the job parameters

    fred_dir = system(sprintf("fred_find -k %s", job_id), intern = TRUE)
    data_dir = file.path(fred_dir,
                         "DATA", "OUT")

    if(!is.null(attr(fred_dir, "status"))){
        return(FALSE)
    }

    fred_status = system(sprintf("fred_status -k %s", job_id), intern = TRUE)
    if(unlist(strsplit(fred_status, "\\s+"))[1] != "FINISHED"){
        print(sprintf("STATUS is %s\n", unlist(strsplit(fred_status, "\\s+"))[1]))
        return(FALSE)
    }    

    fred_csv = sprintf("fred_csv -k %s -v I", job_id)
    fred_csv_st = system(fred_csv, intern = T)
    if(!is.null(attr(fred_csv_st, "status"))){
        ##system(sprintf("rm -rf %s", tmp_out), intern = T)
        return(FALSE)
    }

    job_df = read.csv(text = fred_csv_st, stringsAsFactors = FALSE)
    
    job_df = drop_na(job_df)
    if(!(nrow(job_df) == params_in$days[1] & !is.na(job_df$I_mean[1]))){        
        return(FALSE)
    }
    
    if(!file.exists(file.path(data_dir,'AgeGroupsData_mean.csv')) & nrow(job_df) > 1){
        return(FALSE)
    }
    return(TRUE)
}

process_job_data <- function(job_id){
    fred_dir = system(sprintf("fred_find -k %s", job_id), intern = TRUE)
    data_dir = file.path(fred_dir,
                         "DATA", "OUT")
    fred_csv = sprintf("fred_csv -k %s", job_id)

    fred_csv_st = system(fred_csv, intern = T)
    if(!is.null(attr(fred_csv_st, "status"))){
        ##system(sprintf("rm -rf %s", tmp_out), intern = T)
        return(FALSE)
    }

    job_df = read.csv(text = fred_csv_st, stringsAsFactors = FALSE)
    
    print(sprintf("Processing job %s", job_id))

    ##========================================#
    ## Age stuff-------------
    ##========================================#
    age_df = read.csv(file.path(data_dir,'AgeGroupsData_mean.csv'), stringsAsFactors = F)    
    cols_age = colnames(age_df %>% dplyr::select(-Day))
    job_df = dplyr::left_join(job_df, age_df, by = "Day")        
    job_df$Date = as.Date('2020-01-01') + job_df$Day
    
    cols_output = c('Day', 'Week', 'Year', 'Runs', 'N_mean', 'N_std', 'PrevInf_mean',
                    'C_mean', 'C_std', 'Cs_mean', 'Cs_std','CFR_mean', 'Is_mean', 'I_mean',
                    'CFR_std', 'CF_mean', 'CF_std','Nursing_Home_mean', 
                    'Nursing_Home_CF_mean', 'Wrk_mean', 'Sch_mean', 'H_mean', 'Nbr_mean',
                    'RR_mean', 'RR_std', 'AR_mean', 'AR_std','Chosp_mean','Phosp_mean',
                    'H_sheltering_mean', 'N_sheltering_mean')
    
    if('CF_1_mean' %in% colnames(job_df)){
        cols_output = c(cols_output, "CF_1_mean", "C_1_mean", "Cs_1_mean", "Phosp_1_mean", "Chosp_1_mean", "RR_1_mean")
    }
    if('CF_2_mean' %in% colnames(job_df)){
        cols_output = c(cols_output, 'CF_2_mean',  "C_2_mean", "Cs_2_mean", "Phosp_2_mean", "Chosp_2_mean", "RR_2_mean")
    }
    if('CF_3_mean' %in% colnames(job_df)){
        cols_output = c(cols_output, 'CF_3_mean',  "C_3_mean", "Cs_3_mean", "Phosp_3_mean", "Chosp_3_mean", "RR_3_mean")
    }
    if('CF_4_mean' %in% colnames(job_df)){
        cols_output = c(cols_output, 'CF_4_mean',  "C_4_mean", "Cs_4_mean", "Phosp_4_mean", "Chosp_4_mean", "RR_4_mean")
    }
    if('CF_5_mean' %in% colnames(job_df)){
        cols_output = c(cols_output, 'CF_5_mean',  "C_5_mean", "Cs_5_mean", "Phosp_5_mean", "Chosp_5_mean", "RR_5_mean")
    }
    if('CF_6_mean' %in% colnames(job_df)){
        cols_output = c(cols_output, 'CF_6_mean',  "C_6_mean", "Cs_6_mean", "Phosp_6_mean", "Chosp_6_mean", "RR_6_mean")
    }
    if('CF_7_mean' %in% colnames(job_df)){
        cols_output = c(cols_output, 'CF_7_mean',  "C_7_mean", "Cs_7_mean", "Phosp_7_mean", "Chosp_7_mean", "RR_7_mean")
    }
    if(!('V_mean' %in% colnames(job_df))){
        job_df$V_mean = 0
    }
    if(!('Vtd_mean' %in% colnames(job_df))){
        job_df$Vtd_mean = 0
    }
    
    cols_output = c(cols_output, cols_age, 'V_mean', 'Vtd_mean')
    
    return(list(job_df = job_df[,cols_output]))
    
}

##========================================#
## Collect for covid --------------------
##========================================#
projection_label = 'omicron_lineages_BQX'
state_input = 11001
asymp_infectivity_in = 1.0
face_mask_transmission_efficacy_in = 0.73
kids_susceptibility_in = 10
variants_in = 1
vaccination_in = 70
args = (commandArgs(TRUE))
if(length(args) >= 1){
    projection_label = args[1]
    if(length(args) >= 2){
        state_input = as.numeric(args[2])
        if(length(args) >= 3){
            asymp_infectivity_in = as.numeric(args[3])
            if(length(args) >= 4){
                face_mask_transmission_efficacy_in = as.numeric(args[4])
                if(length(args) >= 5){
                    kids_susceptibility_in = as.numeric(args[5])
                    if(length(args) >= 6){
                        variants_in = as.numeric(args[6])
                        if(length(args) >= 7){
                            vaccination_in = as.numeric(args[7])
                        }
                    }
                }
            }    
        }
    }
}
Sys.setenv(scratch_dir='/shared/home/Azure-ASIS/FRED_Implementation/scratch')
scratch_dir = Sys.getenv('scratch_dir')
simdir_name = sprintf('FRED_%d_projections_asymp_%.2f_fm_%.2f_ksus_%.2f_var_%.0f_vax_%03d_mov_%s',state_input, asymp_infectivity_in,face_mask_transmission_efficacy_in, kids_susceptibility_in, variants_in, vaccination_in, projection_label)
simdir = file.path(scratch_dir, simdir_name)

if(!dir.exists(file.path(getwd(), '../output','SHORT_FORECAST'))){
    dir.create(file.path(getwd(), '../output','SHORT_FORECAST'), recursive = T)
}

fred_results_dir = file.path(simdir,"FRED_RESULTS")
Sys.setenv(FRED_RESULTS=fred_results_dir)


output_dir = file.path(getwd(), '../output','SHORT_FORECAST',sprintf("%s_%s", simdir_name, "out"))
fred_outputs = 'fred_output.csv'
params_file = sprintf("%s/%s",simdir, "FRED_parameters.csv")

fredtools::fred_gather_data(params=params_file, outdir=output_dir, outfile=fred_outputs,
                            FUN=check_finished_jobs, FUN2=process_job_data, rm_out = TRUE, appendToFile = FALSE)    
