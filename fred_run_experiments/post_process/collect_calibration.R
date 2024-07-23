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
library(DirichletReg)
library(dplyr)
library(parallel)
library(readr)

####################################################
## local variables
####################################################
repo_name = 'fred_colombia_implementation'
AGORA_path = '/zine/HPC02S1/ex-dveloza/AGORA/apps'

####################################################
## FRED set cluster enviromental variables
####################################################
setwd(sprintf('%s/%s/fred_run_experiments/post_process', AGORA_path, repo_name))
FRED_results_path = sprintf('%s/FRED_results', AGORA_path)

Sys.setenv(FRED_HOME=sprintf('%s/FRED', AGORA_path))
Sys.setenv(FRED_RESULTS=FRED_results_path)
Sys.setenv(scratch_dir=sprintf('%s/%s/scratch', AGORA_path, repo_name))
Sys.setenv(PATH=sprintf('/zine/HPC02S1/ex-dveloza/mambaforge/bin:/bin/:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/usr/local/bin:/usr/local/sbin:%s/FRED/bin', AGORA_path))

##==============================================#
## load dirs -------------
##==============================================#
fred_home = Sys.getenv('FRED_HOME')
fred_defaults = sprintf("%s/fred_input_files/defaults", fred_home)

##========================================#
## Functions -------------
##========================================#
get_loglikelihood_poisson <- function(data_in, model_in){
    model_in[which(model_in == 0)] = 0.00000001
    ll_out = -sum(dpois(data_in, model_in, log = T), na.rm = T)
    return(ll_out)
}

get_loglikelihood <- function(data_in, model_in, reps_in = 1, lambda = 0.8) {
    # Basic Prior Settings
    shape.prior = 1e-3
    scale.prior = 1e3
    
    # Posterior shape based on the prior and the simulated model values
    shape.post = shape.prior + model_in
    # Adjusting the probability parameter to reflect the expected number of events (deaths)
    scale.post = scale.prior / (reps_in * scale.prior + 1)
    prob = 1 / (1 + scale.post)
    
    # Calculate the negative binomial log-likelihood
    nbinom_ll = -sum(dnbinom(data_in, size = shape.post, prob = prob, log = TRUE), na.rm = TRUE)
    
    # Adding a penalty for differences between actual deaths and simulated deaths
    # The penalty is lambda times the sum of squared differences
    penalty = lambda * sum((data_in - model_in)^2, na.rm = TRUE)
    
    # Total log-likelihood adjusted for the penalty
    total_loglikelihood = nbinom_ll + penalty

    return(total_loglikelihood)
}

get_agegroups_likelihood <- function(data_df){
    tmp_age_groups = unique(data_df$AgeGroup)
    total_age_LL = 0
    for(aa in tmp_age_groups){
        tmp_model_data = filter(data_df, AgeGroup == aa)
        total_age_LL = total_age_LL + get_loglikelihood(tmp_model_data$Deaths, tmp_model_data$CF_mean)
    }
    return(total_age_LL)
}

get_death_age_likelihood <- function(data_in, model_in){    
    return(-dmultinom(x=data_in, prob = model_in/sum(model_in), log = T))
}

get_loglikelihood_dominance <- function(data_in, model_in){
    ## Gather data
    data_ll = left_join(data_in, model_in, by = c('Date','Variant'),suffix = c('_data','_model'))

    model_mat = as.matrix(as.data.frame(dplyr::select(data_ll, Date, Variant, Dominance_model) %>% spread(key = Variant, value = Dominance_model) %>% dplyr::select(-Date)))
    model_mat[model_mat <= 0] = 0.000000000000001

    # Replace NaN values across all columns
    model_mat[is.na(model_mat)] <- 0.0000000000001

    data_mat = as.matrix(as.data.frame(dplyr::select(data_ll, Date, Variant, Dominance_data) %>% spread(key = Variant, value = Dominance_data) %>% dplyr::select(-Date)))
    data_mat[data_mat <= 0] = 0.000000000000001

    nll = -ddirichlet(data_mat, model_mat, log = T, sum.up = T)
    ## if(is.na(nll)){
    ##     browser()
    ## }
    return(nll)
}

##========================================#
## Fix files -------------
##========================================#
check_finished_jobs <- function(params_in){
    job_id = params_in$job_id[1] ## This should be only the job parameters
    numDays = params_in$days[1]
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

    tmp_out = sprintf('tmp_%s.csv', job_id)
    fred_csv = sprintf("fred_csv -k %s > %s", job_id, tmp_out)    
    fred_csv_st = system(fred_csv, intern = T)
    
    if(!is.null(attr(fred_csv_st, "status"))){
        unlink(tmp_out)
        return(FALSE)
    }
    job_df = read.csv(tmp_out, stringsAsFactors = FALSE)
    unlink(tmp_out)

    if(!(file.exists(file.path(data_dir,'AgeGroupsData_mean.csv')) & nrow(job_df) == numDays)){
        return(FALSE)
    }
    return(TRUE)
}

process_job_data <- function(job_id, process_age = FALSE, age_data = data.frame()){
    #browser()
    fred_dir = system(sprintf("fred_find -k %s", job_id), intern = TRUE)
    data_dir = file.path(fred_dir, "DATA", "OUT")   
    fred_csv = sprintf("fred_csv -k %s", job_id)
    ##system(fred_csv, intern = T)

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
    
    ##========================================#
    ## Read data and compute LL --------------
    ##========================================#
    job_df$LL = Inf
    job_df$Date = as.Date('2020-01-01') + job_df$Day

    ## Group by date and get the overall deaths
    ## Near future -> collect by localidad
    # tmp_data = filter(BOG_data, DeptCode == state_input) %>%
    #     group_by(Date) %>%
    #     summarize(Cases = sum(Cases), Deaths = sum(Deaths), ICU_admissions = sum(ICU_admissions)) %>%
    #     ungroup() %>%
    #     left_join(job_df, by = c("Date" = "Date"))    

    tmp_data = filter(BOG_data, DeptCode == state_input) %>%
        group_by(Date) %>%
        summarize(Deaths = sum(Deaths), UCI_prevalence = sum(UCI_prevalence)) %>%
        ungroup() %>%
        left_join(job_df, by = c("Date" = "Date"))    
    

    tmp_data$CF_total = tmp_data$CF_mean
    if('CF_1_mean' %in% colnames(job_df)){
        tmp_data$CF_total = tmp_data$CF_1_mean + tmp_data$CF_mean
    }
    if('CF_2_mean' %in% colnames(job_df)){
        tmp_data$CF_total = tmp_data$CF_2_mean + tmp_data$CF_total
    }
    if('CF_3_mean' %in% colnames(job_df)){
        tmp_data$CF_total = tmp_data$CF_3_mean + tmp_data$CF_total
    }
    if('CF_4_mean' %in% colnames(job_df)){
        tmp_data$CF_total = tmp_data$CF_4_mean + tmp_data$CF_total
    }
    if('CF_5_mean' %in% colnames(job_df)){
        tmp_data$CF_total = tmp_data$CF_5_mean + tmp_data$CF_total
    }
    if('CF_6_mean' %in% colnames(job_df)){
        tmp_data$CF_total = tmp_data$CF_6_mean + tmp_data$CF_total
    }
    if(!('Chosp_mean' %in% colnames(job_df))){
        job_df$Chosp_mean = NA
    }

    tmp_data$Chosp_total = tmp_data$Chosp_mean
    if('Chosp_1_mean' %in% colnames(job_df)){
        tmp_data$Chosp_total = tmp_data$Chosp_total + tmp_data$Chosp_1_mean
    }
    if('Chosp_2_mean' %in% colnames(job_df)){
        tmp_data$Chosp_total = tmp_data$Chosp_total + tmp_data$Chosp_2_mean
    }
    if('Chosp_3_mean' %in% colnames(job_df)){
        tmp_data$Chosp_total = tmp_data$Chosp_total + tmp_data$Chosp_3_mean
    }
    if('Chosp_4_mean' %in% colnames(job_df)){
        tmp_data$Chosp_total = tmp_data$Chosp_total + tmp_data$Chosp_4_mean
    }
    if('Chosp_5_mean' %in% colnames(job_df)){
        tmp_data$Chosp_total = tmp_data$Chosp_total + tmp_data$Chosp_5_mean
    }
    if('Chosp_6_mean' %in% colnames(job_df)){
        tmp_data$Chosp_total = tmp_data$Chosp_total + tmp_data$Chosp_6_mean
    }


    tmp_data$Phosp_total = tmp_data$Phosp_mean
    if('Phosp_1_mean' %in% colnames(job_df)){
        tmp_data$Phosp_total = tmp_data$Phosp_total + tmp_data$Phosp_1_mean
    }
    if('Phosp_2_mean' %in% colnames(job_df)){
        tmp_data$Phosp_total = tmp_data$Phosp_total + tmp_data$Phosp_2_mean
    }
    if('Phosp_3_mean' %in% colnames(job_df)){
        tmp_data$Phosp_total = tmp_data$Phosp_total + tmp_data$Phosp_3_mean
    }
    if('Phosp_4_mean' %in% colnames(job_df)){
        tmp_data$Phosp_total = tmp_data$Phosp_total + tmp_data$Phosp_4_mean
    }
    if('Phosp_5_mean' %in% colnames(job_df)){
        tmp_data$Phosp_total = tmp_data$Phosp_total + tmp_data$Phosp_5_mean
    }
    if('Phosp_6_mean' %in% colnames(job_df)){
        tmp_data$Phosp_total = tmp_data$Phosp_total + tmp_data$Phosp_6_mean
    }

    tmp_data <- tmp_data %>% mutate(CF_total = replace_na(CF_total, 0),
                                    Phosp_total = replace_na(Phosp_total, 0), 
                                    Chosp_total = replace_na(Chosp_total, 0))

    tmp_data_dic = tmp_data %>% 
        filter(Date < as.Date('2021-02-01'))

    # job_df_daily <- job_df

    # job_df$Date <- as.Date(job_df$Date)
    # tmp_data$Date <- as.Date(tmp_data$Date)
    # age_data$Date <- as.Date(age_data$Date)
    

    # #browser()

    # # Assuming `job_df` requires specific numeric columns, adjust `c()` as needed
    # job_df <- job_df %>%
    # select(Date, N_mean, matches("^C(_\\d+)?_mean$"), starts_with('ACF')) %>%
    # mutate(Date = floor_date(Date, unit = "week")) %>%
    # group_by(Date) %>%
    # summarise(N_mean = first(N_mean),  # Retain the first N_mean value for each group
    #             across(where(is.numeric) & !matches("^N_mean$"), sum, na.rm = TRUE),
    #             .groups = 'drop')

    # # Assuming `tmp_data` requires specific numeric columns
    # tmp_data <- tmp_data %>%
    # select(Date, Deaths, ICU_admissions, CF_total, Chosp_total) %>%
    # mutate(Date = floor_date(Date, unit = "week")) %>%
    # group_by(Date) %>%
    # summarise(across(where(is.numeric), sum, na.rm = TRUE), .groups = 'drop')

    # # For `age_data`, assuming the inclusion of 'Deaths' with 'AgeGroup'
    # age_data <- age_data %>%
    # select(Date, AgeGroup, Deaths) %>%
    # mutate(Date = floor_date(Date, unit = "week")) %>%
    # group_by(Date, AgeGroup) %>%
    # summarise(across(where(is.numeric), sum, na.rm = TRUE), .groups = 'drop')

    #browser()
    if(nrow(job_df) > 0){
        tmp_LL = get_loglikelihood(tmp_data$Deaths, tmp_data$CF_total)
        tmp_LL_icu = get_loglikelihood(tmp_data$UCI_prevalence, tmp_data$Phosp_total)
        
        tmp_LL_dic = get_loglikelihood(tmp_data_dic$Deaths, tmp_data_dic$CF_total)
        tmp_LL_icu_dic = get_loglikelihood(tmp_data_dic$UCI_prevalence, tmp_data_dic$Phosp_total)

        job_age_df = job_df %>%
            dplyr::select(Date, N_mean,  starts_with('ACF')) %>%
            gather(key = AgeGroup, value = CF_mean, -c('Date','N_mean')) %>%
            group_by(Date, N_mean,  AgeGroup) %>%
            summarize(CF_mean = sum(CF_mean, na.rm = T),.groups = 'drop') %>%
            left_join(BOG_age_time_data, by = c('Date','AgeGroup')) %>%
            drop_na()
        tmp_job_df = job_df %>%
            dplyr::select(Date, N_mean, C_mean) %>%
            filter(Date <= as.Date('2020-11-07')) %>%
            summarize(AR_mean = sum(C_mean, na.rm = T),.groups = 'drop',
                      AR_data = as.integer(0.25 * N_mean[1])) %>%
            drop_na()
        tmp_LL_age_time = get_agegroups_likelihood(job_age_df)
        
        tmp_LL_AR = get_loglikelihood_poisson(data_in = tmp_job_df$AR_data, model_in = tmp_job_df$AR_mean)
        
        tmp_LL_df = data.frame(LL = tmp_LL, stringsAsFactors = FALSE)
        if(!is.na(tmp_LL)){
            job_df$LL = tmp_LL
            job_df$LL_dic = tmp_LL_dic
            job_df$LL_age_time = tmp_LL_age_time
            tmp_LL_df$LL_age_time = tmp_LL_age_time
            job_df$LL_AR = tmp_LL_AR
            tmp_LL_df$LL_AR = tmp_LL_AR
        }

        if('CF_6_mean' %in% colnames(job_df)){                
            tmp_job = filter(job_df, Date %in% dominance_data$Date) %>%
                mutate(DomAlpha = C_1_mean / (C_1_mean + C_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean + C_6_mean),
                       DomGamma = C_2_mean / (C_1_mean + C_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean + C_6_mean), 
                       DomKappa = C_3_mean / (C_1_mean + C_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean + C_6_mean), 
                       DomDelta = C_4_mean / (C_1_mean + C_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean + C_6_mean),
                       DomOmicron = C_5_mean / (C_1_mean + C_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean + C_6_mean),
                       DomOmicronBAX = C_6_mean / (C_1_mean + C_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean + C_6_mean)) %>%
                dplyr::select(Date, starts_with('Dom')) %>%
                gather(key = Variant, value = Dominance, -Date) %>%
                mutate(Variant = str_replace(Variant,'Dom',''))

            tmp_LL_dom =  get_loglikelihood_dominance(dominance_data, tmp_job)
            tmp_LL_df$LL_dom = tmp_LL_dom} else if('CF_1_mean' %in% colnames(job_df)){                
            tmp_job = filter(job_df, Date %in% dominance_data$Date) %>%
                mutate(DomAlpha = C_1_mean / (C_1_mean + C_mean + C_2_mean + C_3_mean + C_4_mean),
                       DomGamma = C_2_mean / (C_1_mean + C_mean + C_2_mean + C_3_mean + C_4_mean), 
                       DomKappa = C_3_mean / (C_1_mean + C_mean + C_2_mean + C_3_mean + C_4_mean), 
                       DomDelta = C_4_mean / (C_1_mean + C_mean + C_2_mean + C_3_mean + C_4_mean)) %>%
                dplyr::select(Date, starts_with('Dom')) %>%
                gather(key = Variant, value = Dominance, -Date) %>%
                mutate(Variant = str_replace(Variant,'Dom',''))

            tmp_LL_dom =  get_loglikelihood_dominance(dominance_data, tmp_job)
            tmp_LL_df$LL_dom = tmp_LL_dom
        }

        if(nrow(age_data) > 0){
            ## This has to be changed to adjust Bogota data
            
            age_data$TotalDeaths  = sum(age_data$Deaths)
            age_data$NDeaths = age_data$Deaths
            age_data$PDeaths = age_data$NDeaths / age_data$TotalDeaths
            
            ## AGE STRUCTURE DEATHS
            tmp_fred_CF =  job_df %>% filter(Date <= age_data$Date[1])%>%
                summarize_at(vars(starts_with('ACF')), 'sum', rm.na = T) %>%
                gather(key = AgeGroup, value = CF_mean) %>%
                left_join(age_data, by = "AgeGroup")

            LL_CF_age = get_death_age_likelihood(replace_na(tmp_fred_CF$NDeaths, 0), replace_na(tmp_fred_CF$CF_mean, 0))
            
            total_LL = tmp_LL + LL_CF_age + tmp_LL_icu
            job_df$LL = total_LL
            total_LL_dic = tmp_LL_dic + LL_CF_age + tmp_LL_icu_dic
            job_df$LL_dic = total_LL_dic
            tmp_LL_df$LL_deaths = tmp_LL
            tmp_LL_df$LL_icu = tmp_LL_icu
            tmp_LL_df$LL_deaths_dic = tmp_LL_dic
            tmp_LL_df$LL_icu_dic = tmp_LL_icu_dic
            tmp_LL_df$LL = total_LL
            tmp_LL_df$LL_CF_age = LL_CF_age
            tmp_LL_df$LL_dic = total_LL_dic
        }        
        
    }
    if('LL_dom' %in% colnames(tmp_LL_df)){
        print(sprintf("LL: %.2f LLdic: %.2f LLdom: %.2f", job_df$LL[1], job_df$LL_dic[1], tmp_LL_df$LL_dom))
    }else{
        print(sprintf("LL: %.2f LLdic: %.2f", job_df$LL[1], job_df$LL_dic[1]))
    }
    if(is.infinite(job_df$LL[1])){browser()}

    ##browser()

    cols_output = c('Day', 'Week', 'Year', 'Runs', 'N_mean', 'N_std',
                    'C_mean', 'C_std', 'Cs_mean', 'Cs_std','CFR_mean', 'Is_mean', 'I_mean',
                    'CFR_std', 'CF_mean', 'CF_std',
                    'RR_mean', 'RR_std', 'AR_mean', 'AR_std', 'Nursing_Home_mean', 
                    'Nursing_Home_CF_mean',
                    'H_sheltering_mean', 'N_sheltering_mean','LL', 'Chosp_mean')

    cols_output = c(cols_output, cols_age)
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

    # Add columns starting with "VA" to cols_output
    vax_cols <- grep("^V", colnames(job_df), value = TRUE)
    cols_output <- c(cols_output, vax_cols)

    #browser()
    #extra_LL_cols <- c("LL", "LL_dic", "LL_age_time", "LL_AR", "V_mean", "Vtd_mean")
    #job_df_daily <- job_df_daily %>% dplyr::select(-all_of(extra_LL_cols))

    # job_df_daily$LL <- tmp_LL_df$LL[[1]]
    # job_df_daily$LL_dic <- tmp_LL_df$LL_dic[[1]]
    # job_df_daily$LL_age_time <- tmp_LL_df$LL_age_time[[1]]
    # job_df_daily$LL_AR <- tmp_LL_df$LL_AR[[1]]
    # job_df_daily$V_mean <- job_df$V_mean[[1]]
    # job_df_daily$Vtd_mean <- job_df$Vtd_mean[[1]]

    print(sprintf("Sucess processing job %s", job_id))

    return(list(job_df = job_df[,cols_output],
                extra_params_df = tmp_LL_df))
}

fred_gather_data <- function(params, outdir, outfile, rm_out = FALSE, appendToFile = FALSE, process_age = FALSE, age_data = data.frame()){
  # Check and prepare the output directory
  if (dir.exists(outdir) && rm_out) {
    unlink(outdir, recursive = TRUE)
    dir.create(outdir)
  } else if (!dir.exists(outdir)) {    
    dir.create(outdir)
  }

  # Initialize the output files
  outfile_path <- file.path(outdir, outfile)
  if (file.exists(outfile_path)) {
    unlink(outfile_path)
  }
  
  # Copy and read parameter file
  param_file_path <- file.path(outdir, basename(params))
  file.copy(params, param_file_path)
  params_orig <- read_csv(params, show_col_types = FALSE) %>% mutate(Finished = 0)
  
  # Set up parallel cluster
  cl <- makeCluster(detectCores() - 1)
  on.exit(stopCluster(cl))  # Ensures the cluster is stopped when the function exits

  # Export necessary objects to each worker node
  env_vars <- list(params_orig = params_orig, process_age = process_age, age_data = age_data)
  clusterExport(cl=cl, varlist=c("env_vars", "BOG_data", "dominance_data", "ddirichlet", "state_input", "BOG_age_time_data", "get_loglikelihood_dominance", "get_death_age_likelihood",  "get_loglikelihood_poisson", "get_loglikelihood", "get_agegroups_likelihood"), envir=environment())
  clusterExport(cl, c("check_finished_jobs", "process_job_data"))


  # Define progress tracking functions and initialize progress tracking
  clusterEvalQ(cl, {
    library(dplyr)
    library(readr)
    library(tidyverse)

    incrementClusterProgress <- function() {
      progress$counter <- progress$counter + 1
      cat(sprintf("\rProgress: %d/%d", progress$counter, progress$total))
      flush.console()
    }
  })
  progress <- new.env()
  progress$counter <- 0
  progress$total <- nrow(params_orig)
  clusterExport(cl, varlist=c("progress"), envir=environment())

  #browser()
    # params_row <- params_orig[1, , drop = FALSE]
    # job_processed_list <- process_job_data(params_row$job_id, process_age, age_data)
    
    # Parallel processing  seq_len(nrow(params_orig))
    results <- parLapply(cl, seq_len(nrow(params_orig)), function(index) {
        params_row <- params_orig[index, , drop = FALSE]
        if (check_finished_jobs(params_row)) {
            job_processed_list <- process_job_data(params_row$job_id, process_age, age_data)

            # Ensure job_id is present and non-empty; provide a fallback if necessary
            job_id_to_use <- if (!is.na(params_row$job_id) && length(params_row$job_id) > 0) {
                as.character(params_row$job_id)
            } else {
                "default_job_id"  # Fallback job_id
            }

            job_processed <- job_processed_list$job_df %>%
                mutate(job_id = job_id_to_use, Finished = 1)

            if ("extra_params_df" %in% names(job_processed_list)) {
                params_row <- bind_cols(params_row, job_processed_list$extra_params_df)
            }

            return(list(params_tmp = params_row, job_processed = job_processed))
        } else {
            return(list(params_tmp = params_row, job_processed = tibble(job_id = character(), Finished = integer())))
        }
    })

 # Combine results using bind_rows instead of rbind
  #browser()
  params_finished <- bind_rows(lapply(results, '[[', "params_tmp"))
  job_processed_all <- bind_rows(lapply(results, '[[', "job_processed"))

  # Write results to files
  params_outfile <- file.path(outdir, 'FRED_parameters_out.csv')
  write_csv(x = params_finished, path = params_outfile)

  if (appendToFile) {
    write_csv(job_processed_all, path = outfile_path, append = TRUE)
  } else {
    job_processed_all <- filter(job_processed_all, Finished == 1)
    write_csv(job_processed_all, path = outfile_path)
  }

  return(0)
}


##========================================#
## Inputs --------------------
##========================================#
calibration_label = 'production'
state_input = 27
data_days_rm = 4
variants_in = 0
args = (commandArgs(TRUE))

process_age_flag = TRUE
asymp_infectivity_in = 1.0
face_mask_transmission_efficacy_in = 0.73
kids_susceptibility_age_in = 10
vaccination_in = 70
if(length(args) >= 1){
    state_input = as.numeric(args[1])
    if(length(args) >=2){
        data_days_rm = as.numeric(args[2])
        if(length(args) >= 3){
            asymp_infectivity_in = as.numeric(args[3])
            if(length(args) >= 4){
                face_mask_transmission_efficacy_in = as.numeric(args[4])
                if(length(args) >= 5){
                    kids_susceptibility_age_in = as.numeric(args[5])
                    if(length(args) >= 6){
                        variants_in = as.numeric(args[6])
                    }
                }    
            }
        }
    }
}

##========================================#
## Collect for covid --------------------
##========================================#
## Incidence data 
incidence_file      = 'COL_covid_death_data.csv'
incidence_file_sm   = 'smooth_test_27.csv'
age_incidence_file  = 'Age_COL_covid_data.csv'
uci_incidence_file  = 'COL_UCI_prevalence_timeseries.csv'
variance_dom_file   = sprintf('%s_dominance_data.csv', state_input)

if(!file.exists(file.path('../../fred_input_files/covid_data',incidence_file))){
    stop("Incidence data not found")
}

interventions_df = read_csv('../../fred_input_files/interventions/interventions_Colombia.csv') %>% dplyr::filter(State == state_input)
interventions_df$Shelter_in_place <- as.Date(interventions_df$Shelter_in_place, format="%m/%d/%y")
interventions_df$School_closure <- as.Date(interventions_df$School_closure, "%m/%d/%y")

BOG_uci_df = read_csv(file.path('../../fred_input_files/covid_data',uci_incidence_file)) %>%
    dplyr::select(Date, UCI_prevalence) %>%
    filter(Date < (max(Date, na.rm = T) - 2))

# BOG_data = read_csv(file.path('../../fred_input_files/covid_data', incidence_file)) %>%
#     mutate(DeptCode = as.numeric(substr(as.character(MunCode), 1, nchar(MunCode) - 3))) %>%
#     filter(DeptCode %in% interventions_df$State) %>%
#     filter(Date < (max(Date, na.rm = T) - data_days_rm)) %>%
#     left_join(interventions_df, by = c("MunCode" = "State")) %>%
#     right_join(BOG_uci_df, by = 'Date')

BOG_data <- read_csv(file.path('../../fred_input_files/covid_data', incidence_file_sm)) %>% 
            dplyr::select(Date, smoothed_by_parts_integer) %>%
            rename(c("Deaths" = "smoothed_by_parts_integer")) %>%
            mutate(DeptCode = state_input) %>%
            left_join(BOG_uci_df, by = 'Date') %>% 
            replace_na(list(UCI_prevalence = 0))


BOG_age_time_data = read_csv(file.path('../../fred_input_files/covid_data', age_incidence_file)) %>%
    mutate(DeptCode = as.numeric(substr(as.character(MunCode), 1, nchar(MunCode) - 3))) %>%
    mutate(MunCode = as.numeric(MunCode)) %>%
    filter(DeptCode %in% interventions_df$State) %>%
    filter(Date < (max(Date, na.rm = T) - data_days_rm)) %>%
    mutate(MaxDate = max(Date, na.rm = T)) %>%
    left_join(interventions_df, by = c("MunCode" = "State")) %>%  
    group_by(MunCode, Date, AgeGroup) %>%
    summarize(Deaths = sum(Deaths, na.rm = T)) %>%
    ungroup() 

# # For `BOG_age_time_data`, assuming the inclusion of 'Deaths' and other numeric cols
# BOG_age_time_data$Date <- as.Date(BOG_age_time_data$Date)
# BOG_age_time_data <- BOG_age_time_data %>%
#                     select(MunCode, Date, AgeGroup, starts_with('Death')) %>%
#                     mutate(Date = floor_date(Date, unit = "week")) %>%
#                     group_by(MunCode, Date, AgeGroup) %>%
#                     summarise(across(where(is.numeric), sum, na.rm = TRUE), .groups = 'drop')

BOG_age_data =  read_csv(file.path('../../fred_input_files/covid_data', age_incidence_file)) %>%
    mutate(DeptCode = as.numeric(substr(as.character(MunCode), 1, nchar(MunCode) - 3))) %>%
    mutate(MunCode = as.numeric(MunCode)) %>%
    filter(DeptCode %in% interventions_df$State) %>%
    filter(Date < (max(Date, na.rm = T) - data_days_rm)) %>%
    mutate(MaxDate = max(Date, na.rm = T)) %>%
    left_join(interventions_df, by = c("MunCode" = "State")) %>%
    group_by(AgeGroup) %>%
    summarize(Deaths = sum(Deaths), Date = MaxDate[1]) %>%
    ungroup()

Bog_dom_data = read_csv(file.path('../../fred_input_files/dominance_data', variance_dom_file)) %>%
    rename(Date = date, Kappa = Mu)

names(Bog_dom_data)

col_variants = 1:5
variants_list = c('Alpha','Gamma','Kappa','Delta', 'Omicron')
Bog_dom_data <- Bog_dom_data %>% dplyr::select('Date','days','Alpha','Gamma','Kappa','Delta', 'Omicron')

dominance_data = Bog_dom_data %>% gather(key = Variant, value = Dominance, -Date, -days) %>%
    left_join(data.frame(Color = col_variants[1:length(variants_list)], Variant = variants_list), by = 'Variant') %>%
    drop_na()

## TODO:
## [ ] USE INS data instead of Bogotá data -> More data    
simdir_name = sprintf('FRED_%.d_calibration_asymp_%.2f_fm_%.2f_ksus_%.2f_var_%.0f_vax_%03d_mov_%s',state_input, asymp_infectivity_in,face_mask_transmission_efficacy_in,kids_susceptibility_age_in, variants_in, vaccination_in, calibration_label)
simdir = file.path(Sys.getenv('scratch_dir'), simdir_name)

if(!dir.exists(file.path(getwd(), '../output','CALIBRATION'))){
    dir.create(file.path(getwd(), '../output','CALIBRATION'), recursive = T)
}

fred_results_dir = file.path(simdir,"FRED_RESULTS")
Sys.setenv(FRED_RESULTS=fred_results_dir)


job_id_list = list.dirs(file.path(fred_results_dir, "JOB"), recursive=FALSE, full.names=FALSE)
system(paste('rm -rf ', file.path(fred_results_dir, "KEY")))
for(job_id in job_id_list){
    cat(sprintf("FRED_%s_calibration_%s %s", state_input, job_id, job_id), file=file.path(fred_results_dir, "KEY"), sep="\n", append = TRUE)
}

output_dir = file.path(getwd(), '../output','CALIBRATION', sprintf("%s_%s", simdir_name, "out"))
fred_outputs = 'fred_output.csv'
params_file = sprintf("%s/%s",simdir, "FRED_parameters.csv")
params_limit_file = sprintf("%s/%s",simdir, "FRED_parameters_limits.csv")


fred_gather_data(params=params_file, outdir=output_dir, outfile=fred_outputs,
                            rm_out = TRUE,
                            appendToFile = FALSE,
                            process_age = process_age_flag,
                            age_data = BOG_age_data)



if(!file.exists(file.path(params_limit_file))){
    stop("Limits data not found")
}

print(file.path(output_dir,'FRED_parameters_limits.csv'))
file.copy(params_limit_file, file.path(output_dir,'FRED_parameters_limits.csv'), overwrite = T)


## Save a copy of the file used for calibration: don't save today's or yesterday's data
epi_data_out = read_csv(file.path('../../fred_input_files/covid_data', incidence_file)) %>%
    mutate(DeptCode = as.numeric(substr(as.character(MunCode), 1, nchar(MunCode) - 3))) %>%
    mutate(MunCode = as.numeric(MunCode)) %>%
    filter(DeptCode %in% interventions_df$State) %>%
    filter(Date < (max(Date, na.rm = T) - data_days_rm)) 

write.csv(epi_data_out, file.path(output_dir, incidence_file))
write.csv(BOG_age_data, file.path(output_dir, age_incidence_file))