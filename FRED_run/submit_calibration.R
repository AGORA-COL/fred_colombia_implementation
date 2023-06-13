##==============================================#
## Author: Guido Espana
## Simulate COVID-19 in BOGOTA
## Year: 2019
##==============================================#
## User's input----------------
##==============================================#
library(pomp)
library(lubridate)
library(tidyverse)
library(fredtools)
library(doParallel)

parallel::detectCores()
n.cores <- parallel::detectCores() - 1

my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "FORK"
  )

print(my.cluster)
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

####################################
## FRED enviromental variables
####################################

Sys.setenv(FRED_HOME='/shared/home/Azure-ASIS/FRED')
Sys.setenv(FRED_RESULTS='/shared/home/Azure-ASIS/FRED_Implementation/FRED_results')
Sys.setenv(PATH='/bin/:/shared/home/Azure-ASIS/anaconda3/condabin:/usr/local/cuda/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/usr/local/bin:/usr/local/sbin:/shared/home/Azure-ASIS/.local/bin:/shared/home/Azure-ASIS/bin:/shared/home/Azure-ASIS/FRED/bin:/shared/home/Azure-ASIS/.local/bin:/shared/home/Azure-ASIS/bin')
Sys.setenv(scratch_dir='/shared/home/Azure-ASIS/FRED_Implementation/scratch')

fred_home = Sys.getenv('FRED_HOME')
fred_defaults = sprintf("%s/input_files/defaults", fred_home)

##==============================================#
## Functions severity-------------------
##==============================================#
get_fred_dis_prog_params <- function(symp.df, ifr.df, dis_name){
    ##======================================#
    ## Recovery from symptoms ------------
    ##======================================#
    ## Recovery (I think this comes from Bi et al. Need to check 
    ## This first bit calculates the probabilty you have recovered by a certain day
    median = 20.8
    CI95 = 0.7
    n = 225 

    variance = n*(CI95/1.96)^2
    meanlog = log(median)
    sdlog = sqrt(log((1+sqrt(1+4*variance/exp(2*meanlog)))/2))

    recover.cdf = plnorm(seq(0,100,1),meanlog,sdlog)
    trunc.val = which(recover.cdf > 0.999)[1]
    rec.cdf = recover.cdf[1:trunc.val]/recover.cdf[trunc.val]

    ##======================================#
    ## isolation rate ------------
    ##======================================#
    ## Needs to be after recovery
    isolation.func <- function(p) {
        isolation.vec = c(0,(1-p)^seq(0,length(recover.cdf)-2,1))
        isolate.p = 1 - sum(isolation.vec*c(0,diff(recover.cdf)))
        return(abs(isolate.p - isolation.prob))
    }

    isolation.prob <- 0.68
    ## Using MMWR paper
    isolation.prob.daily <- optimize(
        interval = c(0,1),maximum = F,    
        f= isolation.func)$minimum

    ## Print daily isolation prob:
    isolation.prob.daily

    ## Check it matches the literature value
    isolation.vec = (1-isolation.prob.daily)^seq(0,100,1)
    1 - sum(isolation.vec*c(0,diff(recover.cdf)))

    isolation.fred.rate = isolation.prob.daily


    ##======================================#
    ## IFR ------------
    ##======================================#
    ##Death timing ##imperial lancet global health
    symp_to_death_mean = 17.8
    symp_to_death_mode = 14.5
    symp_to_death_rate = 1/(symp_to_death_mean-symp_to_death_mode)
    symp_to_death_shape = symp_to_death_mean*symp_to_death_rate
    rows <- findInterval(ifr.df$ages,symp.df$ages,left.open=TRUE)+1

    ifr.df$symp.probs <- symp.df$probs[rows]
    ifr.df$cfr.probs <- 0.01*ifr.df$probs/ifr.df$symp.probs

    cfr.vec = ifr.df$cfr.probs
    cfr = cfr.vec[length(cfr.vec)]
    days = 1:length(rec.cdf)-1
    rec.pdf = diff(rec.cdf)

    death.cdf = pgamma(days,symp_to_death_shape,symp_to_death_rate)
    death.cdf = cfr*death.cdf/death.cdf[length(death.cdf)]
    death.pdf = diff(death.cdf)

    ## calculate daily death probability
    fred.death.pdf = rep(0,length(days)-1)
    for (t in 1:(length(days)-1)){
        fred.death.pdf[t] = (death.pdf[t])/sum(rec.pdf[t:length(rec.pdf)])/prod(1-fred.death.pdf[0:(t-1)])
    }

    fred.death.pdf=pmax(pmin(fred.death.pdf,1),0)

    test.death.pdf = fred.death.pdf*cumprod(1-c(0,fred.death.pdf[-length(fred.death.pdf)]))*rev(cumsum(rev(rec.pdf)))

    ## check
    abs(death.pdf - test.death.pdf) < 1e-10

    ## get cfr multipliers
    cfr.multipliers <- cfr.vec/cfr
    for (i in 1:length(cfr.vec)) {
        cfr.mult.temp = cfr.multipliers[i]
        cfr.multipliers[i] = optim(par=cfr.mult.temp,
                                   fn=function(par){
                                       test.death.pdf = par*fred.death.pdf*cumprod(1-c(0,par*fred.death.pdf[-length(fred.death.pdf)]))*rev(cumsum(rev(rec.pdf)))
                                       return(abs(cfr.vec[i]-sum(test.death.pdf)))
                                   }, method="Brent",lower=0,upper=1
                                   )$par
    }

    ##======================================#
    ## Print to FRED format ------------
    ##======================================#
    days_symptomatic_str = sprintf("%s_days_symptomatic = %s",dis_name, paste(c(trunc.val,sprintf("%.9f",rec.cdf)), collapse = " "))

    case_fatality_values_str = sprintf("%s_case_fatality_values = %s", dis_name, paste(c(length(cfr.multipliers),sprintf("%.9f",cfr.multipliers)), collapse = " "))

    case_fatality_age_groups_str = sprintf("%s_case_fatality_age_groups = %s", dis_name, paste(c(nrow(ifr.df),sprintf("%.0f",ifr.df$ages)), collapse = " "))

    prob_symptoms_values_str = sprintf("%s_prob_symptoms_values = %s", dis_name, paste(c(nrow(symp.df),sprintf("%.9f",symp.df$probs)), collapse = " "))

    prob_symptoms_age_groups_str = sprintf("%s_prob_symptoms_age_groups = %s", dis_name,paste(c(nrow(symp.df),sprintf("%.0f",symp.df$ages)), collapse = " "))

    fred.death.pdf[length(fred.death.pdf)] = 1.0
    max.death.day = which(fred.death.pdf > 0.9)[1]
    if(max.death.day < length(fred.death.pdf)){
        ind_post_max = max.death.day:length(fred.death.pdf)    
        fred.death.pdf = fred.death.pdf[-ind_post_max[fred.death.pdf[ind_post_max] == 0]]
    }

    ## print the death pdf by day
    case_fatality_prob_by_day_str = sprintf("%s_case_fatality_prob_by_day = %s", dis_name, paste(c(length(fred.death.pdf),sprintf("%.9f",fred.death.pdf)), collapse = " "))

    fred_output_lines = c(
        days_symptomatic_str,
        prob_symptoms_age_groups_str,
        prob_symptoms_values_str,
        case_fatality_age_groups_str,        
        case_fatality_values_str,
        case_fatality_prob_by_day_str)
    
    return(fred_output_lines)
}

get_ifr_param_string <- function(severity_factor, var_name){
    davies.symp.df <- data.frame(
        ages=c(seq(10,70,10),120),
        probs=c(0.29, 0.21,0.27,0.33,0.4,0.49,0.63,0.69)) 
    verity.ifr.df <- data.frame(
        ages=c(9,19,29,39,49,59,69,79,120),
        probs=c(0.00161,0.00695,0.0309,0.0844,0.161,0.595,1.93,4.28,7.8))
    combined_params_df = data.frame(stringsAsFactors = FALSE)
    for(vv in 1:length(severity_factor)){
        variant.ifr.df = verity.ifr.df
        variant.ifr.df$probs =  100 * (1 - exp(-log(1/(1- 0.01*variant.ifr.df$probs)) * severity_factor[vv]))

        
        variant_params = get_fred_dis_prog_params(symp.df = davies.symp.df,
                                                  ifr.df = variant.ifr.df, dis_name = sprintf("variant%s",var_name))

        variant_params_df = lapply(str_split(variant_params, "[:space:]*=[:space:]*"), function(x){x[2]})
        names(variant_params_df) = lapply(str_split(variant_params, "[:space:]*=[:space:]*"), function(x){x[1]})
        combined_params_df = bind_rows(combined_params_df, as.data.frame(variant_params_df))
    }
    return(combined_params_df)
}
##==============================================#
## Custom functions-------------
##==============================================#
write_cmd_function <- function(scalars_in, tmpfile){
    ## Delte jobs!!!
    cat("Deleting jobs\n")
    if(dir.exists(file.path(Sys.getenv('FRED_RESULTS'),'JOB'))){
        for(s in 1:nrow(scalars_in)){
            system(sprintf("fred_delete -f -k %s",scalars_in$job_id[s]), intern = T)
        }
        cat("All jobs deleted\n")
    }else{
        cat("FRED RESULTS doesn't exist, nothing to delete",Sys.getenv('FRED_RESULTS'),"\n")
    }

    job_cmd_str = sprintf("fred_job -k %s -I %d -p %s -n %.0f; Rscript ./post_process_fred_calibration.R %s %.0f",
                          scalars_in$job_id,
                          scalars_in$run_id,
                          scalars_in$params_file,
                          scalars_in$reps,
                          scalars_in$job_id,
                          scalars_in$reps
                          )

    fileConn<-file(tmpfile)
    writeLines(job_cmd_str, fileConn)
    close(fileConn)
}

##==============================================#
## Set STATE-------------------
##==============================================#
state_code = 11001
reps = 2000
reps_per_job = 1
fit_date = as.Date('2022-01-01')
asymp_infectivity_in = 1.0
face_mask_transmission_efficacy_in = 0.73
kids_susceptibility_age_in = 10
variants_in = 1
vaccination_in = 70
args = (commandArgs(TRUE))
if(length(args) >= 1){
    state_code = as.numeric(args[1])
    if(length(args) >= 2){
        reps = as.numeric(args[2])
        if(length(args) >= 3){
            reps_per_job = as.numeric(args[3])
            if(length(args) >= 4){
                fit_date = as.Date(args[4])
                if(length(args) >= 5){
                    asymp_infectivity_in = as.numeric(args[5])
                    if(length(args) >= 6){
                        face_mask_transmission_efficacy_in = as.numeric(args[6])
                        if(length(args) >= 7){
                            kids_susceptibility_age_in = as.numeric(args[7])
                            if(length(args) >= 8){
                                variants_in = as.numeric(args[8])
                                if(length(args) >= 9){
                                    vaccination_in = as.numeric(args[9])
                                }
                            }
                        }                        
                    }
                }
            }
        }
    }
}

##==============================================#
## FRED setup-------------
##==============================================#
fred_home = Sys.getenv('FRED_HOME')
fred_defaults = sprintf("%s/input_files/defaults", fred_home)
output.dir = file.path(Sys.getenv('scratch_dir'),sprintf('FRED_%d_calibration_asymp_%.2f_fm_%.2f_ksus_%.2f_var_%.0f_vax_%03d_mov_old_files_v2',state_code, asymp_infectivity_in,face_mask_transmission_efficacy_in,kids_susceptibility_age_in, variants_in,vaccination_in))
print(output.dir)
if(file.exists(output.dir)){
    system(paste('rm -rf ', output.dir,sep = ''))
}
system(paste('mkdir -p ', output.dir,sep = ''))

file.copy('../scripts/post_process_fred_calibration.R',output.dir)
file.copy("./input_files/infection_hospitalization_risk.csv", output.dir)
# file.copy('../input_files/params_covid.txt','./input_files/params_covid.txt', overwrite = T)
# file.copy('../input_files/params_covid_variants.txt','./input_files/params_covid_variants.txt', overwrite = T)
# file.copy('../input_files/params_covid_alpha.txt','./input_files/params_covid_alpha.txt', overwrite = T)
# file.copy('../input_files/params_covid_gamma.txt','./input_files/params_covid_gamma.txt', overwrite = T)
# file.copy('../input_files/params_covid_kappa.txt','./input_files/params_covid_kappa.txt', overwrite = T)
# file.copy('../input_files/params_covid_delta.txt','./input_files/params_covid_delta.txt', overwrite = T)
# file.copy('../input_files/params_covid_vaccine.txt','./input_files/params_covid_vaccine.txt', overwrite = T)
# if(vaccination_in > 0){
#     file.copy(sprintf('../input_files/params_covid_vaccine_%d.txt', vaccination_in),sprintf('./input_files/params_covid_vaccine_%d.txt',vaccination_in), overwrite = T)
# }
# file.copy('../input_files/COL_hosp_duration_time.txt','./input_files/COL_hosp_duration_time.txt', overwrite = T)
# file.copy('../input_files/COL_covid_death_data.csv','./input_files/COL_covid_death_data.csv', overwrite = T)
# file.copy('../input_files/Age_COL_covid_data.csv','./input_files/Age_COL_covid_data.csv', overwrite = T)
# file.copy('../input_files/COL_variant_imports.csv','./input_files/COL_variant_imports.csv', overwrite = T)
# file.copy('../input_files/BOG_covid_death_data.csv','./input_files/BOG_covid_death_data.csv', overwrite = T)
# file.copy('../input_files/BOG_UCI_timeseries.csv','./input_files/BOG_UCI_timeseries.csv', overwrite = T)
# file.copy('../input_files/Age_BOG_covid_data.csv','./input_files/Age_BOG_covid_data.csv', overwrite = T)
# ##file.copy('../input_files/11001_mobility_trends.csv','./input_files/11001_mobility_trends.csv', overwrite = T)
# file.copy('../input_files/facemask_timeseries_compliance.csv','./input_files/facemask_timeseries_compliance.csv', overwrite = T)
# file.copy('../input_files/Localidad_Unidad_Catastral.csv','./input_files/Localidad_Unidad_Catastral.csv', overwrite = T)
file.copy('./input_files/Localidad_Unidad_Catastral.csv',output.dir, overwrite = T)
# file.copy('../input_files/Pars-variants.xlsx','./input_files/Pars-variants.xlsx', overwrite = T)

# file.copy(list.files(path='../input_files',pattern='11001_schools_open_gps.*.csv',full.names=T),'./input_files/', overwrite = T)

import_files = list.files('../input_files/','*imports*.csv', full.names = T)
file.copy(import_files,'./input_files/', overwrite = T)

##==============================================#
## Sweep fixed parameters-------------------
##==============================================#
start_date = '2020-01-01' # YYYY-MM-DD
numDays = as.numeric(fit_date - as.Date(start_date)) # We don't need many days yet
track_infection_events_in = 4 # No need to track infections in the calibrations right now
track_fatality_events_in = 1 
report_age_of_infection_in = 4 # Track individual ages of infections
report_incidence_by_county_in = 0 # We don't need this for the calibration either
report_place_of_infection_in = 1 # We might want to track place of infection
synthetic_population_id = sprintf('synthetic_population_id = colombia_%d', state_code)

isolation_rate = 0.05175 # Specific for the US
shelter_in_place_students = 0 # Don't shelter students, use school closure for that

advance_seeding = 'exposed'
school_vacation_end = as.Date('2021-01-25')
current_open_date = as.Date('2020-10-15')

## Shelter in place
interventions_st_df = read_csv('./input_files/interventions_Colombia.csv')

enable_shelter_in_place = 1
enable_shelter_in_place_timeseries = 1
shelter_in_place_delay_mean = as.integer(interventions_st_df$Shelter_in_place[interventions_st_df$State == state_code] - as.Date(start_date))

## Facemask usage
enable_face_mask_timeseries_in = 1
enable_face_mask_usage_in = 1
min_age_face_masks_in = 8

## Age specific susceptibility
susceptibility_params = read_csv('./input_files/age_susceptibility_fit.csv')
enable_age_specific_susceptibility_in = 1

## School closure policies
## School closure
early_school_vacation_end = as.Date('2021-01-25')
early_closure_duration = as.integer(early_school_vacation_end - interventions_st_df$School_closure[interventions_st_df$State == state_code])
early_school_reopening_day = as.numeric(early_school_vacation_end - as.Date(start_date))

school_closure_policy_in = 'global_schedule'
school_closure_day_in = as.integer(interventions_st_df$School_closure[interventions_st_df$State == state_code] - as.Date(start_date))
school_closure_duration_in = as.integer(school_vacation_end - interventions_st_df$School_closure[interventions_st_df$State == state_code])

## School parameters
## Reduced capacity
enable_school_reduced_capacity_in = 0
school_reduced_capacity_in = 0.35
school_reduced_capacity_day_in =  as.integer(early_school_vacation_end - as.Date(start_date))
school_student_teacher_ratio_in = 18

## Nursing home importations
enable_nursing_homes_importations_in = 1

## Contacts
neighborhood_contacts_in = 0.7883
workplace_contacts_in = 0.0686
office_contacts_in = 0.1372
household_contacts_in = 0.1402

## School contacts
school_contacts_in = 0.6295
classroom_contacts_in = 1.2590

## Holiday contacts
enable_community_contact_timeseries_in = 1
enable_holiday_contacts_in = 0
holiday_start_in = "11-16"
holiday_end_in = "01-07"

## Age bias
enable_transmission_bias_in = 1
neighborhood_same_age_bias_in = 0.1

##==============================================#
## PARAMETER SWEEP--------------
##==============================================#
## Parameters to calibrate!!!!
imports_factor = c(0.08,1.0)
influenza_transmissibility = c(0.3,1.5)
shelter_in_place_compliance = c(0.3,1.0)
facemask_compliance = c(0.1,0.96)
workplace_contact_factor_in = c(1.0,1.2)
neighborhood_contact_factor_in = c(0.8,1.2)
household_contact_factor_in = c(0.8,1.5)
community_contact_rate_1_in = c(0.5,2)
influenza_susceptibility_by_age_offset_in = susceptibility_params %>% filter(estimate != "mean") %>% pull(offset_low)
influenza_susceptibility_by_age_rate_in = susceptibility_params %>% filter(estimate != "mean") %>% pull(rate_in)
influenza_susceptibility_by_age_cutoff_in = susceptibility_params %>% filter(estimate != "mean") %>% pull(cutoff)
influenza_susceptibility_by_age_high_in =  susceptibility_params %>% filter(estimate != "mean") %>% pull(high)
nursing_home_incidence_importations_factor_in = c(0.01, 0.15)
neighborhood_same_age_bias_in = c(0.01,0.1)
   
    
## VARIANTS 
variantalpha_transmissibility_factor_in = c(1,1)
variantalpha_cross_protection_prob_in = c(1,1)
variantalpha_introduction_day_in = c(as.numeric(as.Date('2021-03-10') - as.Date('2020-01-01')),as.numeric(as.Date('2021-02-14') - as.Date('2020-01-01')))
variantalpha_imports_factor_in = c(0.5,10)

variantgamma_severity_factor_in = c(1,2.5)
variantgamma_transmissibility_factor_in = c(1,1)
variantgamma_cross_protection_prob_in = c(1,1)
variantgamma_introduction_day_in = c(as.numeric(as.Date('2020-12-10') - as.Date('2020-01-01')),as.numeric(as.Date('2020-12-14') - as.Date('2020-01-01')))
variantgamma_imports_factor_in = c(0.5,10)

variantkappa_severity_factor_in = c(1,2.5)
variantkappa_transmissibility_factor_in = c(1,1)
variantkappa_cross_protection_prob_in = c(1,1)
variantkappa_introduction_day_in = c(as.numeric(as.Date('2020-11-01') - as.Date('2020-01-01')),as.numeric(as.Date('2021-02-14') - as.Date('2020-01-01')))

variantdelta_severity_factor_in = c(1,1.5)
variantdelta_transmissibility_factor_in = c(1,1)
variantdelta_cross_protection_prob_in = c(1,1)
variantdelta_introduction_day_in = c(as.numeric(as.Date('2021-07-08') - as.Date('2020-01-01')),as.numeric(as.Date('2021-07-09') - as.Date('2020-01-01')))
variantdelta_imports_factor_in = c(0.5,10)


if(variants_in == 1){    
    variantalpha_transmissibility_factor_in = c(1.1,1.5)
    variantalpha_cross_protection_prob_in = c(0.95,0.95)
    variantalpha_introduction_day_in = c(as.numeric(as.Date('2020-12-15') - as.Date('2020-01-01')),as.numeric(as.Date('2021-03-10') - as.Date('2020-01-01')))
    variantalpha_imports_factor_in = c(0.5,10)

    variantgamma_severity_factor_in = c(1,2.5)
    variantgamma_transmissibility_factor_in = c(1.1,2.4)
    variantgamma_cross_protection_prob_in = c(0.54,0.79)
    variantgamma_introduction_day_in = c(as.numeric(as.Date('2020-11-04') - as.Date('2020-01-01')),as.numeric(as.Date('2021-03-24') - as.Date('2020-01-01')))
    variantgamma_imports_factor_in = c(0.5,10)

    variantkappa_severity_factor_in = c(1,2.5)
    variantkappa_transmissibility_factor_in = c(1.0,2.0)
    variantkappa_cross_protection_prob_in = c(0.54,0.90)
    variantkappa_introduction_day_in = c(as.numeric(as.Date('2020-11-11') - as.Date('2020-01-01')),as.numeric(as.Date('2021-03-14') - as.Date('2020-01-01')))

    variantdelta_severity_factor_in = c(1,1.5)
    variantdelta_cross_protection_prob_in = c(0.3,1.0) ## This is temporary
    variantdelta_transmissibility_factor_in = c(1.0,2.1)
    variantdelta_introduction_day_in = c(as.numeric(as.Date('2021-06-08') - as.Date('2020-01-01')),as.numeric(as.Date('2021-06-09') - as.Date('2020-01-01')))
    variantdelta_imports_factor_in = c(0.5,10)
}

enable_age_specific_susceptibility_min_in = 0

influenza_susceptibility_by_age_minage_in = 10
influenza_susceptibility_by_age_minvalue_in = 1

## Read the susceptibility cutoff from an input
influenza_susceptibility_by_age_cutoff_in = c(kids_susceptibility_age_in, kids_susceptibility_age_in)

if(variants_in == 0){
    scalars_sobol_df = sobol_design(
        lower = c(seed=1,imports_factor=imports_factor[1],
                  influenza_transmissibility=influenza_transmissibility[1],
                  shelter_in_place_compliance=shelter_in_place_compliance[1],
                  facemask_compliance = facemask_compliance[1],
                  influenza_susceptibility_by_age_offset = influenza_susceptibility_by_age_offset_in[1],
                  influenza_susceptibility_by_age_rate = influenza_susceptibility_by_age_rate_in[1],
                  influenza_susceptibility_by_age_cutoff = influenza_susceptibility_by_age_cutoff_in[1],
                  influenza_susceptibility_by_age_high = influenza_susceptibility_by_age_high_in[1],
                  nursing_home_incidence_importations_factor = nursing_home_incidence_importations_factor_in[1],
                  neighborhood_same_age_bias = neighborhood_same_age_bias_in[1],
                  workplace_contact_factor = workplace_contact_factor_in[1],
                  neighborhood_contact_factor = neighborhood_contact_factor_in[1],
                  community_contact_rate_1 = community_contact_rate_1_in[1]
                  ),
        upper = c(seed=as.integer(Sys.time()),
                  imports_factor=imports_factor[2],
                  influenza_transmissibility=influenza_transmissibility[2],
                  shelter_in_place_compliance=shelter_in_place_compliance[2],
                  facemask_compliance = facemask_compliance[2],
                  influenza_susceptibility_by_age_offset = influenza_susceptibility_by_age_offset_in[2],
                  influenza_susceptibility_by_age_rate = influenza_susceptibility_by_age_rate_in[2],
                  influenza_susceptibility_by_age_cutoff = influenza_susceptibility_by_age_cutoff_in[2],
                  influenza_susceptibility_by_age_high = influenza_susceptibility_by_age_high_in[2],
                  nursing_home_incidence_importations_factor = nursing_home_incidence_importations_factor_in[2],
                  neighborhood_same_age_bias = neighborhood_same_age_bias_in[2],
                  workplace_contact_factor = workplace_contact_factor_in[2],
                  neighborhood_contact_factor = neighborhood_contact_factor_in[2],
                  community_contact_rate_1 = community_contact_rate_1_in[2]
                  ),
        reps)
    scalars_sobol_df$variantalpha_transmissibility_factor = variantalpha_transmissibility_factor_in[1]
    scalars_sobol_df$variantalpha_cross_protection_prob = variantalpha_cross_protection_prob_in[1]
    scalars_sobol_df$variantgamma_transmissibility_factor = variantgamma_transmissibility_factor_in[1]
    scalars_sobol_df$variantgamma_cross_protection_prob = variantgamma_cross_protection_prob_in[1]
    scalars_sobol_df$variantkappa_transmissibility_factor = variantkappa_transmissibility_factor_in[1]
    scalars_sobol_df$variantkappa_cross_protection_prob = variantkappa_cross_protection_prob_in[1]
    scalars_sobol_df$variantalpha_introduction_day = variantalpha_introduction_day_in[1]
    scalars_sobol_df$variantkappa_introduction_day = variantkappa_introduction_day_in[1]
    scalars_sobol_df$variantgamma_introduction_day = variantgamma_introduction_day_in[1]
    scalars_sobol_df$variantdelta_introduction_day = variantdelta_introduction_day_in[1]
    scalars_sobol_df$variantdelta_transmissibility_factor = mean(variantdelta_transmissibility_factor_in)
    scalars_sobol_df$variantdelta_cross_protection_prob = mean(variantdelta_cross_protection_prob_in)
    scalars_sobol_df$variantkappa_introduction_day = variantkappa_introduction_day_in[1]
    scalars_sobol_df$variantalpha_imports_factor = variantalpha_imports_factor_in[1]
    scalars_sobol_df$variantgamma_imports_factor = variantgamma_imports_factor_in[1]
    scalars_sobol_df$variantdelta_imports_factor = variantdelta_imports_factor_in[1]
    scalars_sobol_df$variantgamma_severity_factor = variantgamma_severity_factor_in[1]
    scalars_sobol_df$variantkappa_severity_factor = variantkappa_severity_factor_in[1]
    scalars_sobol_df$variantdelta_severity_factor = variantdelta_severity_factor_in[1]
}


if(variants_in == 1){
    state_code_in = state_code
    calibration_simdir = sprintf('FRED_%.0f_calibration_asymp_%.2f_fm_%.2f_ksus_%.2f_var_%.0f_vax_%03d_mov',
                                 state_code, 
                                 asymp_infectivity_in,
                                 face_mask_transmission_efficacy_in,
                                 kids_susceptibility_age_in, 
                                 0,
                                 0)

    calibration_dir = file.path(getwd(), 'output','CALIBRATION',sprintf("%s_%s", calibration_simdir, "out"))
    params_df = read_csv(file.path(calibration_dir, 'FRED_parameters_out.csv'))
    fred_sweep_df = read_csv(file.path(calibration_dir, 'fred_output.csv'))
    params_sweep_ll = params_df %>%
        filter(state_code == state_code_in) %>%
        mutate(LL_total = LL_deaths)


    ## 2. Sample based on their likelihood to ensure 10% sampled
    # particles_w = 0.02
    # prob_array = exp(-params_sweep_ll$LL_total*particles_w)
    # if(sum(prob_array) > 0){
    #     indx_sampled = sample.int(n = length(prob_array), size = reps, prob = prob_array, replace = T)
    #     n_sampled = length(unique(indx_sampled))
    # }    
    
    indx_sampled = c(rep(1544, 667), rep(1683,667), rep(1712,666))
    print(params_sweep_ll$job_id[unique(indx_sampled)])
    n_sampled = length(indx_sampled)
    print(n_sampled)
    scalars_sampled = params_sweep_ll[indx_sampled,] %>%
        dplyr::select(imports_factor, influenza_transmissibility,
                      shelter_in_place_compliance, facemask_compliance,                  
                      influenza_susceptibility_by_age_offset,influenza_susceptibility_by_age_rate,
                      influenza_susceptibility_by_age_cutoff,influenza_susceptibility_by_age_high,
                      nursing_home_incidence_importations_factor,neighborhood_same_age_bias,
                      workplace_contacts, office_contacts, workplace_contact_factor,
                      neighborhood_contacts, neighborhood_contact_factor,holiday_contact_rate,
                      community_contact_rate_1,
                      school_contacts, classroom_contacts,school_contact_factor
                      ) 

    scalars_sobol_df = sobol_design(
        lower = c(seed=1,
                  variantalpha_transmissibility_factor = variantalpha_transmissibility_factor_in[1],
                  variantalpha_cross_protection_prob = variantalpha_cross_protection_prob_in[1],
                  variantgamma_transmissibility_factor = variantgamma_transmissibility_factor_in[1],
                  variantgamma_severity_factor = variantgamma_severity_factor_in[1],
                  variantgamma_cross_protection_prob = variantgamma_cross_protection_prob_in[1],
                  variantkappa_transmissibility_factor = variantkappa_transmissibility_factor_in[1],
                  variantkappa_severity_factor = variantkappa_severity_factor_in[1],
                  variantkappa_cross_protection_prob = variantkappa_cross_protection_prob_in[1],
                  variantalpha_imports_factor = variantalpha_imports_factor_in[1],
                  variantkappa_introduction_day = variantkappa_introduction_day_in[1],
                  variantgamma_imports_factor = variantgamma_imports_factor_in[1],
                  variantdelta_severity_factor = variantdelta_severity_factor_in[1],
                  variantdelta_transmissibility_factor = variantdelta_transmissibility_factor_in[1],
                  variantdelta_imports_factor = variantdelta_imports_factor_in[1]
                  ),
        upper = c(seed=as.integer(Sys.time()),
                  variantalpha_transmissibility_factor = variantalpha_transmissibility_factor_in[2],
                  variantalpha_cross_protection_prob = variantalpha_cross_protection_prob_in[2],
                  variantgamma_transmissibility_factor = variantgamma_transmissibility_factor_in[2],
                  variantgamma_severity_factor = variantgamma_severity_factor_in[2],
                  variantgamma_cross_protection_prob = variantgamma_cross_protection_prob_in[2],
                  variantkappa_transmissibility_factor = variantkappa_transmissibility_factor_in[2],
                  variantkappa_severity_factor = variantkappa_severity_factor_in[2],
                  variantkappa_cross_protection_prob = variantkappa_cross_protection_prob_in[2],
                  variantalpha_imports_factor = variantalpha_imports_factor_in[2],
                  variantkappa_introduction_day = variantkappa_introduction_day_in[2],
                  variantgamma_imports_factor = variantgamma_imports_factor_in[2],
                  variantdelta_severity_factor = variantdelta_severity_factor_in[2],
                  variantdelta_transmissibility_factor = variantdelta_transmissibility_factor_in[2],
                  variantdelta_imports_factor = variantdelta_imports_factor_in[2]),
        reps)
    ##variantdelta_cross_protection_prob = variantdelta_cross_protection_prob_in[2],    
    ##scalars_sobol_df$variantdelta_transmissibility_factor = mean(variantdelta_transmissibility_factor_in)
    ##scalars_sobol_df$variantdelta_imports_factor = mean(variantdelta_imports_factor_in)

    
    scalars_sobol_df$variantdelta_cross_protection_prob = 1 - rbeta(reps, 0.9, 8.2)
    scalars_sobol_df$variantdelta_cross_protection_prob[scalars_sobol_df$variantdelta_cross_protection_prob < variantdelta_cross_protection_prob_in[1]] = variantdelta_cross_protection_prob_in[1]

    ## variantdelta_cross_protection_prob_in[2]
    scalars_sobol_df = bind_cols(scalars_sobol_df, scalars_sampled)
}
scalars_sobol_df$variantalpha_transmissibility =  scalars_sobol_df$influenza_transmissibility * scalars_sobol_df$variantalpha_transmissibility_factor
scalars_sobol_df$variantkappa_transmissibility =  scalars_sobol_df$influenza_transmissibility * scalars_sobol_df$variantkappa_transmissibility_factor
scalars_sobol_df$variantgamma_transmissibility =  scalars_sobol_df$influenza_transmissibility * scalars_sobol_df$variantgamma_transmissibility_factor
scalars_sobol_df$variantdelta_transmissibility =  scalars_sobol_df$influenza_transmissibility * scalars_sobol_df$variantdelta_transmissibility_factor    

scalars_sobol_df$variantalpha_susceptibility_by_age_offset  = scalars_sobol_df$influenza_susceptibility_by_age_offset 
scalars_sobol_df$variantalpha_susceptibility_by_age_rate = scalars_sobol_df$influenza_susceptibility_by_age_rate 
scalars_sobol_df$variantalpha_susceptibility_by_age_cutoff = scalars_sobol_df$influenza_susceptibility_by_age_cutoff
scalars_sobol_df$variantalpha_susceptibility_by_age_high = scalars_sobol_df$influenza_susceptibility_by_age_high

scalars_sobol_df$variantkappa_susceptibility_by_age_offset  = scalars_sobol_df$influenza_susceptibility_by_age_offset 
scalars_sobol_df$variantkappa_susceptibility_by_age_rate = scalars_sobol_df$influenza_susceptibility_by_age_rate 
scalars_sobol_df$variantkappa_susceptibility_by_age_cutoff = scalars_sobol_df$influenza_susceptibility_by_age_cutoff
scalars_sobol_df$variantkappa_susceptibility_by_age_high = scalars_sobol_df$influenza_susceptibility_by_age_high

scalars_sobol_df$variantgamma_susceptibility_by_age_offset  = scalars_sobol_df$influenza_susceptibility_by_age_offset 
scalars_sobol_df$variantgamma_susceptibility_by_age_rate = scalars_sobol_df$influenza_susceptibility_by_age_rate 
scalars_sobol_df$variantgamma_susceptibility_by_age_cutoff = scalars_sobol_df$influenza_susceptibility_by_age_cutoff
scalars_sobol_df$variantgamma_susceptibility_by_age_high = scalars_sobol_df$influenza_susceptibility_by_age_high

scalars_sobol_df$variantdelta_susceptibility_by_age_offset  = scalars_sobol_df$influenza_susceptibility_by_age_offset 
scalars_sobol_df$variantdelta_susceptibility_by_age_rate = scalars_sobol_df$influenza_susceptibility_by_age_rate 
scalars_sobol_df$variantdelta_susceptibility_by_age_cutoff = scalars_sobol_df$influenza_susceptibility_by_age_cutoff
scalars_sobol_df$variantdelta_susceptibility_by_age_high = scalars_sobol_df$influenza_susceptibility_by_age_high

scalars_sobol_df$school_contact_factor = 1.0
scalars_sobol_df$holiday_contact_rate = 1.0

scalars_sobol_df = bind_cols(scalars_sobol_df, get_ifr_param_string(scalars_sobol_df$variantgamma_severity_factor, 'gamma'))
scalars_sobol_df = bind_cols(scalars_sobol_df, get_ifr_param_string(scalars_sobol_df$variantkappa_severity_factor, 'kappa'))
scalars_sobol_df = bind_cols(scalars_sobol_df, get_ifr_param_string(scalars_sobol_df$variantdelta_severity_factor, 'delta'))

## get_ifr_param_string(scalars_sobol_df$variantdelta_severity_factor[4], 'delta')

##==============================================#
## Initial conditions-------------------
##==============================================#
##initial_inf_file = sprintf('input_files/%d_imports_alternative.csv', state_code)
initial_inf_file = sprintf('input_files/%d_imports_combined.csv', state_code)

primary_cases_file = file.path(output.dir, sprintf('initial_cases_%d_%d.txt',state_code,1:reps))

initial_df = read_csv(initial_inf_file)
initial_df$day = as.numeric(difftime(initial_df$Date, as.Date(start_date), units='days'))

variants_imp_file = 'input_files/COL_variant_imports.csv'
variants_imp_df = read_csv(variants_imp_file)

## Sample 'reps' from the initial conditions
## replicate_init = sample.int(n=length(unique(initial_df$Replicate)), size = reps, replace=T)

for(nn in 1:reps){
    ## tmp_df = filter(initial_df, Replicate == replicate_init[nn]) %>%
    ##     arrange(day) %>% dplyr::select(Imports, day)

    ## For now, choose the mean
    tmp_df = initial_df %>% group_by(day) %>% summarize(Imports = ceiling(mean(Imports))) %>% ungroup()
    
    tmp_df$Imports = round(tmp_df$Imports * scalars_sobol_df$imports_factor[nn])
    month_imports = ceiling(mean(tmp_df$Imports[(tmp_df$day <= (max(tmp_df$day) - 30)) & (tmp_df$Imports > 0)]))

    ## Now, let's add cases by week after initial cases            
    extra_rows = numDays - nrow(tmp_df)
    weekly_imports = month_imports * 7
    if(extra_rows > 0){
        week_imp = data.frame(Imports = rep(weekly_imports, floor(extra_rows/7)), stringsAsFactors = F)
        week_imp$day = seq(from=max(tmp_df$day) + 1, by = 7, length.out = nrow(week_imp))
        tmp_df = bind_rows(tmp_df, week_imp)
    }
    if(variants_in == 1){
        alpha_imports = filter(variants_imp_df, variant == "20I (Alpha, V1)") %>%
            mutate(day = as.numeric(Date - as.Date('2020-01-01'))) %>%
            mutate(Imports = round(scalars_sobol_df$variantalpha_imports_factor[nn] * TotalVariantImports)) %>%
            filter(Imports >= 1)

        gamma_imports = filter(variants_imp_df, variant == "20J (Gamma, V3)") %>%
            mutate(day = as.numeric(Date - as.Date('2020-01-01'))) %>%
            mutate(Imports = round(scalars_sobol_df$variantgamma_imports_factor[nn] * TotalVariantImports)) %>%
            filter(Imports >= 1,Date <= as.Date('2021-04-01'))

        delta_imports = filter(variants_imp_df, variant == "21A (Delta)") %>%
            mutate(day = as.numeric(Date - as.Date('2020-01-01'))) %>%
            mutate(Imports = round(scalars_sobol_df$variantdelta_imports_factor[nn] * TotalVariantImports)) %>%
            filter(Imports >= 1, Date >= as.Date('2021-04-01')) %>%
            mutate(day = day + 60 + 15) ## Temporary delay
        
        tmp_df$Date = tmp_df$day + as.Date('2020-01-01')
        tmp_df_main = tmp_df[tmp_df$Date < as.Date('2020-12-10'),]

        tmp_df_variantkappa = tmp_df[tmp_df$Date >= as.Date('2020-01-01') + round(scalars_sobol_df$variantkappa_introduction_day[nn]) & tmp_df$Date <= as.Date('2021-06-01'),]
    
        init_cases_base_lines = sprintf('%.0f %.0f %.0f 0 1 %.0f', tmp_df_main$day, tmp_df_main$day, tmp_df_main$Imports, tmp_df_main$Imports)
    
        init_cases_alpha_lines = sprintf('%.0f %.0f %.0f 1 1 %.0f', alpha_imports$day, alpha_imports$day, alpha_imports$Imports, alpha_imports$Imports)
        
        init_cases_gamma_lines = sprintf('%.0f %.0f %.0f 2 1 %.0f', gamma_imports$day, gamma_imports$day, gamma_imports$Imports, gamma_imports$Imports)

        init_cases_kappa_lines = sprintf('%.0f %.0f %.0f 3 1 %.0f', tmp_df_variantkappa$day, tmp_df_variantkappa$day, tmp_df_variantkappa$Imports, tmp_df_variantkappa$Imports)
        
        init_cases_delta_lines = sprintf('%.0f %.0f %.0f 4 1 %.0f', delta_imports$day, delta_imports$day, delta_imports$Imports, delta_imports$Imports)
        
        init_cases_lines = c(init_cases_base_lines, init_cases_alpha_lines,init_cases_kappa_lines,init_cases_gamma_lines, init_cases_delta_lines)
    }else{
        init_cases_lines = sprintf('%.0f %.0f %.0f 0 1 %.0f', tmp_df$day, tmp_df$day, tmp_df$Imports, tmp_df$Imports)
    }
    
    fileConn<-file(primary_cases_file[nn])
    writeLines(init_cases_lines, fileConn)
    close(fileConn)
}


##==============================================#
## School closure-------------------
##==============================================#
schools_open_list = read_csv('./input_files/11001_schools_open_gps.csv')
schools_open_list$PropOpen[schools_open_list$PropOpen > 1] = 1.0
localidad_esc = read_csv('./input_files/Localidad_Unidad_Catastral.csv')
localidad_list = 1:19

school_schedule_closed_file = file.path(output.dir, sprintf('school_schedule_closed_%d.txt',state_code))

## Parse updated list of schools reopened
## Process school reopen files and transform into Unidad Catastral capacity
sim_start_date = as.Date(start_date)
school_reopen_list_file = './input_files/11001_schools_open_gps_2021.csv'
school_reopen_df = read_csv(school_reopen_list_file) %>%
    group_by(Grade, zipcode, start_date, end_date) %>%
    summarize(Capacity = ifelse(InPerson<Total_students, InPerson/Total_students, 1.0)) %>%
    ungroup() %>%
    mutate(MinAge = ifelse(Grade == 'PREK', 0,
                    ifelse(Grade == 'PRIMARY', 6,
                    ifelse(Grade == 'SECONDARY',11,18))),
           MaxAge = ifelse(Grade == 'PREK', 5,
                    ifelse(Grade == 'PRIMARY', 10,
                    ifelse(Grade == 'SECONDARY',17,20)))) %>%
    mutate(start_day = as.numeric(start_date - sim_start_date),
           end_day = as.numeric(end_date - sim_start_date)) %>%
    mutate(end_day = ifelse(end_day == max(end_day), numDays, end_day)) %>%
    replace_na(list(zipcode = "")) %>%
    filter(zipcode != "Bogota")

    
## Write files
{
    initial_closure = sprintf("%d %d 1 20 0", school_closure_day_in,  as.integer(as.Date('2020-12-15') - as.Date(start_date)))
    ## Make sure this is working properly
    vacation_closure = sprintf("%d %d 1 20 0", as.integer(as.Date('2020-12-15') - as.Date(start_date)), numDays)
    school_baseline_reopen_date = as.integer(early_school_vacation_end - as.Date(start_date)) + 1
    school_current_early_open_date = as.integer(current_open_date - as.Date(start_date)) + 1
    ## Schools already open in 2020
    current_open_school_lines = c(initial_closure, vacation_closure)
    for(ll in 1:length(localidad_list)){
        localidad_tmp = filter(schools_open_list, Localidad_ID == localidad_list[ll])        
        esc_in_loc = sprintf("11001%s", localidad_esc$SCACODIGO[localidad_esc$Localidad == localidad_list[ll]])
        
        current_open_school_lines = c(current_open_school_lines,
                                      sprintf("%d %d 18 20 %0.4f 0 %s",
                                              school_current_early_open_date,
                                              as.integer(as.Date('2020-12-15') - as.Date(start_date)),
                                              localidad_tmp[localidad_tmp$Grade == 'university','PropOpen'], esc_in_loc),
                                      sprintf("%d %d 1 17 %.4f 0 %s",
                                              school_current_early_open_date,
                                              as.integer(as.Date('2020-12-15') - as.Date(start_date)),
                                              localidad_tmp[localidad_tmp$Grade == 'basic_ed','PropOpen'], esc_in_loc))
    }   
    ## 2021 reopened schools
    tmp_school_open = filter(school_reopen_df, Grade != 'UNIVERSITY')
    tmp_college_open = filter(school_reopen_df, Grade == 'UNIVERSITY')
    current_open_school_lines = c(current_open_school_lines, sprintf("%d %d %d %d %.4f 0 11001%s",tmp_school_open$start_day,tmp_school_open$end_day, tmp_school_open$MinAge, tmp_school_open$MaxAge, tmp_school_open$Capacity, tmp_school_open$zipcode))
    current_open_school_lines = c(current_open_school_lines,
                                  sprintf("%d %d 18 20 %.4f", tmp_college_open$start_day, tmp_college_open$end_day, tmp_college_open$Capacity))

    
    ## Write reopening files    
    fileConn<-file(school_schedule_closed_file)
    writeLines(current_open_school_lines, fileConn)
    close(fileConn)
}

##==============================================#
## Community increase-------------------
##==============================================#
if(variants_in == 0){
    community_timeseries_df = read_csv('./input_files/interventions_covid_timevarying_community.csv')
}else{
    community_timeseries_df = read_csv('./input_files/interventions_covid_timevarying_community_baseline.csv')
}
community_timeseries_file = file.path(output.dir, sprintf("community_timeseries_%d.txt", 1:reps))
community_timeseries_df$day = as.numeric(community_timeseries_df$date - as.Date(start_date))
community_timeseries_df$community_trend[community_timeseries_df$date < as.Date('2020-05-15')] = community_timeseries_df$community_trend[community_timeseries_df$date == as.Date('2020-05-15')]

community_timeseries_df$community_trend = community_timeseries_df$community_trend - community_timeseries_df$community_trend[1]
for(nn in 1:reps){        
    tmp_comm_df = community_timeseries_df
    ##community_timeseries_df$contact_rate = 1 - (1 - community_timeseries_df$community_trend) * scalars_sobol_df$community_contact_rate_1[nn]
    tmp_comm_df$contact_rate = tmp_comm_df$community_trend * scalars_sobol_df$community_contact_rate_1[nn] + 1
    
    tmp_comm_df$contact_rate[tmp_comm_df$contact_rate < 0] = 0

    community_time_lines = sprintf('%.0f %.0f %.4f', tmp_comm_df$day, tmp_comm_df$day, tmp_comm_df$contact_rate)

    fileConn<-file(community_timeseries_file[nn])
    writeLines(community_time_lines, fileConn)
    close(fileConn)
}

##==============================================#
## Vaccine timeseries-------------------
##==============================================#
vaccine_daily_capacity_file = file.path(output.dir, "vaccination_daily_capacity_timeseries.txt")
vaccine_stock_file = file.path(output.dir, "vaccination_stock_timeseries.txt")

vaccine_stock_df = read.csv('./input_files/11001_vaccine_stock_timeseries.csv') %>%
    drop_na()
vaccine_stock_df$day = as.numeric(as.Date(vaccine_stock_df$Date) - as.Date(start_date))
vaccine_daily_df = read.csv('./input_files/11001_vaccine_capacity_timeseries.csv')
vaccine_daily_df$day = as.numeric(as.Date(vaccine_daily_df$Date) - as.Date(start_date))
vaccine_daily_df$VaccinesApplied = round(vaccine_daily_df$VaccinesApplied)

vaccine_stock_lines = sprintf("%.0f %.0f %.0f %.0f",vaccine_stock_df$day, vaccine_stock_df$day,vaccine_stock_df$TotalVaccines,vaccine_stock_df$ID)

fileConn<-file(vaccine_stock_file)
writeLines(vaccine_stock_lines, fileConn)
close(fileConn)

vaccine_capacity_lines = c(sprintf("%.0f %.0f",vaccine_daily_df$day,vaccine_daily_df$VaccinesApplied))

fileConn<-file(vaccine_daily_capacity_file)
writeLines(vaccine_capacity_lines, fileConn)
close(fileConn)

##==============================================#
## Facemask compliance-------------------
##==============================================#
school_facemask_compliance = 0.75
facemask_time_file = file.path(output.dir, sprintf('facemask_compliance_timeseries_%d_%d.txt',state_code,1:reps))

facemask_timeseries_df = read_csv('./input_files/facemask_timeseries_compliance.csv') %>%
    dplyr::select(-Day)
facemask_timeseries_df$day = as.numeric(difftime(facemask_timeseries_df$Date, as.Date(start_date), units='days'))

for(nn in 1:reps){
    ## For now, choose the mean
    tmp_df = facemask_timeseries_df %>% filter(State == tolower(interventions_st_df$state_name[interventions_st_df$State == state_code])) %>%
        mutate(FacemaskTrends = FacemaskTrends * scalars_sobol_df$facemask_compliance[nn]) %>%
        dplyr::select(day, FacemaskTrends)
    if(max(tmp_df$day) < numDays){
        tmp_df2 = data.frame(day = seq(from=max(tmp_df$day) + 1, to = numDays),
                             FacemaskTrends = tail(tmp_df$FacemaskTrends,1), stringsAsFactors = F)
        tmp_df = bind_rows(tmp_df, tmp_df2)
    }

    compliance_time_lines = c()
    for(loc in c("workplace", "office", "other")){
        compliance_time_lines = c(compliance_time_lines,
                                  sprintf('%.0f %.0f %.4f %s', tmp_df$day, tmp_df$day, tmp_df$FacemaskTrends, loc)
                                  )
    }

    ## Schools start facemask compliance when they open
    facemask_school_day = as.integer(current_open_date - as.Date(start_date))
    school_high_lines = c()
    
    for(loc in c("school", "classroom")){
        school_high_lines = c(school_high_lines, sprintf('%.0f %.0f %.4f %s', facemask_school_day, numDays, school_facemask_compliance, loc))
    }

    school_high_lines = c(compliance_time_lines, school_high_lines)
    
    ## Schools start facemask compliance when they open    
    fileConn<-file(facemask_time_file[nn])
    writeLines(school_high_lines, fileConn)
    close(fileConn)    
}


##==============================================#
## Shelter in place-------------------
##==============================================#
post_lockdown_mobility_in = 0.1
shelter_time_file = file.path(output.dir, sprintf('shelter_timeseries_%d_%d.txt',state_code,1:reps))

## Moving away from google's data to Grandata census-tract specific
## shelter_timeseries_df = read_csv('./input_files/interventions_covid_timevarying_shelter.csv')
shelter_timeseries_df = read_csv('./input_files/11001_mobility_trends.csv')
shelter_timeseries_df$day = as.numeric(difftime(shelter_timeseries_df$date, as.Date(start_date), units='days'))


state_shelter_df = shelter_timeseries_df %>%
    dplyr::filter(State == state_code, replicate == 1) %>%
    mutate(StateCod = sprintf('%d%s',State,SCACODIGO),
           minAge = 0, maxAge = 120) %>%
    dplyr::select(day,shelter_trend,minAge,maxAge,StateCod)
        
max_tmp_day = max(state_shelter_df$day)
if(max_tmp_day > numDays){max_tmp_day = numDays}
start_t = Sys.time()
##foreach(nn = 1:reps) %dopar% {

for(nn in 1:reps){    
    ## For now, choose one replicate
    ##start_t = Sys.time()
    tmp_df = state_shelter_df %>% 
        mutate(shelter_trend = shelter_trend * scalars_sobol_df$shelter_in_place_compliance[nn])
    tmp_df$day2 = tmp_df$day
    
    tmp_df_tail = tmp_df[tmp_df$day == max_tmp_day,]
    tmp_df_tail$day = max_tmp_day + 1
    tmp_df_tail$shelter_trend = post_lockdown_mobility_in * tmp_df_tail$shelter_trend
    tmp_df_tail$day2 = tmp_df_tail$day
    
    tmp_df$shelter_trend[tmp_df$shelter_trend < 0.0001] = 0.0
    tmp_df$shelter_trend[tmp_df$shelter_trend > 1.0] = 1.0
    tmp_df = bind_rows(tmp_df, tmp_df_tail)
    
    write_delim(tmp_df[c('day','day2','shelter_trend','minAge','maxAge','StateCod')],file = shelter_time_file[nn],col_names = F)
    
    ##shelter_time_lines = sprintf('%.0f %.0f %.4f 0 120 %s', tmp_df$day, tmp_df$day, tmp_df$shelter_trend, tmp_df$StateCod)    
    ##ashelter_time_lines = sapply(1:nrow(tmp_df),function(x){sprintf('%.0f %.0f %.4f 0 120 %d%s', tmp_df$day[x], tmp_df$day[x], tmp_df$shelter_trend[x], tmp_df$State[x], tmp_df$SCACODIGO[x])})
    ##a = do.call(sprintf,c('%.0f %.0f %.4f 0 120 %d%s',
    ##shelter_time_tail_lines = sprintf('%.0f %.0f %.4f 0 120 %s', tmp_df_tail$day, numDays, tmp_df_tail$shelter_trend, tmp_df_tail$StateCod)

    ## Add shelter by age
    ## shelter_age_lines = sprintf("82 212 %.4f 60 120", scalars_sobol_df$shelter_in_place_age_compliance[nn])
    ## fileConn<-file(shelter_time_file[nn])
    ## writeLines(c(shelter_time_lines, shelter_time_tail_lines), fileConn)
    ## close(fileConn)
}

end_t = Sys.time()
print(sprintf("Creating shelter files took %.04f", as.numeric(end_t - start_t)))

##==============================================#
## fixed parameters-------------------
##==============================================#
## IMPORTANT: The offset delays the epidemic. Does not include added days of imported cases
epidemic_offset = 0
num_demes = 1

##scalars_intervention = data.frame(stringsAsFactors=F)

scalars_intervention = scalars_sobol_df %>%
    mutate(
        seed = floor(scalars_sobol_df$seed),
        primary_cases_file = primary_cases_file,
        imports_factor = scalars_sobol_df$imports_factor,
        shelter_in_place_file = shelter_time_file,
        shelter_in_place_compliance = scalars_sobol_df$shelter_in_place_compliance,
        face_mask_timeseries_file = facemask_time_file,
        facemask_compliance = scalars_sobol_df$facemask_compliance,
        influenza_transmissibility = scalars_sobol_df$influenza_transmissibility,
        influenza_susceptibility_by_age_offset = scalars_sobol_df$influenza_susceptibility_by_age_offset,
        influenza_susceptibility_by_age_rate = scalars_sobol_df$influenza_susceptibility_by_age_rate,
        influenza_susceptibility_by_age_cutoff = scalars_sobol_df$influenza_susceptibility_by_age_cutoff,
        influenza_susceptibility_by_age_high = scalars_sobol_df$influenza_susceptibility_by_age_high,
        nursing_home_incidence_importations_factor = scalars_sobol_df$nursing_home_incidence_importations_factor,
        neighborhood_same_age_bias = scalars_sobol_df$neighborhood_same_age_bias,
        school_contacts = school_contacts_in * scalars_sobol_df$school_contact_factor,
        classroom_contacts = classroom_contacts_in * scalars_sobol_df$school_contact_factor,
        school_contact_factor = scalars_sobol_df$school_contact_factor,
        workplace_contacts = workplace_contacts_in * scalars_sobol_df$workplace_contact_factor,
        neighborhood_contacts = neighborhood_contacts_in * scalars_sobol_df$neighborhood_contact_factor,
        neighborhood_contact_factor = scalars_sobol_df$neighborhood_contact_factor,
        office_contacts = office_contacts_in * scalars_sobol_df$workplace_contact_factor,
        workplace_contact_factor = scalars_sobol_df$workplace_contact_factor,
        holiday_contact_rate = scalars_sobol_df$holiday_contact_rate,
        community_contact_rate_1 = scalars_sobol_df$community_contact_rate_1,
        community_contact_timeseries_file = community_timeseries_file,
        
        variantalpha_transmissibility =  scalars_sobol_df$variantalpha_transmissibility,
        variantalpha_transmissibility_factor =  scalars_sobol_df$variantalpha_transmissibility_factor,
        variantalpha_cross_protection_prob = scalars_sobol_df$variantalpha_cross_protection_prob,
        variantalpha_imports_factor = scalars_sobol_df$variantalpha_imports_factor,
        variantalpha_susceptibility_by_age_offset  = scalars_sobol_df$variantalpha_susceptibility_by_age_offset,
        variantalpha_susceptibility_by_age_rate = scalars_sobol_df$variantalpha_susceptibility_by_age_rate,
        variantalpha_susceptibility_by_age_cutoff = scalars_sobol_df$variantalpha_susceptibility_by_age_cutoff,
        variantalpha_susceptibility_by_age_high = scalars_sobol_df$variantalpha_susceptibility_by_age_high,
        variantgamma_transmissibility =  scalars_sobol_df$variantgamma_transmissibility,
        
        variantgamma_transmissibility_factor =  scalars_sobol_df$variantgamma_transmissibility_factor,
        variantgamma_cross_protection_prob = scalars_sobol_df$variantgamma_cross_protection_prob,
        variantgamma_imports_factor = scalars_sobol_df$variantgamma_imports_factor,
        variantgamma_susceptibility_by_age_offset  = scalars_sobol_df$variantgamma_susceptibility_by_age_offset,
        variantgamma_susceptibility_by_age_rate = scalars_sobol_df$variantgamma_susceptibility_by_age_rate,
        variantgamma_susceptibility_by_age_cutoff = scalars_sobol_df$variantgamma_susceptibility_by_age_cutoff,
        variantgamma_susceptibility_by_age_high = scalars_sobol_df$variantgamma_susceptibility_by_age_high,
        
        variantkappa_transmissibility =  scalars_sobol_df$variantkappa_transmissibility,
        variantkappa_transmissibility_factor =  scalars_sobol_df$variantkappa_transmissibility_factor,
        variantkappa_cross_protection_prob = scalars_sobol_df$variantkappa_cross_protection_prob,
        variantkappa_introduction_day = scalars_sobol_df$variantkappa_introduction_day,
        variantkappa_susceptibility_by_age_offset  = scalars_sobol_df$variantkappa_susceptibility_by_age_offset,
        variantkappa_susceptibility_by_age_rate = scalars_sobol_df$variantkappa_susceptibility_by_age_rate,
        variantkappa_susceptibility_by_age_cutoff = scalars_sobol_df$variantkappa_susceptibility_by_age_cutoff,
        variantkappa_susceptibility_by_age_high = scalars_sobol_df$variantkappa_susceptibility_by_age_high,
        
        variantdelta_transmissibility =  scalars_sobol_df$variantdelta_transmissibility,
        variantdelta_transmissibility_factor =  scalars_sobol_df$variantdelta_transmissibility_factor,
        variantdelta_cross_protection_prob = scalars_sobol_df$variantdelta_cross_protection_prob,
        variantdelta_imports_factor = scalars_sobol_df$variantdelta_imports_factor,
        variantdelta_susceptibility_by_age_offset  = scalars_sobol_df$variantdelta_susceptibility_by_age_offset,
        variantdelta_susceptibility_by_age_rate = scalars_sobol_df$variantdelta_susceptibility_by_age_rate,
        variantdelta_susceptibility_by_age_cutoff = scalars_sobol_df$variantdelta_susceptibility_by_age_cutoff,
        variantdelta_susceptibility_by_age_high = scalars_sobol_df$variantdelta_susceptibility_by_age_high)



##==============================================#
## Create parameters to sweep-----------------
##==============================================#
## TODO:
## 1. Create DF with parameters
scalars = scalars_intervention %>%
    mutate(days = numDays,
           track_infection_events = track_infection_events_in,
           enable_age_specific_susceptibility = enable_age_specific_susceptibility_in,
           influenza_asymp_infectivity = asymp_infectivity_in,
           enable_age_specific_susceptibility_min = enable_age_specific_susceptibility_min_in,
           influenza_susceptibility_by_age_minage = influenza_susceptibility_by_age_minage_in,
           influenza_susceptibility_by_age_minvalue = influenza_susceptibility_by_age_minvalue_in,
           influenza_face_mask_transmission_efficacy = face_mask_transmission_efficacy_in,
           
           variantalpha_asymp_infectivity = asymp_infectivity_in,
           variantalpha_susceptibility_by_age_minage = influenza_susceptibility_by_age_minage_in,
           variantalpha_susceptibility_by_age_minvalue = influenza_susceptibility_by_age_minvalue_in,
           variantalpha_face_mask_transmission_efficacy = face_mask_transmission_efficacy_in,
           
           variantgamma_asymp_infectivity = asymp_infectivity_in,
           variantgamma_susceptibility_by_age_minage = influenza_susceptibility_by_age_minage_in,
           variantgamma_susceptibility_by_age_minvalue = influenza_susceptibility_by_age_minvalue_in,
           variantgamma_face_mask_transmission_efficacy = face_mask_transmission_efficacy_in,
           
           variantkappa_asymp_infectivity = asymp_infectivity_in,
           variantkappa_susceptibility_by_age_minage = influenza_susceptibility_by_age_minage_in,
           variantkappa_susceptibility_by_age_minvalue = influenza_susceptibility_by_age_minvalue_in,
           variantkappa_face_mask_transmission_efficacy = face_mask_transmission_efficacy_in,
           
           variantdelta_asymp_infectivity = asymp_infectivity_in,
           variantdelta_susceptibility_by_age_minage = influenza_susceptibility_by_age_minage_in,
           variantdelta_susceptibility_by_age_minvalue = influenza_susceptibility_by_age_minvalue_in,
           variantdelta_face_mask_transmission_efficacy = face_mask_transmission_efficacy_in, 
           
           enable_nursing_homes_importations = enable_nursing_homes_importations_in,
           enable_transmission_bias = enable_transmission_bias_in,
           track_fatality_events = track_fatality_events_in,
           report_age_of_infection = report_age_of_infection_in,
           report_incidence_by_county = report_incidence_by_county_in,
           report_place_of_infection = report_place_of_infection_in,
           isolation_rate = isolation_rate,
           shelter_in_place_students = shelter_in_place_students,
           num_demes = num_demes,
           synthetic_population_id = synthetic_population_id,
           advance_seeding = advance_seeding,
           epidemic_offset = epidemic_offset,
           enable_shelter_in_place = enable_shelter_in_place,
           enable_shelter_in_place_timeseries = enable_shelter_in_place_timeseries,
           shelter_in_place_delay_mean = shelter_in_place_delay_mean,
           school_closure_policy = school_closure_policy_in,
           school_reduced_capacity = school_reduced_capacity_in,
           enable_school_reduced_capacity = enable_school_reduced_capacity_in,
           school_global_schedule_file = school_schedule_closed_file,
           school_student_teacher_ratio = school_student_teacher_ratio_in,
           school_closure_day = school_closure_day_in,
           school_closure_duration = school_closure_duration_in,
           enable_face_mask_timeseries = enable_face_mask_timeseries_in,
           min_age_face_masks = min_age_face_masks_in,
           enable_face_mask_usage = enable_face_mask_usage_in,
           holiday_start = holiday_start_in,
           holiday_end = holiday_end_in,
           enable_holiday_contacts = enable_holiday_contacts_in,
           enable_community_contact_timeseries = enable_community_contact_timeseries_in,
           start_date = start_date)
if(variants_in == 1){
    scalars$diseases = 5
    scalars$enable_disease_cross_protection = 1
    scalars$influenza_cross_protection_prob = 1.0
    scalars$disease_names = "5 influenza variantalpha variantgamma variantkappa variantdelta"

    ## Vaccination
    scalars$enable_vaccination = 1
    scalars$vaccination_capacity_file = vaccine_daily_capacity_file
    scalars$vaccine_stock_timeseries_file = vaccine_stock_file
    scalars$vaccination_phases_enable_discrete_timing = 1
}else{
    scalars$enable_vaccination = 0
}

##===============================================##
## Write the parameters to files---------------
##===============================================##
defaults_params = './input_files/params_covid.txt'
defaults_covid_params = './input_files/params_covid.txt'

defaults_alpha_params = './input_files/params_covid_alpha.txt'
defaults_gamma_params = './input_files/params_covid_gamma.txt'
defaults_kappa_params = './input_files/params_covid_kappa.txt'
defaults_delta_params = './input_files/params_covid_delta.txt'
defaults_vaccine_params = './input_files/params_covid_vaccine.txt'

if(variants_in == 1){
    if(vaccination_in > 0){
        defaults_vaccine_params = sprintf('./input_files/params_covid_vaccine_%d.txt', vaccination_in)
    }
    defaults_params = './input_files/params_covid_combined.txt'
    system(sprintf("cat %s %s %s %s %s %s > %s",defaults_covid_params, 
                   defaults_alpha_params,
                   defaults_gamma_params, 
                   defaults_kappa_params,
                   defaults_delta_params, 
                   defaults_vaccine_params, 
                   defaults_params, intern = T))
}

basename_params = sprintf('covid_%s_params',state_code)
basename_jobs = sprintf('FRED_%s_calibration',state_code)
write_fred_parameters(scalars, defaults_params, output.dir,basename.in=basename_params, fred_defaults = fred_defaults)

## print report scalars parameters file with IDs
report_scalars = dplyr::select(
                            scalars, influenza_transmissibility,
                            influenza_asymp_infectivity,
                            shelter_in_place_compliance,influenza_face_mask_transmission_efficacy,
                            enable_age_specific_susceptibility_min,
                            enable_age_specific_susceptibility,
                            influenza_susceptibility_by_age_minage,
                            influenza_susceptibility_by_age_minvalue,
                            enable_face_mask_usage, enable_face_mask_timeseries,
                            facemask_compliance, min_age_face_masks,
                            enable_shelter_in_place_timeseries, start_date,
                            community_contact_rate_1, enable_community_contact_timeseries,
                            nursing_home_incidence_importations_factor,
                            neighborhood_same_age_bias,
                            influenza_susceptibility_by_age_rate,
                            influenza_susceptibility_by_age_offset,
                            influenza_susceptibility_by_age_cutoff,
                            influenza_susceptibility_by_age_high,
                            variantalpha_transmissibility,
                            variantalpha_transmissibility_factor,
                            variantalpha_cross_protection_prob, 
                            variantalpha_imports_factor,
                            variantgamma_transmissibility,
                            variantgamma_transmissibility_factor,
                            variantgamma_cross_protection_prob, 
                            variantgamma_imports_factor,
                            variantgamma_severity_factor,
                            variantkappa_transmissibility,
                            variantkappa_transmissibility_factor,
                            variantkappa_cross_protection_prob, 
                            variantkappa_introduction_day,
                            variantkappa_severity_factor,
                            variantdelta_transmissibility, 
                            variantdelta_transmissibility_factor,
                            variantdelta_cross_protection_prob, 
                            variantdelta_imports_factor,
                            variantdelta_severity_factor,
                            imports_factor, days, seed, primary_cases_file, school_closure_policy,
                            school_closure_duration, school_closure_day, school_student_teacher_ratio,
                            shelter_in_place_delay_mean,
                            school_contacts,
                            classroom_contacts,
                            school_contact_factor,
                            workplace_contacts, office_contacts, workplace_contact_factor,
                            neighborhood_contacts, neighborhood_contact_factor,
                            enable_holiday_contacts,
                            holiday_contact_rate,
                            enable_vaccination
                            ) %>%
    mutate(job_id = sprintf("%s_%d", basename_jobs, row_number()), run_id = row_number(),
           params_file = sprintf('%s_%d.txt',basename_params,row_number()),
           reps = reps_per_job,
           state_code = state_code,
           advance_seeding = advance_seeding,
           epidemic_offset = epidemic_offset
           )

write.csv(report_scalars, file.path(output.dir, 'FRED_parameters.csv'), row.names= F, quote = F)

##===============================================##
## submit to CRC---------------
##===============================================##
fred_results_dir = file.path(output.dir,"FRED_RESULTS")

Sys.setenv(FRED_RESULTS=fred_results_dir)

if(!dir.exists(fred_results_dir)){
    dir.create(fred_results_dir)
}

# submit_jobs(experiment_supername_in = sprintf('FRED_CALIB-ASYMP-%.2fFM%.2f-KSUS%.2f-M-V%.0f-VX-%d',asymp_infectivity_in,face_mask_transmission_efficacy_in,kids_susceptibility_age_in, variants_in, vaccination_in),
#             experiment_name_in = as.character(state_code),
#             experiment_dir_in = output.dir,
#             params_base = basename_params,
#             job_base = basename_jobs,
#             reps = reps_per_job,
#             scalars = report_scalars,
#             FUN = write_cmd_function,cores_in=5,walltime_in = "3:00:00",
#             subsys="UGE",
#             fred_home_dir_in=fred_home, fred_results_in=fred_results_dir)