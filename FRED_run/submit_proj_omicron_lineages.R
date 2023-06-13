
##==============================================#
## Author: Guido Espana
## Simulate COVID-19 in BOG
## Year: 2019
##==============================================#
## User's input----------------
##==============================================#
#set.seed(123456)
library(pomp)
library(lubridate)
library(tidyverse)
library(fredtools)
## [x] Change 150 to 100
## [x] Change early delta to +15
## [x] Change start of 150 to august 1st

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
## Custom functions-------------
##==============================================#
write_cmd_function <- function(scalars_in, tmpfile){
    ## Delete jobs!!!
    cat("Deleting jobs\n")
    if(dir.exists(file.path(Sys.getenv('FRED_RESULTS'),'JOB'))){
        for(s in 1:nrow(scalars_in)){
            system(sprintf("fred_delete -f -k %s",scalars_in$job_id[s]), intern = T)
        }
        cat("All jobs deleted\n")
    }else{
        cat("FRED RESULTS doesn't exist, nothing to delete",Sys.getenv('FRED_RESULTS'),"\n")
    }
    job_cmd_str = sprintf("fred_job -f -k %s -I %d -p %s -n %.0f; Rscript ./post_process_fred_projection_lineages.R %s %.0f",
                          scalars_in$job_id,
                          scalars_in$run_id,
                          scalars_in$params_file,
                          scalars_in$reps,
                          scalars_in$job_id,
                          scalars_in$reps
                          )
    fileConn<-file(sprintf("run_files/%s",tmpfile))
    writeLines(job_cmd_str, fileConn)
    close(fileConn)
}

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
## Set STATE-------------------
##==============================================#
projection_label = 'omicron_lineages_import_test_3'
state_code = 11001
reps = 5
reps_per_job = 1
forecast_date = "2023-06-30"
fit_date = "2023-01-10"
asymp_infectivity_in = 1.0
face_mask_transmission_efficacy_in = 0.73 
kids_susceptibility_age_in = 10
variants_in = 1
vaccination_in = 70
subm_jobs = FALSE
args = (commandArgs(TRUE))
if(length(args) >= 1){
    projection_label = args[1]
    if(length(args) >= 2){
        state_name_input = args[2]
        if(length(args) >= 3){
            reps = as.numeric(args[3])
            if(length(args) >= 4){
                reps_per_job = as.numeric(args[4])
                if(length(args) >=5){
                    forecast_date = args[5]
                    if(length(args) >=6){
                        fit_date = args[6]
                        if(length(args) >=7){
                            asymp_infectivity_in = as.numeric(args[7])
                            if(length(args) >=8){
                                face_mask_transmission_efficacy_in = as.numeric(args[8])
                                if(length(args) >=9){
                                    kids_susceptibility_age_in = as.numeric(args[9])
                                    if(length(args) >=10){
                                        variants_in = as.numeric(args[10])
                                        if(length(args) >=11){
                                            vaccination_in = as.numeric(args[11])
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
forecast_date = as.Date(forecast_date)
fit_date = as.Date(fit_date)

##==============================================#
## FRED setup-------------
##==============================================#
fred_home = Sys.getenv('FRED_HOME')
scratch_dir = Sys.getenv('scratch_dir')

fred_defaults = sprintf("%s/input_files/defaults", fred_home)

output.dir = file.path(scratch_dir, sprintf('FRED_%.0f_projections_asymp_%.2f_fm_%.2f_ksus_%.2f_var_%.0f_vax_%03d_mov_%s',state_code, asymp_infectivity_in, face_mask_transmission_efficacy_in, kids_susceptibility_age_in, variants_in, vaccination_in, projection_label))

fred_results_dir = file.path(output.dir,"FRED_RESULTS")
Sys.setenv(FRED_RESULTS=fred_results_dir)
print(sprintf("Checking if %s exists", output.dir))
if(file.exists(output.dir)){
    print(sprintf("File %s exists", output.dir))
    ##system(paste('/bin/rm -rf ', output.dir,sep = ''))
    unlink(output.dir,recursive = T)
}
system(paste('/bin/mkdir -p ', output.dir, sep = ''))

file.copy('../scripts/post_process_fred_projection_lineages.R',output.dir)
# file.copy('../input_files/params_covid.txt','./input_files/params_covid.txt', overwrite = T)
# file.copy('../input_files/11001_schools_open_gps.csv','./input_files/11001_schools_open_gps.csv', overwrite = T)
file.copy("./input_files/infection_hospitalization_risk.csv", output.dir)
file.copy("./input_files/infection_hospitalization_risk_5.csv", output.dir)
file.copy('./input_files/Localidad_Unidad_Catastral.csv',output.dir, overwrite = T)
file.copy('./input_files/COL_hosp_duration_time.txt',output.dir,overwrite = T)

##==============================================#
## Sweep fixed parameters-------------------
##==============================================#
start_date = '2020-01-01' # YYYY-MM-DD
numDays = as.numeric(forecast_date - as.Date(start_date))
track_infection_events_in = 4 # No need to track infections in the calibrations right now
isolation_rate = 0.05175 # Specific for the US
report_place_of_infection_in = 1 # We might want to track place of infection
track_fatality_events_in = 1 
report_age_of_infection_in = 4 # Track individual ages of infections
report_incidence_by_county_in = 1 # We don't need this for the calibration either
shelter_in_place_students = 0 # Don't shelter students, use school closure for that
report_place_of_infection_in = 1 # We might want to track place of infection
synthetic_population_id = sprintf('synthetic_population_id = colombia_%d', state_code)

advance_seeding = 'exposed'
school_vacation_end = as.Date('2021-04-25')

## Shelter in place
interventions_st_df = read_csv('./input_files/interventions_Colombia.csv')
enable_shelter_in_place_in = 1

enable_shelter_in_place_timeseries_in = 1
shelter_in_place_delay_mean_in = as.integer(interventions_st_df$Shelter_in_place[interventions_st_df$State == state_code] - as.Date(start_date)) ## Just for book-keeping. The parameter should not have an effect in shelter timeseries

## School closure
early_school_vacation_end = as.Date('2021-04-25')
early_closure_duration = as.integer(early_school_vacation_end - interventions_st_df$School_closure[interventions_st_df$State == state_code])
early_school_reopening_day = as.numeric(early_school_vacation_end - as.Date(start_date))
current_open_date = as.Date('2020-10-15')

school_closure_policy_in = 'global_schedule'
school_closure_day_in = as.integer(interventions_st_df$School_closure[interventions_st_df$State == state_code] - as.Date(start_date))
school_closure_duration_in = as.integer(school_vacation_end - interventions_st_df$School_closure[interventions_st_df$State == state_code])

## Nursing home importations
enable_nursing_homes_importations_in = 1

## Age specific susceptibility
susceptibility_params = read_csv('./input_files/age_susceptibility_fit.csv')
enable_age_specific_susceptibility_in = 1
influenza_susceptibility_by_age_rate_in = susceptibility_params$rate_in[1]
influenza_susceptibility_by_age_cutoff_in = susceptibility_params$cutoff[1]

## Facemask usage
enable_face_mask_timeseries_in = 1
enable_face_mask_usage_in = 1
min_age_face_masks_in = 8

## Reduced capacity
enable_school_reduced_capacity_in = 0
school_reduced_capacity_in = 0.35
school_reduced_capacity_day_in =  as.integer(early_school_vacation_end - as.Date(start_date))
school_student_teacher_ratio_in = 18

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
## Filter and sample particles--------------
##==============================================#
state_code_in = state_code
## 1. Read parameters with LL
calibration_simdir = sprintf('FRED_%.0f_calibration_asymp_%.2f_fm_%.2f_ksus_%.2f_var_%.0f_vax_%03d_mov_old_files',state_code, asymp_infectivity_in,face_mask_transmission_efficacy_in,kids_susceptibility_age_in, variants_in, vaccination_in)

calibration_dir = file.path(getwd(), 'post_process/output','CALIBRATION',sprintf("%s_%s", calibration_simdir, "out"))
params_df = read_csv(file.path(calibration_dir, 'FRED_parameters_out.csv'))
fred_sweep_df = read_csv(file.path(calibration_dir, 'fred_output.csv'))
params_sweep_ll = params_df %>%
    filter(state_code == state_code_in) %>%
    mutate(LL_total = LL)


# ## 2. Sample based on their likelihood to ensure 10% sampled
# particles_w = 0.007
# prob_array = exp(-params_sweep_ll$LL_total*particles_w)
# if(variants_in == 0){
#     prob_array = exp(-(params_sweep_ll$LL_deaths*1 + params_sweep_ll$LL_age_time*1 + params_sweep_ll$LL_AR*0.005 + params_sweep_ll$LL_CF_age*0 + params_sweep_ll$LL_icu*0)*particles_w)

# }

# if(variants_in >= 1){
#     particles_w = 0.003
#     ##prob_array = exp(-(params_sweep_ll$LL_deaths*1 + params_sweep_ll$LL_CF_age*0 + params_sweep_ll$LL_icu*0 + params_sweep_ll$LL_dom * 2)*particles_w)
#     ##prob_array = exp(-(params_sweep_ll$LL_deaths + params_sweep_ll$LL_dom)*particles_w)
#     prob_array = exp(-(params_sweep_ll$LL_deaths)*particles_w)
# }

# if(sum(prob_array) > 0){
#     indx_sampled = sample.int(n = length(prob_array), size = reps, prob = prob_array, replace = T)
#     n_sampled = length(unique(indx_sampled))
#     print(params_sweep_ll$job_id[unique(indx_sampled)])
# }

## 7
#indx_sampled = c(13, 107, 192, 242, 394, 689, 760)

## 4
indx_sampled = c(192, 242, 394, 689, 760)

## 2
#indx_sampled = c(689, 760)

## 1
#indx_sampled = c(689)


n_sampled = length(indx_sampled)
print(params_sweep_ll$job_id[indx_sampled])

print(n_sampled)
scalars_sampled = params_sweep_ll[indx_sampled,] %>%
    dplyr::select(seed,
                  job_id,
                  imports_factor, influenza_transmissibility,
                  shelter_in_place_compliance, facemask_compliance,                  
                  influenza_susceptibility_by_age_offset,influenza_susceptibility_by_age_rate,
                  influenza_susceptibility_by_age_cutoff,influenza_susceptibility_by_age_high,
                  nursing_home_incidence_importations_factor,
                  workplace_contacts, office_contacts, workplace_contact_factor,
                  neighborhood_contacts, neighborhood_contact_factor,holiday_contact_rate,
                  community_contact_rate_1,
                  school_contacts, classroom_contacts,school_contact_factor,
                  neighborhood_same_age_bias,
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
                  variantdelta_severity_factor) 

# if(reps == 1){
#     scalars_sobol_df = data.frame('seed' = 139580060)
# }else{
#     scalars_sobol_df = sobol_design(
#         lower = c(seed=1),
#         upper = c(seed=as.integer(Sys.time())),
#         reps)
# }

#scalars_sobol_df = bind_cols(scalars_sobol_df, scalars_sampled)
scalars_sobol_df = scalars_sampled
## scalars_sobol_df$variantinf_transmissibility = scalars_sobol_df$variantinf_transmissibility * 0.9
##scalars_sobol_df$workplace_contact_factor = 1.0
##scalars_sobol_df$workplace_contacts = scalars_sobol_df$workplace_contacts * scalars_sobol_df$workplace_contact_factor
##scalars_sobol_df$office_contacts = scalars_sobol_df$office_contacts * scalars_sobol_df$workplace_contact_factor

scalars_sobol_df$variantalpha_susceptibility_by_age_offset  = scalars_sobol_df$influenza_susceptibility_by_age_offset 
scalars_sobol_df$variantalpha_susceptibility_by_age_rate = scalars_sobol_df$influenza_susceptibility_by_age_rate 
scalars_sobol_df$variantalpha_susceptibility_by_age_cutoff = scalars_sobol_df$influenza_susceptibility_by_age_cutoff
scalars_sobol_df$variantalpha_susceptibility_by_age_high = scalars_sobol_df$influenza_susceptibility_by_age_high

scalars_sobol_df$variantgamma_susceptibility_by_age_offset  = scalars_sobol_df$influenza_susceptibility_by_age_offset 
scalars_sobol_df$variantgamma_susceptibility_by_age_rate = scalars_sobol_df$influenza_susceptibility_by_age_rate 
scalars_sobol_df$variantgamma_susceptibility_by_age_cutoff = scalars_sobol_df$influenza_susceptibility_by_age_cutoff
scalars_sobol_df$variantgamma_susceptibility_by_age_high = scalars_sobol_df$influenza_susceptibility_by_age_high

scalars_sobol_df$variantkappa_susceptibility_by_age_offset  = scalars_sobol_df$influenza_susceptibility_by_age_offset 
scalars_sobol_df$variantkappa_susceptibility_by_age_rate = scalars_sobol_df$influenza_susceptibility_by_age_rate 
scalars_sobol_df$variantkappa_susceptibility_by_age_cutoff = scalars_sobol_df$influenza_susceptibility_by_age_cutoff
scalars_sobol_df$variantkappa_susceptibility_by_age_high = scalars_sobol_df$influenza_susceptibility_by_age_high

scalars_sobol_df$variantdelta_susceptibility_by_age_offset  = scalars_sobol_df$influenza_susceptibility_by_age_offset
scalars_sobol_df$variantdelta_susceptibility_by_age_rate = scalars_sobol_df$influenza_susceptibility_by_age_rate 
scalars_sobol_df$variantdelta_susceptibility_by_age_cutoff = scalars_sobol_df$influenza_susceptibility_by_age_cutoff
scalars_sobol_df$variantdelta_susceptibility_by_age_high = scalars_sobol_df$influenza_susceptibility_by_age_high

# scalars_sobol_df$variantdelta_transmissibility_factor = 2.1
# scalars_sobol_df$variantdelta_transmissibility = scalars_sobol_df$variantdelta_transmissibility_factor * scalars_sobol_df$influenza_transmissibility
#scalars_sobol_df$variantdelta_cross_protection_prob = 0.83 

### Omicron 
scalars_sobol_df$variantomicron_susceptibility_by_age_offset  = scalars_sobol_df$influenza_susceptibility_by_age_offset
scalars_sobol_df$variantomicron_susceptibility_by_age_rate = scalars_sobol_df$influenza_susceptibility_by_age_rate 
scalars_sobol_df$variantomicron_susceptibility_by_age_cutoff = scalars_sobol_df$influenza_susceptibility_by_age_cutoff
scalars_sobol_df$variantomicron_susceptibility_by_age_high = scalars_sobol_df$influenza_susceptibility_by_age_high

scalars_sobol_df$variantomicron_transmissibility_factor = 2.1
scalars_sobol_df$variantomicron_transmissibility = scalars_sobol_df$variantomicron_transmissibility_factor * scalars_sobol_df$influenza_transmissibility
scalars_sobol_df$variantomicron_cross_protection_prob = 0.7 
scalars_sobol_df$variantomicron_imports_factor = 2

### Omicron BAX
scalars_sobol_df$variantomicronBAX_susceptibility_by_age_offset  = scalars_sobol_df$influenza_susceptibility_by_age_offset
scalars_sobol_df$variantomicronBAX_susceptibility_by_age_rate = scalars_sobol_df$influenza_susceptibility_by_age_rate 
scalars_sobol_df$variantomicronBAX_susceptibility_by_age_cutoff = scalars_sobol_df$influenza_susceptibility_by_age_cutoff
scalars_sobol_df$variantomicronBAX_susceptibility_by_age_high = scalars_sobol_df$influenza_susceptibility_by_age_high

scalars_sobol_df$variantomicronBAX_transmissibility_factor = 2.1
scalars_sobol_df$variantomicronBAX_transmissibility = scalars_sobol_df$variantomicronBAX_transmissibility_factor * scalars_sobol_df$influenza_transmissibility
scalars_sobol_df$variantomicronBAX_cross_protection_prob = 0.7 
scalars_sobol_df$variantomicronBAX_imports_factor = 2

##scalars_sobol_df$community_contact_rate_1 = 0.8895
##scalars_sobol_df$neighborhood_contacts = 1.65

scalars_sobol_df$household_contact_factor = 1.0
scalars_sobol_df$household_contacts = household_contacts_in * scalars_sobol_df$household_contact_factor

scalars_sobol_df = bind_cols(scalars_sobol_df, get_ifr_param_string(scalars_sobol_df$variantgamma_severity_factor, 'gamma'))
scalars_sobol_df = bind_cols(scalars_sobol_df, get_ifr_param_string(scalars_sobol_df$variantkappa_severity_factor, 'kappa'))
scalars_sobol_df = bind_cols(scalars_sobol_df, get_ifr_param_string(scalars_sobol_df$variantdelta_severity_factor, 'delta'))

##==============================================#
## Initial conditions-------------------
##==============================================#
## initial_inf_file = sprintf('input_files/%d_imports_alternative.csv', state_code)
initial_inf_file = sprintf('input_files/%d_imports_combined.csv', state_code)
primary_cases_file_in = file.path(output.dir, sprintf('initial_cases_%d_%d.txt',state_code,1:reps))
#primary_cases_verylatedelta_file_in = file.path(output.dir, sprintf('initial_cases_verylatedelta_%d_%d.txt',state_code,1:reps))
#primary_cases_latedelta_file_in = file.path(output.dir, sprintf('initial_cases_latedelta_%d_%d.txt',state_code,1:reps))

# primary_cases_verylatedelta_omicron_file_in = file.path(output.dir, sprintf('initial_cases_verylatedelta_omicron_%d_%d.txt',state_code,1:reps))
# primary_cases_latedelta_omicron_file_in = file.path(output.dir, sprintf('initial_cases_latedelta_omicron_%d_%d.txt',state_code,1:reps))

primary_cases_plus_omicron_lineage_file_in = file.path(output.dir, sprintf('initial_cases_plus_omicron_lineage_%d_%d.txt',state_code,1:reps))
primary_cases_normal_omicron_lineage_file_in = file.path(output.dir, sprintf('initial_cases_normal_omicron_lineagen_%d_%d.txt',state_code,1:reps))

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
    tmp_df$Imports[tmp_df$Imports < 0] = 0
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
            mutate(Imports = round(0.5*scalars_sobol_df$variantalpha_imports_factor[nn] * TotalVariantImports)) %>%
            filter(Imports >= 1)

        gamma_imports = filter(variants_imp_df, variant == "20J (Gamma, V3)") %>%
            mutate(day = as.numeric(Date - as.Date('2020-01-01'))) %>%
            mutate(Imports = round(scalars_sobol_df$variantgamma_imports_factor[nn] * TotalVariantImports)) %>%
            filter(Imports >= 1,Date <= as.Date('2021-04-01'))

        ########################################################################################################
        delta_imports = filter(variants_imp_df, variant == "21A (Delta)") %>%
            mutate(day = as.numeric(Date - as.Date('2020-01-01'))) %>%
            mutate(Imports = round(scalars_sobol_df$variantdelta_imports_factor[nn] * TotalVariantImports)) %>%
            filter(Imports >= 1, Date > as.Date('2021-04-01')) %>%
            mutate(day = day + 15) %>%
            dplyr::select(day, Imports)
        
#         max_delta = delta_imports %>% filter(day == max(delta_imports$day))
#         max_delta = max_delta[rep(1:nrow(max_delta), numDays - max(max_delta$day)),]
#         max_delta$day = max_delta$day[1] + 1:nrow(max_delta)

#         delta_imports = bind_rows(delta_imports, max_delta)
        ########################################################################################################
        verylate_delta_imports = filter(variants_imp_df, variant == "21A (Delta)") %>%
            mutate(day = as.numeric(Date - as.Date('2020-01-01'))) %>%
            mutate(Imports = round(scalars_sobol_df$variantdelta_imports_factor[nn] * TotalVariantImports)) %>%
            filter(Imports >= 1, Date > as.Date('2021-04-01')) %>%
            mutate(day = day + 40) %>%
            dplyr::select(day, Imports)
        
        # max_delta = verylate_delta_imports %>% filter(day == max(verylate_delta_imports$day))
        # max_delta = max_delta[rep(1:nrow(max_delta), numDays - max(max_delta$day)),]
        # max_delta$day = max_delta$day[1] + 1:nrow(max_delta)
        # verylate_delta_imports = bind_rows(verylate_delta_imports, max_delta)
        ########################################################################################################
        late_delta_imports = filter(variants_imp_df, variant == "21A (Delta)") %>%
            mutate(day = as.numeric(Date - as.Date('2020-01-01'))) %>%
            mutate(Imports = round(scalars_sobol_df$variantdelta_imports_factor[nn] * TotalVariantImports)) %>%
            filter(Imports >= 1, Date > as.Date('2021-04-01')) %>%
            mutate(day = day + 30 + 15) %>%
            dplyr::select(day, Imports)
        
        # max_delta = late_delta_imports %>% filter(day == max(late_delta_imports$day))
        # max_delta = max_delta[rep(1:nrow(max_delta), numDays - max(max_delta$day)),]
        # max_delta$day = max_delta$day[1] + 1:nrow(max_delta)
        # late_delta_imports = bind_rows(late_delta_imports, max_delta)
        ########################################################################################################
        omicron_imports = filter(variants_imp_df, variant == "21K (Omicron)" | variant == "21L (Omicron)" | variant == "22A (Omicron)") %>%
            mutate(day = as.numeric(Date - as.Date('2020-01-01'))) %>%
            mutate(Imports = round(scalars_sobol_df$variantomicron_imports_factor[nn] * TotalVariantImports)) %>%
            filter(Imports >= 1, Date > as.Date('2021-11-01')) %>%
            mutate(day = day + 9) %>%
            dplyr::select(day, Imports)
                
        # max_omicron = omicron_imports %>% filter(day == max(omicron_imports$day))
        # max_omicron = max_omicron[rep(1:nrow(max_omicron), numDays - max(max_omicron$day)),]
        # max_omicron$day = max_omicron$day[1] + 1:nrow(max_omicron)
        # omicron_imports = bind_rows(omicron_imports, max_omicron)
        ########################################################################################################
        omicronBAX_imports = filter(variants_imp_df, variant == "22E (Omicron)") %>%
            mutate(day = as.numeric(Date - as.Date('2020-01-01'))) %>%
            mutate(Imports = round(scalars_sobol_df$variantomicronBAX_imports_factor[nn] * TotalVariantImports)) %>%
            #filter(Imports >= 1, Date > as.Date('2021-12-15'), Date <= as.Date('2022-03-20')) %>%
            filter(Imports >= 1, Date > as.Date('2022-09-01')) %>%
            mutate(day = day - 90) %>%
            dplyr::select(day, Imports)
        
        # max_omicron = omicronBAX_imports %>% filter(day == max(omicronBAX_imports$day))
        # max_omicron = max_omicron[rep(1:nrow(max_omicron), numDays - max(max_omicron$day)),]
        # max_omicron$day = max_omicron$day[1] + 1:nrow(max_omicron)
        # omicronBAX_imports = bind_rows(omicronBAX_imports, max_omicron)
        ########################################################################################################
        plus_omicronBAX_imports = filter(variants_imp_df, variant == "22E (Omicron)") %>%
            mutate(day = as.numeric(Date - as.Date('2020-01-01'))) %>%
            mutate(Imports = round(scalars_sobol_df$variantomicronBAX_imports_factor[nn] * TotalVariantImports)) %>%
            #filter(Imports >= 1, Date > as.Date('2021-12-15'), Date <= as.Date('2022-03-20')) %>%
            filter(Imports >= 1, Date > as.Date('2022-09-01')) %>%
            mutate(day = day - 75) %>%
            dplyr::select(day, Imports)
        
        # max_omicron = plus_omicronBAX_imports %>% filter(day == max(plus_omicronBAX_imports$day))
        # max_omicron = max_omicron[rep(1:nrow(max_omicron), numDays - max(max_omicron$day)),]
        # max_omicron$day = max_omicron$day[1] + 1:nrow(max_omicron)
        # omicronBAX_imports = bind_rows(plus_omicronBAX_imports, max_omicron)
        ########################################################################################################

        tmp_df$Date = tmp_df$day + as.Date('2020-01-01')
        tmp_df_main = tmp_df[tmp_df$Date < as.Date('2020-12-10'),]               

        tmp_df_variantkappa = tmp_df[tmp_df$Date >= as.Date('2020-01-01') + round(scalars_sobol_df$variantkappa_introduction_day[nn]) & tmp_df$Date <= as.Date('2021-06-01'),]
    
        ## Base
        init_cases_base_lines = sprintf('%.0f %.0f %.0f 0 1 %.0f', tmp_df_main$day, tmp_df_main$day, tmp_df_main$Imports, tmp_df_main$Imports)

        ## Alpha
        init_cases_alpha_lines = sprintf('%.0f %.0f %.0f 1 1 %.0f', alpha_imports$day, alpha_imports$day, alpha_imports$Imports, alpha_imports$Imports)

        ## Gamma
        init_cases_gamma_lines = sprintf('%.0f %.0f %.0f 2 1 %.0f', gamma_imports$day, gamma_imports$day, gamma_imports$Imports, gamma_imports$Imports)

        ## Mu
        init_cases_kappa_lines = sprintf('%.0f %.0f %.0f 3 1 %.0f', tmp_df_variantkappa$day, tmp_df_variantkappa$day, tmp_df_variantkappa$Imports, tmp_df_variantkappa$Imports)

        ## Delta
        init_cases_delta_lines = sprintf('%.0f %.0f %.0f 4 1 %.0f', delta_imports$day, delta_imports$day, delta_imports$Imports, delta_imports$Imports)
        init_cases_verylatedelta_lines = sprintf('%.0f %.0f %.0f 4 1 %.0f', verylate_delta_imports$day, verylate_delta_imports$day, verylate_delta_imports$Imports, verylate_delta_imports$Imports)
        init_cases_latedelta_lines = sprintf('%.0f %.0f %.0f 4 1 %.0f', late_delta_imports$day, late_delta_imports$day, late_delta_imports$Imports, late_delta_imports$Imports)
        
        ## Omicron
        init_cases_omicron_lines = sprintf('%.0f %.0f %.0f 5 1 %.0f', omicron_imports$day, omicron_imports$day, omicron_imports$Imports, omicron_imports$Imports)  
        
        ## OmicronBAX
        init_cases_omicronBAX_plus_lines = sprintf('%.0f %.0f %.0f 6 1 %.0f', plus_omicronBAX_imports$day, plus_omicronBAX_imports$day, plus_omicronBAX_imports$Imports, plus_omicronBAX_imports$Imports)  
        init_cases_omicronBAX_normal_lines = sprintf('%.0f %.0f %.0f 6 1 %.0f', omicronBAX_imports$day, omicronBAX_imports$day, omicronBAX_imports$Imports, omicronBAX_imports$Imports)  
 
        init_cases_lines = c(init_cases_base_lines, init_cases_alpha_lines, init_cases_gamma_lines, init_cases_kappa_lines, init_cases_delta_lines, init_cases_omicron_lines, init_cases_omicronBAX_plus_lines)
        init_cases_plus_lines = c(init_cases_base_lines, init_cases_alpha_lines, init_cases_gamma_lines, init_cases_kappa_lines, init_cases_verylatedelta_lines, init_cases_omicron_lines, init_cases_omicronBAX_plus_lines)   
        init_cases_normal_lines = c(init_cases_base_lines, init_cases_alpha_lines, init_cases_gamma_lines, init_cases_kappa_lines, init_cases_verylatedelta_lines, init_cases_omicron_lines, init_cases_omicronBAX_normal_lines)   
        
        fileConn<-file(primary_cases_plus_omicron_lineage_file_in[nn])
        writeLines(init_cases_plus_lines, fileConn)
        close(fileConn)
        
        fileConn<-file(primary_cases_normal_omicron_lineage_file_in[nn])
        writeLines(init_cases_normal_lines, fileConn)
        close(fileConn)
        
        # init_cases_lines = c(init_cases_base_lines, init_cases_alpha_lines, init_cases_gamma_lines, init_cases_kappa_lines, init_cases_delta_lines, init_cases_omicron_lines, init_cases_omicronBAX_plus_lines)
        # init_cases_late_lines = c(init_cases_base_lines, init_cases_alpha_lines, init_cases_gamma_lines,init_cases_kappa_lines, init_cases_latedelta_lines, init_cases_omicron_lines, init_cases_omicronBAX_normal_lines)
        # init_cases_verylate_lines = c(init_cases_base_lines, init_cases_alpha_lines, init_cases_gamma_lines, init_cases_kappa_lines, init_cases_verylatedelta_lines, init_cases_omicron_lines, init_cases_omicronBAX_normal_lines)
        
#         fileConn<-file(primary_cases_verylatedelta_omicron_file_in[nn])
#         writeLines(init_cases_verylate_lines, fileConn)
#         close(fileConn)
        
#         fileConn<-file(primary_cases_latedelta_omicron_file_in[nn])
#         writeLines(init_cases_late_lines, fileConn)
#         close(fileConn)
    }else{
        init_cases_lines = sprintf('%.0f %.0f %.0f 0 1 %.0f', tmp_df$day, tmp_df$day, tmp_df$Imports, tmp_df$Imports)
    }    
    
    fileConn<-file(primary_cases_file_in[nn])
    writeLines(init_cases_lines, fileConn)
    close(fileConn)
    
}

##==============================================#
## Shelter in place-------------------
##==============================================#
post_lockdown_mobility_in = 0.1
shelter_time_file = file.path(output.dir, sprintf('shelter_timeseries_%s_%d.txt',state_code,1:reps))
shelter_time_lockdown_file = file.path(output.dir, sprintf('shelter_timeseries_lockdown_%s_%d.txt',state_code,1:reps))

data_mov_day = as.integer(as.Date(fit_date) - as.Date(start_date))

## Moving away from google's data to Grandata census-tract specific
## shelter_timeseries_df = read_csv('./input_files/interventions_covid_timevarying_shelter.csv')
shelter_timeseries_df = read_csv('./input_files/11001_mobility_trends.csv')
shelter_timeseries_df$day = as.numeric(difftime(shelter_timeseries_df$date, as.Date(start_date), units='days'))

shelter_timeseries_lockdown_df = read_csv('./input_files/11001_mobility_trends_lockdown.csv')
shelter_timeseries_lockdown_df$day = as.numeric(difftime(shelter_timeseries_lockdown_df$date, as.Date(start_date), units='days'))

#for(nn in indx_sampled){
for(nn in 1:reps){
   
    ## For now, choose the mean
    tmp_df = shelter_timeseries_df %>% filter(State == state_code, replicate == 1) %>%
        mutate(shelter_trend = shelter_trend * scalars_sobol_df$shelter_in_place_compliance[nn])
    max_tmp_day = max(tmp_df$day)

    tmp_low_df = shelter_timeseries_lockdown_df %>% filter(State == state_code, replicate == 1) %>%
        mutate(shelter_trend = shelter_trend * scalars_sobol_df$shelter_in_place_compliance[nn])

    
    if(variants_in == 0){
        ## if variants are not enabled, push mobility to the max
        tmp_df$shelter_trend[tmp_df$date > as.Date('2021-02-15')] = 0.0
    }
    
    
    tmp_df$shelter_trend[tmp_df$shelter_trend < 0.0001] = 0.0
    tmp_df$shelter_trend[tmp_df$shelter_trend > 1.0] = 1.0

    tmp_low_df$shelter_trend[tmp_low_df$shelter_trend < 0.0001] = 0.0
    tmp_low_df$shelter_trend[tmp_low_df$shelter_trend > 1.0] = 1.0   
    
    tmp_df_tail = tmp_df[tmp_df$day == max_tmp_day,]
    tmp_df_tail$day = max_tmp_day + 1
    tmp_df_tail$shelter_trend = post_lockdown_mobility_in * tmp_df_tail$shelter_trend    

    tmp_low_df_tail = tmp_low_df[tmp_low_df$day == max_tmp_day,]
    tmp_low_df_tail$day = max_tmp_day + 1
    tmp_low_df_tail$shelter_trend = post_lockdown_mobility_in * tmp_low_df_tail$shelter_trend    
    
    shelter_time_lines = sprintf('%.0f %.0f %.4f 0 120 %d%s', tmp_df$day, tmp_df$day, tmp_df$shelter_trend, tmp_df$State, tmp_df$SCACODIGO)
    
    shelter_time_lockdown_lines = sprintf('%.0f %.0f %.4f 0 120 %d%s', tmp_low_df$day, tmp_low_df$day, tmp_low_df$shelter_trend, tmp_low_df$State, tmp_low_df$SCACODIGO)

    numDays_ = 731
    if(numDays > max(tmp_df_tail$day)){
        shelter_time_tail_lines = sprintf('%.0f %.0f %.4f 0 120 %d%s', tmp_df_tail$day, numDays_, tmp_df_tail$shelter_trend, tmp_df_tail$State, tmp_df_tail$SCACODIGO)
        shelter_time_tail_low_lines = sprintf('%.0f %.0f %.4f 0 120 %d%s', tmp_low_df_tail$day, numDays_, tmp_low_df_tail$shelter_trend, tmp_low_df_tail$State, tmp_low_df_tail$SCACODIGO)
        shelter_time_lines = c(shelter_time_lines,shelter_time_tail_lines)
        shelter_time_lockdown_lines = c(shelter_time_lockdown_lines,shelter_time_tail_low_lines)
    }
    ## Add shelter by age
    ## shelter_age_lines = sprintf("82 212 %.4f 60 120", scalars_sobol_df$shelter_in_place_age_compliance[nn])
    fileConn<-file(shelter_time_file[nn])
    writeLines(shelter_time_lines, fileConn)
    close(fileConn)

    fileConn<-file(shelter_time_lockdown_file[nn])
    writeLines(shelter_time_lockdown_lines, fileConn)
    close(fileConn)        
}

##==============================================#
## holiday increase-------------------
##==============================================#
if(variants_in == 0){
    community_timeseries_df = read_csv('./input_files/interventions_covid_timevarying_community.csv')
}else{
    community_timeseries_df = read_csv('./input_files/interventions_covid_timevarying_community_baseline.csv')
}
community_timeseries_file = file.path(output.dir, sprintf("community_timeseries_%d.txt", 1:reps))
community_timeseries_high_file = file.path(output.dir, sprintf("community_timeseries_high_%d.txt", 1:reps))
community_timeseries_high_20_file = file.path(output.dir, sprintf("community_timeseries_high_20_%d.txt", 1:reps))
community_timeseries_high_40_file = file.path(output.dir, sprintf("community_timeseries_high_40_%d.txt", 1:reps))
community_timeseries_high_60_file = file.path(output.dir, sprintf("community_timeseries_high_60_%d.txt", 1:reps))

community_timeseries_low_file = file.path(output.dir, sprintf("community_timeseries_low_%d.txt", 1:reps))
community_timeseries_df$day = as.numeric(community_timeseries_df$date - as.Date(start_date))
community_timeseries_df$community_trend[community_timeseries_df$date < as.Date('2020-05-15')] = community_timeseries_df$community_trend[community_timeseries_df$date == as.Date('2020-05-15')]

# ##holiweek
# community_timeseries_df$community_trend[community_timeseries_df$date > as.Date('2021-03-27') & community_timeseries_df$date <= as.Date('2021-04-03')] = community_timeseries_df$community_trend[community_timeseries_df$date > as.Date('2021-03-27') & community_timeseries_df$date <= as.Date('2021-04-03')] * 1.1

community_timeseries_df$community_trend = community_timeseries_df$community_trend - community_timeseries_df$community_trend[1]
#for(nn in indx_sampled){
for(nn in 1:reps){
    print(scalars_sobol_df$community_contact_rate_1[nn])
    print(scalars_sobol_df$job_id[nn])
    tmp_comm_df = community_timeseries_df    
    tmp_comm_df$contact_rate = tmp_comm_df$community_trend * scalars_sobol_df$community_contact_rate_1[nn] + 1
    
    tmp_comm_df$contact_rate[tmp_comm_df$contact_rate < 0] = 0
    
    high_comm_df = tmp_comm_df
    high_comm_df$contact_rate[high_comm_df$date >= as.Date('2021-12-18')] = high_comm_df$contact_rate[high_comm_df$date >= as.Date('2021-12-18')] * 1.1
    
    high_20_comm_df = tmp_comm_df
    high_20_comm_df$contact_rate[high_20_comm_df$date >= as.Date('2021-12-18')] = high_20_comm_df$contact_rate[high_20_comm_df$date >= as.Date('2021-12-18')] * 1.2
    
    high_40_comm_df = tmp_comm_df
    high_40_comm_df$contact_rate[high_40_comm_df$date >= as.Date('2021-12-18')] = high_40_comm_df$contact_rate[high_40_comm_df$date >= as.Date('2021-12-18')] * 1.4
    
    high_60_comm_df = tmp_comm_df
    high_60_comm_df$contact_rate[high_60_comm_df$date >= as.Date('2021-12-18')] = high_60_comm_df$contact_rate[high_60_comm_df$date >= as.Date('2021-12-18')] * 1.6

    low_comm_df = tmp_comm_df
    low_comm_df$contact_rate[low_comm_df$date >= as.Date('2021-12-18')] = low_comm_df$contact_rate[low_comm_df$date >= as.Date('2021-12-18')] * 0.5

    community_time_lines = sprintf('%.0f %.0f %.4f', tmp_comm_df$day, tmp_comm_df$day, tmp_comm_df$contact_rate)
    community_time_high_lines = sprintf('%.0f %.0f %.4f', high_comm_df$day, high_comm_df$day, high_comm_df$contact_rate)
    ########################################
    community_time_high_20_lines = sprintf('%.0f %.0f %.4f', high_20_comm_df$day, high_20_comm_df$day, high_20_comm_df$contact_rate)
    community_time_high_40_lines = sprintf('%.0f %.0f %.4f', high_40_comm_df$day, high_40_comm_df$day, high_40_comm_df$contact_rate)
    community_time_high_60_lines = sprintf('%.0f %.0f %.4f', high_60_comm_df$day, high_60_comm_df$day, high_60_comm_df$contact_rate)
    #######################################
    community_time_low_lines = sprintf('%.0f %.0f %.4f', low_comm_df$day, low_comm_df$day, low_comm_df$contact_rate)

    fileConn<-file(community_timeseries_file[nn])
    writeLines(community_time_lines, fileConn)
    close(fileConn)

    fileConn<-file(community_timeseries_high_file[nn])
    writeLines(community_time_high_lines, fileConn)
    close(fileConn)
    ######################################
    fileConn<-file(community_timeseries_high_20_file[nn])
    writeLines(community_time_high_20_lines, fileConn)
    close(fileConn)
    
    fileConn<-file(community_timeseries_high_40_file[nn])
    writeLines(community_time_high_40_lines, fileConn)
    close(fileConn)
    
    fileConn<-file(community_timeseries_high_60_file[nn])
    writeLines(community_time_high_60_lines, fileConn)
    close(fileConn)
    ######################################
    fileConn<-file(community_timeseries_low_file[nn])
    writeLines(community_time_low_lines, fileConn)
    close(fileConn)
}

##==============================================#
## Vaccine timeseries-------------------
##==============================================#
vaccine_daily_capacity_file = file.path(output.dir, "vaccination_daily_capacity_timeseries.txt")
vaccine_stock_file = file.path(output.dir, "vaccination_stock_timeseries.txt")

vaccine_daily_high_capacity_file = file.path(output.dir, "vaccination_daily_high_capacity_timeseries.txt")
vaccine_stock_high_capacity_file = file.path(output.dir, "vaccination_stock_timeseries_high.txt")

vaccine_stock_df = read.csv('./input_files/11001_vaccine_stock_timeseries.csv')
vaccine_stock_df$day = as.numeric(as.Date(vaccine_stock_df$Date) - as.Date(start_date))
vaccine_daily_df = read.csv('./input_files/11001_vaccine_capacity_timeseries.csv')
vaccine_daily_df$day = as.numeric(as.Date(vaccine_daily_df$Date) - as.Date(start_date))

vaccine_stock_high_df = vaccine_stock_df
vaccine_stock_high_df$ID[vaccine_stock_high_df$ID == 2 & vaccine_stock_high_df$Date == max(vaccine_stock_high_df$Date)] = 6

vaccine_daily_high_df = vaccine_daily_df
vaccine_daily_high_df$VaccinesApplied[vaccine_daily_high_df$Date > as.Date('2021-10-24')] = 60000

vaccine_stock_lines = sprintf("%.0f %.0f %.0f %.0f",vaccine_stock_df$day, vaccine_stock_df$day,vaccine_stock_df$TotalVaccines,vaccine_stock_df$ID)

fileConn<-file(vaccine_stock_file)
writeLines(vaccine_stock_lines, fileConn)
close(fileConn)

vaccine_stock_high_lines = sprintf("%.0f %.0f %.0f %.0f",vaccine_stock_high_df$day, vaccine_stock_high_df$day,vaccine_stock_high_df$TotalVaccines,vaccine_stock_high_df$ID)

fileConn<-file(vaccine_stock_high_capacity_file)
writeLines(vaccine_stock_high_lines, fileConn)
close(fileConn)

vaccine_capacity_lines = c(sprintf("%.0f %.0f",vaccine_daily_df$day,vaccine_daily_df$VaccinesApplied))

fileConn<-file(vaccine_daily_capacity_file)
writeLines(vaccine_capacity_lines, fileConn)
close(fileConn)

vaccine_high_capacity_lines = c(sprintf("%.0f %.0f",vaccine_daily_high_df$day,vaccine_daily_high_df$VaccinesApplied))

fileConn<-file(vaccine_daily_high_capacity_file)
writeLines(vaccine_high_capacity_lines, fileConn)
close(fileConn)

##==============================================#
## Facemask compliance-------------------
##==============================================#
facemask_compliance_file = file.path(output.dir, sprintf('facemask_compliance_%d_%d.txt',state_code,1:reps))
facemask_compliance_off_file = file.path(output.dir, sprintf('facemask_compliance_%d_%d.txt',state_code,1:reps))
facemask_community_df  =  read_csv('./input_files/facemask_timeseries_compliance_community.csv')  %>%
    mutate(state_name = tolower(state_name))

facemask_timeseries_df = read_csv('./input_files/facemask_timeseries_compliance.csv') %>%
    dplyr::select(-Day)
facemask_timeseries_df$day = as.numeric(difftime(facemask_timeseries_df$Date, as.Date(start_date), units='days'))
facemask_community_df$day = as.numeric(difftime(facemask_community_df$Date, as.Date(start_date), units='days'))

school_facemask_compliance = 0.75
community_scaling = 1.0

for(nn in 1:reps){
    ## FM compliance stays at the last value in the data
    tmp_df = facemask_timeseries_df %>% filter(State == tolower(interventions_st_df$state_name[interventions_st_df$State == state_code])) %>%
        mutate(FacemaskTrends = FacemaskTrends * scalars_sobol_df$facemask_compliance[nn]) %>%
        dplyr::select(day, FacemaskTrends)

    tmp_comm_df = filter(facemask_community_df, state_name == tolower(interventions_st_df$state_name[interventions_st_df$State == state_code])) %>%
        mutate(FacemaskTrends = FacemaskTrends * community_scaling) %>%
        dplyr::select(day, FacemaskTrends) %>%
        filter(day > max(tmp_df$day)) %>%
        mutate(FacemaskTrends = ifelse(FacemaskTrends < 1.0,FacemaskTrends, 1.0))

    tmp_comm_df = bind_rows(tmp_df, tmp_comm_df)
    tmp_comm_off = tmp_comm_df
    
    if(max(tmp_df$day) < numDays){
        tmp_df2 = data.frame(day = seq(from=max(tmp_df$day) + 1, to = numDays),
                             FacemaskTrends = tail(tmp_df$FacemaskTrends,1), stringsAsFactors = F)
        tmp_df = bind_rows(tmp_df, tmp_df2)
    }

    off_day = as.Date('2022-05-01') - as.Date(start_date)
    tmp_off_df = tmp_df
    tmp_off_df$FacemaskTrends[tmp_off_df$day > off_day] = 0
    tmp_comm_off$FacemaskTrends[tmp_comm_off$day > off_day] = 0
    
    ## Schools start facemask compliance when they open
    facemask_school_day = as.integer(current_open_date - as.Date(start_date))
    compliance_time_lines = c(sprintf('%.0f %.0f %.4f %s', facemask_school_day, facemask_school_day, school_facemask_compliance, "school"),
                              sprintf('%.0f %.0f %.4f %s', facemask_school_day, facemask_school_day, school_facemask_compliance, "school"),
                              sprintf('%.0f %.0f %.4f %s', facemask_school_day, facemask_school_day, school_facemask_compliance, "classroom"))
    
    compliance_off_time_lines = c(sprintf('%.0f %.0f %.4f %s', facemask_school_day, facemask_school_day, school_facemask_compliance, "school"),
                                  sprintf('%.0f %.0f %.4f %s', facemask_school_day, facemask_school_day, school_facemask_compliance, "school"),
                                  sprintf('%.0f %.0f %.4f %s', facemask_school_day, facemask_school_day, school_facemask_compliance, "classroom"))
    
    for(loc in c("workplace", "office", "other")){
        if(!(loc %in% c("other_community"))){ ## Says other community just to ignore the community facemask for now go with other
            compliance_time_lines = c(compliance_time_lines,
                                      sprintf('%.0f %.0f %.4f %s', tmp_df$day, tmp_df$day,  tmp_df$FacemaskTrends, loc)
                                      )
            compliance_off_time_lines = c(compliance_off_time_lines,
                                      sprintf('%.0f %.0f %.4f %s', tmp_off_df$day, tmp_off_df$day,  tmp_off_df$FacemaskTrends, loc)
                                      )
        }else{
            compliance_time_lines = c(compliance_time_lines,
                                      sprintf('%.0f %.0f %.4f %s', tmp_comm_df$day, tmp_comm_df$day,  tmp_comm_df$FacemaskTrends, loc)
                                      )
            compliance_off_time_lines = c(compliance_off_time_lines,
                                      sprintf('%.0f %.0f %.4f %s', tmp_comm_off$day, tmp_comm_off$day,  tmp_comm_off$FacemaskTrends, loc)
                                      )
        }
    }

    
    ## Schools start facemask compliance when they open
    # facemask_school_day = as.integer(current_open_date - as.Date(start_date))
    # school_high_lines = c()
    # for(loc in c("school", "classroom")){
    #     school_high_lines = c(school_high_lines, sprintf('%.0f %.0f %.4f %s', facemask_school_day, numDays, school_facemask_compliance, loc))
    # }    
    # school_high_lines = c(compliance_time_lines, school_high_lines)
    # school_high_off_lines = c(compliance_off_time_lines, school_high_lines)
    
    school_high_lines = compliance_time_lines
    school_high_off_lines = compliance_off_time_lines
    
    ## Schools start facemask compliance when they open    
    fileConn<-file(facemask_compliance_file[nn])
    writeLines(school_high_lines, fileConn)
    close(fileConn)

    fileConn<-file(facemask_compliance_off_file[nn])
    writeLines(school_high_off_lines, fileConn)
    close(fileConn)        
}

##==============================================#
## School closure-------------------
##==============================================#
schools_open_list = read_csv('./input_files/11001_schools_open_gps.csv')
schools_open_list$PropOpen[schools_open_list$PropOpen > 1] = 1.0
current_open_date = as.Date('2020-10-15')
localidad_esc = read_csv('./input_files/Localidad_Unidad_Catastral.csv')
localidad_list = 1:19

school_schedule_closed_file = file.path(output.dir, sprintf('school_schedule_closed_%d.txt',state_code))
school_schedule_open_all_file = file.path(output.dir, sprintf('school_schedule_open_all_%d.txt',state_code))

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

forecast_reopen_day = as.numeric(as.Date('2021-06-07') - as.Date(start_date))
school_reopen_to_date = read_csv(school_reopen_list_file) %>%
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
    mutate(end_day = ifelse(end_day == max(end_day), forecast_reopen_day, end_day)) %>%
    replace_na(list(zipcode = "")) %>%
    filter(zipcode != "Bogota")

## Write files
{
    initial_closure = sprintf("%d %d 1 20 0", school_closure_day_in,  as.integer(as.Date('2020-12-15') - as.Date(start_date)))
    ## Make sure this is working properly
    vacation_closure = sprintf("%d %d 1 20 0", as.integer(as.Date('2020-12-15') - as.Date(start_date)),numDays)
    
    ##school_baseline_reopen_date = as.integer(early_school_vacation_end - as.Date(start_date)) + 1
    school_current_early_open_date = as.integer(current_open_date - as.Date(start_date)) + 1
    ## Schools already open
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
    school_open_to_date_lines = current_open_school_lines
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

    ## 2021 reopened schools to date
    tmp_school_open_to_date = filter(school_reopen_to_date, Grade != 'UNIVERSITY') 
    tmp_college_open_to_date = filter(school_reopen_to_date, Grade == 'UNIVERSITY')
    school_open_to_date_lines = c(school_open_to_date_lines,
                                  sprintf("%d %d %d %d %.4f 0 11001%s",
                                          tmp_school_open_to_date$start_day,
                                          tmp_school_open_to_date$end_day,
                                          tmp_school_open_to_date$MinAge,
                                          tmp_school_open_to_date$MaxAge,
                                          tmp_school_open_to_date$Capacity,
                                          tmp_school_open_to_date$zipcode))
    
    school_open_to_date_lines = c(school_open_to_date_lines,
                                  sprintf("%d %d 18 20 %.4f",
                                          tmp_college_open_to_date$start_day,
                                          tmp_college_open_to_date$end_day, tmp_college_open_to_date$Capacity))
    
    ## School open at full capacity on February
    school_open_all_lines = c(   
        school_open_to_date_lines,
        sprintf("%d %d 1 17 1.0",  as.integer(as.Date('2022-02-28') - as.Date(start_date)) + 1, numDays),
        sprintf("%d %d 18 20 1.0",  as.integer(as.Date('2022-02-28') - as.Date(start_date)) + 1, numDays))
        
    fileConn<-file(school_schedule_open_all_file)
    writeLines(school_open_all_lines, fileConn)
    close(fileConn)
}

##==============================================#
## Intervention parameters-------------------
##==============================================#
epidemic_offset = 0
num_demes_in = 1
intervention_base = data.frame(
    enable_shelter_in_place = 1,
    enable_shelter_in_place_timeseries = 1,
    shelter_in_place_delay_mean = as.integer(interventions_st_df$Shelter_in_place[interventions_st_df$State == state_code] - as.Date(start_date)),    
    school_closure_policy = school_closure_policy_in,
    school_closure_day = as.integer(interventions_st_df$School_closure[interventions_st_df$State == state_code] - as.Date(start_date)),
    school_closure_duration = as.integer(school_vacation_end - interventions_st_df$School_closure[interventions_st_df$State == state_code]),
    school_student_teacher_ratio = school_student_teacher_ratio_in,
    stringsAsFactors = F)

intervention_df = intervention_base[rep(1,reps),]

intervention_df$enable_community_contact_timeseries = enable_community_contact_timeseries_in

## Vaccination
intervention_df$enable_transmission_bias = enable_transmission_bias_in
intervention_df$enable_vaccination = 0
intervention_df$vaccination_capacity_file = vaccine_daily_capacity_file
intervention_df$vaccine_stock_timeseries_file = vaccine_stock_file

intervention_df$vaccination_phases_names = "12 age ltc essentialworkers age essentialworkers comorbidity teachers age age age age age"
intervention_df$vaccination_phases_age_low = "12 80 16 16 60 16 16 16 50 40 30 20 16"
intervention_df$vaccination_phases_age_high = "12 120 120 120 79 59 59 59 59 49 39 29 19"
intervention_df$vaccination_phases_id = "12 1 1 1 2 2 3 3 4 5 6 7 8"
intervention_df$vaccination_phases_pop_prob = "12 0.0 0.0 0.0001 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0"

##HOSPITALIZATION DURATION VARYING IN TIME
hosp_duration_time_file_in = file.path(output.dir, 'COL_hosp_duration_time.txt')
intervention_df$enable_hospitalization_duration_timeseries = 0
intervention_df$hospitalization_duration_timeseries_file = hosp_duration_time_file_in

intervention_df_default = intervention_df %>%    
    mutate(seed = floor(scalars_sobol_df$seed),
           primary_cases_file = primary_cases_file_in,
           imports_factor = scalars_sobol_df$imports_factor,
           shelter_in_place_file = shelter_time_file,
           shelter_in_place_compliance = scalars_sobol_df$shelter_in_place_compliance,
           influenza_transmissibility = scalars_sobol_df$influenza_transmissibility,
           school_closure_duration = school_closure_duration_in,
           influenza_susceptibility_by_age_offset = scalars_sobol_df$influenza_susceptibility_by_age_offset,
           influenza_susceptibility_by_age_rate  = scalars_sobol_df$influenza_susceptibility_by_age_rate,
           influenza_susceptibility_by_age_cutoff = scalars_sobol_df$influenza_susceptibility_by_age_cutoff,
           influenza_susceptibility_by_age_high = scalars_sobol_df$influenza_susceptibility_by_age_high,
           nursing_home_incidence_importations_factor = scalars_sobol_df$nursing_home_incidence_importations_factor,
           
           face_mask_timeseries_file = facemask_compliance_off_file,
           facemask_compliance = scalars_sobol_df$facemask_compliance,
           
           enable_school_reduced_capacity = enable_school_reduced_capacity_in,
           school_reduced_capacity = school_reduced_capacity_in,
           school_reduced_capacity_day = school_reduced_capacity_day_in,
           enable_face_mask_usage = enable_face_mask_usage_in,
           school_contacts = scalars_sobol_df$school_contacts,
           classroom_contacts = scalars_sobol_df$classroom_contacts,
           school_contact_factor = scalars_sobol_df$school_contact_factor,
           neighborhood_contacts = scalars_sobol_df$neighborhood_contacts,
           neighborhood_contact_factor = scalars_sobol_df$neighborhood_contact_factor,
           holiday_contact_rate = scalars_sobol_df$holiday_contact_rate,
           workplace_contacts = scalars_sobol_df$workplace_contacts,
           household_contacts = scalars_sobol_df$household_contacts,
           office_contacts = scalars_sobol_df$office_contacts,
           workplace_contact_factor = scalars_sobol_df$workplace_contact_factor,
           household_contact_factor = scalars_sobol_df$household_contact_factor,
           school_global_schedule_file = school_schedule_closed_file,
           enable_holiday_contacts = enable_holiday_contacts_in,
           community_contact_rate_1 = scalars_sobol_df$community_contact_rate_1,
           community_contact_timeseries_file = community_timeseries_file,        
           holiday_start = holiday_start_in,
           holiday_end = holiday_end_in,
           neighborhood_same_age_bias = scalars_sobol_df$neighborhood_same_age_bias,
           
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
           variantgamma_severity_factor = scalars_sobol_df$variantgamma_severity_factor,
           
           variantkappa_transmissibility =  scalars_sobol_df$variantkappa_transmissibility,
           variantkappa_transmissibility_factor =  scalars_sobol_df$variantkappa_transmissibility_factor,
           variantkappa_cross_protection_prob = scalars_sobol_df$variantkappa_cross_protection_prob,
           variantkappa_introduction_day = scalars_sobol_df$variantkappa_introduction_day,
           variantkappa_susceptibility_by_age_offset  = scalars_sobol_df$variantkappa_susceptibility_by_age_offset,
           variantkappa_susceptibility_by_age_rate = scalars_sobol_df$variantkappa_susceptibility_by_age_rate,
           variantkappa_susceptibility_by_age_cutoff = scalars_sobol_df$variantkappa_susceptibility_by_age_cutoff,
           variantkappa_susceptibility_by_age_high = scalars_sobol_df$variantkappa_susceptibility_by_age_high,
           variantkappa_severity_factor = scalars_sobol_df$variantkappa_severity_factor,

           variantdelta_transmissibility =  1.5*scalars_sobol_df$variantdelta_transmissibility,
           variantdelta_transmissibility_factor =  scalars_sobol_df$variantdelta_transmissibility_factor,
           variantdelta_cross_protection_prob = scalars_sobol_df$variantdelta_cross_protection_prob,
           variantdelta_imports_factor = scalars_sobol_df$variantdelta_imports_factor,
           variantdelta_susceptibility_by_age_offset  = scalars_sobol_df$variantdelta_susceptibility_by_age_offset,
           variantdelta_susceptibility_by_age_rate = scalars_sobol_df$variantdelta_susceptibility_by_age_rate,
           variantdelta_susceptibility_by_age_cutoff = scalars_sobol_df$variantdelta_susceptibility_by_age_cutoff,
           variantdelta_susceptibility_by_age_high = scalars_sobol_df$variantdelta_susceptibility_by_age_high,   
           variantdelta_severity_factor = scalars_sobol_df$variantdelta_severity_factor,

           variantomicron_transmissibility =  scalars_sobol_df$variantomicron_transmissibility,
           variantomicron_transmissibility_factor =  scalars_sobol_df$variantomicron_transmissibility_factor,
           variantomicron_cross_protection_prob = scalars_sobol_df$variantomicron_cross_protection_prob,
           variantomicron_imports_factor = scalars_sobol_df$variantomicron_imports_factor,
           variantomicron_susceptibility_by_age_offset  = scalars_sobol_df$variantomicron_susceptibility_by_age_offset,
           variantomicron_susceptibility_by_age_rate = scalars_sobol_df$variantomicron_susceptibility_by_age_rate,
           variantomicron_susceptibility_by_age_cutoff = scalars_sobol_df$variantomicron_susceptibility_by_age_cutoff,
           variantomicron_susceptibility_by_age_high = scalars_sobol_df$variantomicron_susceptibility_by_age_high,   
           
           variantomicronBAX_transmissibility =  scalars_sobol_df$variantomicronBAX_transmissibility,
           variantomicronBAX_transmissibility_factor =  scalars_sobol_df$variantomicronBAX_transmissibility_factor,
           variantomicronBAX_cross_protection_prob = scalars_sobol_df$variantomicronBAX_cross_protection_prob,
           variantomicronBAX_imports_factor = scalars_sobol_df$variantomicronBAX_imports_factor,
           variantomicronBAX_susceptibility_by_age_offset  = scalars_sobol_df$variantomicronBAX_susceptibility_by_age_offset,
           variantomicronBAX_susceptibility_by_age_rate = scalars_sobol_df$variantomicronBAX_susceptibility_by_age_rate,
           variantomicronBAX_susceptibility_by_age_cutoff = scalars_sobol_df$variantomicronBAX_susceptibility_by_age_cutoff,
           variantomicronBAX_susceptibility_by_age_high = scalars_sobol_df$variantomicronBAX_susceptibility_by_age_high,   

           intervention_id = "default")

if(variants_in == 0){
    intervention_df_default$diseases = 1
    intervention_df_default$enable_disease_cross_protection = 0
    intervention_df_default$influenza_cross_protection_prob = 1.0
    intervention_df_default$disease_names = "1 influenza"
}

if(variants_in == 1){
    intervention_df_default$diseases = 7
    intervention_df_default$enable_disease_cross_protection = 1
    intervention_df_default$influenza_cross_protection_prob = 1.0
    intervention_df_default$disease_names = "7 influenza variantalpha variantgamma variantkappa variantdelta variantomicron variantomicronBAX"
}


intervention_vax_df = intervention_df_default %>%
    mutate(intervention_id = "default_vax",
           shelter_in_place_file = shelter_time_file,
           enable_vaccination = 1
           )

scalars_intervention = bind_rows(intervention_df_default,
                                 intervention_vax_df)

vaccine_list = c(1)

## VACCINE
intervention_vac_priority_df = data.frame(
    vaccination_phases_names = "18 essentialworkers ltc age age age age teachers age age age age age age age age age age age ",
    vaccination_phases_age_low = "18 18 18 80 70 60 55 18 50 45 40 35 30 25 20 15 12 10 3",
    vaccination_phases_age_high = "18 120 120 120 79 69 59 120 54 49 44 39 34 29 24 19 14 11 9",
    vaccination_phases_id = "18 1 2 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17",
    vaccination_phases_pop_prob = "18 0.0001 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0", stringsAsFactors = F)

intervention_vac_noprior_df = data.frame(
    vaccination_phases_names = "10 age ltc essentialworkers age essentialworkers comorbidity teachers age age age",
    vaccination_phases_age_low = "10 80 16 16 60 16 16 16 50 40 16",
    vaccination_phases_age_high = "10 120 120 120 79 59 59 59 59 49 39",
    vaccination_phases_id = "10 1 1 1 2 2 3 3 4 5 6",
    vaccination_phases_pop_prob = "10 0.0 0.0 0.0001 0.0 0.0 0.0 0.0 0.0 0.0 0.0", stringsAsFactors = F)

# ## SCENARIOS
# lockdown_file_list = list(
#     shelter_time_file, shelter_time_file,
#     shelter_time_file,shelter_time_file,
#     shelter_time_file,shelter_time_file,
#     shelter_time_file,shelter_time_file,
#     shelter_time_file,shelter_time_file,
#     shelter_time_file,shelter_time_file,
#     shelter_time_file,shelter_time_file)

# community_file_list = list(
#     community_timeseries_file,     community_timeseries_high_file,
#     community_timeseries_high_file,community_timeseries_high_file,
#     community_timeseries_high_file,community_timeseries_high_file,
#     community_timeseries_high_file,community_timeseries_high_file,
#     community_timeseries_file,     community_timeseries_file,
#     community_timeseries_file,     community_timeseries_file,
#     community_timeseries_file,     community_timeseries_file)

# school_file_list = list(
#     school_schedule_closed_file, school_schedule_open_all_file,
#     school_schedule_open_all_file,school_schedule_open_all_file,
#     school_schedule_open_all_file,school_schedule_open_all_file,
#     school_schedule_open_all_file,school_schedule_open_all_file,
#     school_schedule_closed_file,school_schedule_closed_file,
#     school_schedule_closed_file,school_schedule_closed_file,
#     school_schedule_closed_file,school_schedule_closed_file)

# vaccination_capacity_file_list = list(
#     vaccine_daily_capacity_file, vaccine_daily_capacity_file,
#     vaccine_daily_high_capacity_file,vaccine_daily_high_capacity_file,
#     vaccine_daily_high_capacity_file,vaccine_daily_high_capacity_file,
#     vaccine_daily_capacity_file, vaccine_daily_capacity_file,
#     vaccine_daily_high_capacity_file,vaccine_daily_high_capacity_file,
#     vaccine_daily_high_capacity_file,vaccine_daily_high_capacity_file,
#     vaccine_daily_capacity_file, vaccine_daily_capacity_file)

# vaccine_stock_timeseries_file_list = list(
#     vaccine_stock_file,               vaccine_stock_file,
#     vaccine_stock_file,               vaccine_stock_file,
#     vaccine_stock_high_capacity_file, vaccine_stock_high_capacity_file,
#     vaccine_stock_file,               vaccine_stock_high_capacity_file,
#     vaccine_stock_file,               vaccine_stock_file,
#     vaccine_stock_high_capacity_file, vaccine_stock_high_capacity_file,
#     vaccine_stock_file,               vaccine_stock_high_capacity_file)

# vaccination_priority_list = list(
#     intervention_vac_priority_df, intervention_vac_priority_df,
#     intervention_vac_noprior_df,  intervention_vac_priority_df,
#     intervention_vac_noprior_df,  intervention_vac_priority_df,
#     intervention_vac_noprior_df,  intervention_vac_priority_df,
#     intervention_vac_noprior_df,  intervention_vac_priority_df,
#     intervention_vac_noprior_df,  intervention_vac_priority_df,
#     intervention_vac_noprior_df,  intervention_vac_priority_df)


# lockdown_list = c('base','highcontacts',
#                   'highvac','highvacpriority',
#                   'highvaclong','highvaclongpriority',
#                   'highcontactsnopriority', 'highcontactslong',
#                   'basehighvac','basehighvacpriority',
#                   'basehighvaclong','basehighvaclongpriority',
#                   'basehighcontactsnopriority', 'basehighcontactslong')


# lockdown_file_list = list(
#     shelter_time_file, shelter_time_file,
#     shelter_time_file, shelter_time_file, shelter_time_file,
#     shelter_time_file, shelter_time_file, shelter_time_file)

# community_file_list = list(
#     community_timeseries_file,      community_timeseries_file,
#     community_timeseries_high_20_file, community_timeseries_high_40_file, community_timeseries_high_60_file,
#     community_timeseries_high_20_file, community_timeseries_high_40_file, community_timeseries_high_60_file)

# school_file_list = list(
#     school_schedule_closed_file,  school_schedule_closed_file,
#     school_schedule_open_all_file,school_schedule_open_all_file,school_schedule_open_all_file,
#     school_schedule_open_all_file,school_schedule_open_all_file,school_schedule_open_all_file)

# vaccination_capacity_file_list = list(
#     vaccine_daily_capacity_file,      vaccine_daily_high_capacity_file,
#     vaccine_daily_capacity_file,      vaccine_daily_capacity_file, vaccine_daily_capacity_file,
#     vaccine_daily_high_capacity_file, vaccine_daily_high_capacity_file, vaccine_daily_high_capacity_file)

# vaccine_stock_timeseries_file_list = list(
#     vaccine_stock_file,               vaccine_stock_high_capacity_file,
#     vaccine_stock_file,               vaccine_stock_file,               vaccine_stock_file,
#     vaccine_stock_high_capacity_file, vaccine_stock_high_capacity_file, vaccine_stock_high_capacity_file)

# vaccination_priority_list = list(
#     intervention_vac_priority_df, intervention_vac_priority_df,
#     intervention_vac_noprior_df,  intervention_vac_noprior_df,  intervention_vac_noprior_df,
#     intervention_vac_priority_df,  intervention_vac_priority_df, intervention_vac_priority_df)

# lockdown_list = c('base','highvac',
#                   'highcontacts_20', 'highcontacts_40', 'highcontacts_60',
#                   'highvac_contacts_20','highvac_contacts_40', 'highvac_contacts_60')

lockdown_file_list = list(
    shelter_time_file, 
    shelter_time_file)

community_file_list = list(
    community_timeseries_file, 
    community_timeseries_file)

school_file_list = list(
    school_schedule_closed_file,  
    school_schedule_open_all_file)

vaccination_capacity_file_list = list(
    vaccine_daily_capacity_file,
    vaccine_daily_capacity_file
    )

vaccine_stock_timeseries_file_list = list(
    vaccine_stock_file,               
    vaccine_stock_file)

vaccination_priority_list = list(
    intervention_vac_priority_df, 
    intervention_vac_priority_df)

#lockdown_list = c('base','open_school')
lockdown_list = c('open_school')

for(ll in 1:length(lockdown_list)){
    intervention_tmp_df = intervention_df_default %>%
        mutate(
            shelter_in_place_file = lockdown_file_list[[ll]],
            school_global_schedule_file = school_file_list[[ll]],
            community_contact_timeseries_file = community_file_list[[ll]],
            vaccine_stock_timeseries_file = vaccine_stock_timeseries_file_list[[ll]],
            vaccination_capacity_file = vaccination_capacity_file_list[[ll]])
    intervention_tmp_df[,colnames(vaccination_priority_list[[ll]])] = vaccination_priority_list[[ll]]

    for(v in 1:length(vaccine_list)){
        vv = vaccine_list[v]

        intervention_default_df = intervention_tmp_df %>%
            mutate(intervention_id = sprintf("default-vax_%d-mov_%s", vv, lockdown_list[ll]),                   
                   enable_vaccination = vv
                   )
        
        #############################################################################################         
        #####################################  Omicron   ############################################
        #############################################################################################      
        ## Works with oscillation
#  plus_omicronBAX2_cross_060_transexcess_0_50_sev_100_omicron_cross_040_transexcess_0_40_severity_400_df =  intervention_tmp_df %>%
#           mutate(intervention_id = sprintf("plus_omicronBAX2_cross_040_transexcess_0_40_sev_100_omicron_cross_020_transexcess_0_50_severity_400-vax_%d-mov_%s", vv, lockdown_list[ll]),
#                  primary_cases_file = primary_cases_plus_omicron_lineage_file_in,
#                  enable_vaccination = vv,
#                  vaccination_phases_enable_discrete_timing = 1,
#                  variantdelta_severity_factor = 0.025*variantdelta_severity_factor,
                 
#                  variantomicron_transmissibility = 0.60*variantdelta_transmissibility, ## Trans Excess 0.9, 1.2, 1.5, 2
#                  variantomicron_cross_protection_prob = 0.40, ## Escape 0.5, 0.8 (Cross 0.5, 0.2)
#                  variantomicron_severity_factor = 4.0*variantdelta_severity_factor,
#                  get_ifr_param_string(variantomicron_severity_factor, 'omicron'),
                 
#                  variantomicronBAX_severity_factor = variantomicron_severity_factor,
#                  variantomicronBAX_cross_protection_prob = 0.6,
#                  variantomicronBAX_transmissibility = 0.5*variantomicron_transmissibility,
#                  get_ifr_param_string(variantomicronBAX_severity_factor, 'omicronBAX')
#                  )  

#  plus_omicronBAX2_cross_050_transexcess_0_70_sev_070_omicron_cross_030_transexcess_0_58_severity_500_df =  intervention_tmp_df %>%
#           mutate(intervention_id = sprintf("plus_omicronBAX2_cross_050_transexcess_0_70_sev_070_omicron_cross_030_transexcess_0_58_severity_500-vax_%d-mov_%s", vv, lockdown_list[ll]),
#                  primary_cases_file = primary_cases_plus_omicron_lineage_file_in,
#                  enable_vaccination = vv,
#                  vaccination_phases_enable_discrete_timing = 1,
#                  variantdelta_severity_factor = 0.025*variantdelta_severity_factor,
                 
#                  variantomicron_transmissibility = 0.58*variantdelta_transmissibility, ## Trans Excess 0.9, 1.2, 1.5, 2
#                  variantomicron_cross_protection_prob = 0.30, ## Escape 0.5, 0.8 (Cross 0.5, 0.2)
#                  variantomicron_severity_factor = 5.0*variantdelta_severity_factor,
#                  get_ifr_param_string(variantomicron_severity_factor, 'omicron'),
                 
#                  variantomicronBAX_severity_factor = 0.7*variantomicron_severity_factor,
#                  variantomicronBAX_cross_protection_prob = 0.5,
#                  variantomicronBAX_transmissibility = 0.7*variantomicron_transmissibility,
#                  get_ifr_param_string(variantomicronBAX_severity_factor, 'omicronBAX')
#                  )  
        
#  plus_omicronBAX2_cross_050_transexcess_0_70_sev_070_omicron_cross_030_transexcess_0_58_severity_500_df =  intervention_tmp_df %>%
#           mutate(intervention_id = sprintf("plus_omicronBAX2_cross_048_transexcess_0_75_sev_070_omicron_cross_030_transexcess_0_58_severity_500-vax_%d-mov_%s", vv, lockdown_list[ll]),
#                  primary_cases_file = primary_cases_plus_omicron_lineage_file_in,
#                  enable_vaccination = vv,
#                  vaccination_phases_enable_discrete_timing = 1,
#                  variantdelta_severity_factor = 0.025*variantdelta_severity_factor,
                 
#                  variantomicron_transmissibility = 0.58*variantdelta_transmissibility, ## Trans Excess 0.9, 1.2, 1.5, 2
#                  variantomicron_cross_protection_prob = 0.30, ## Escape 0.5, 0.8 (Cross 0.5, 0.2)
#                  variantomicron_severity_factor = 5.0*variantdelta_severity_factor,
#                  get_ifr_param_string(variantomicron_severity_factor, 'omicron'),
                 
#                  variantomicronBAX_severity_factor = 0.7*variantomicron_severity_factor,
#                  variantomicronBAX_cross_protection_prob = 0.48,
#                  variantomicronBAX_transmissibility = 0.75*variantomicron_transmissibility,
#                  get_ifr_param_string(variantomicronBAX_severity_factor, 'omicronBAX')
#                  )  
#   plus_omicronBAX2_cross_048_transexcess_0_75_sev_050_omicron_cross_030_transexcess_0_58_severity_500_df =  intervention_tmp_df %>%
#           mutate(intervention_id = sprintf("plus_omicronBAX2_cross_048_transexcess_0_75_sev_050_omicron_cross_030_transexcess_0_58_severity_500-vax_%d-mov_%s", vv, lockdown_list[ll]),
#                  primary_cases_file = primary_cases_plus_omicron_lineage_file_in,
#                  enable_vaccination = vv,
#                  vaccination_phases_enable_discrete_timing = 1,
#                  variantdelta_severity_factor = 0.025*variantdelta_severity_factor,
                 
#                  variantomicron_transmissibility = 0.58*variantdelta_transmissibility, ## Trans Excess 0.9, 1.2, 1.5, 2
#                  variantomicron_cross_protection_prob = 0.30, ## Escape 0.5, 0.8 (Cross 0.5, 0.2)
#                  variantomicron_severity_factor = 5.0*variantdelta_severity_factor,
#                  get_ifr_param_string(variantomicron_severity_factor, 'omicron'),
                 
#                  variantomicronBAX_severity_factor = 0.5*variantomicron_severity_factor,
#                  variantomicronBAX_cross_protection_prob = 0.48,
#                  variantomicronBAX_transmissibility = 0.75*variantomicron_transmissibility,
#                  get_ifr_param_string(variantomicronBAX_severity_factor, 'omicronBAX')
#                  )  
        
 normal_omicronBAX2_cross_050_transexcess_0_85_sev_060_omicron_cross_030_transexcess_0_5_severity_420_df =  intervention_tmp_df %>%
          mutate(intervention_id = sprintf("normal_omicronBAX2_cross_050_transexcess_0_85_sev_060_omicron_cross_030_transexcess_0_5_severity_420-vax_%d-mov_%s", vv, lockdown_list[ll]),
                 primary_cases_file = primary_cases_normal_omicron_lineage_file_in,
                 enable_vaccination = vv,
                 vaccination_phases_enable_discrete_timing = 1,
                 variantdelta_severity_factor = 0.025*variantdelta_severity_factor,
                 
                 variantomicron_transmissibility = 0.5*variantdelta_transmissibility, ## Trans Excess 0.9, 1.2, 1.5, 2
                 variantomicron_cross_protection_prob = 0.30, ## Escape 0.5, 0.8 (Cross 0.5, 0.2)
                 variantomicron_severity_factor = 4.2*variantdelta_severity_factor,
                 get_ifr_param_string(variantomicron_severity_factor, 'omicron'),
                 
                 variantomicronBAX_severity_factor = 0.6*variantomicron_severity_factor,
                 variantomicronBAX_cross_protection_prob = 0.5,
                 variantomicronBAX_transmissibility = 0.85*variantomicron_transmissibility,
                 get_ifr_param_string(variantomicronBAX_severity_factor, 'omicronBAX')
                 )  
        
 plus_omicronBAX2_cross_050_transexcess_0_85_sev_060_omicron_cross_030_transexcess_0_5_severity_420_df =  intervention_tmp_df %>%
          mutate(intervention_id = sprintf("plus_omicronBAX2_cross_050_transexcess_0_85_sev_060_omicron_cross_030_transexcess_0_5_severity_420-vax_%d-mov_%s", vv, lockdown_list[ll]),
                 primary_cases_file = primary_cases_plus_omicron_lineage_file_in,
                 enable_vaccination = vv,
                 vaccination_phases_enable_discrete_timing = 1,
                 variantdelta_severity_factor = 0.025*variantdelta_severity_factor,
                 
                 variantomicron_transmissibility = 0.5*variantdelta_transmissibility, ## Trans Excess 0.9, 1.2, 1.5, 2
                 variantomicron_cross_protection_prob = 0.30, ## Escape 0.5, 0.8 (Cross 0.5, 0.2)
                 variantomicron_severity_factor = 4.2*variantdelta_severity_factor,
                 get_ifr_param_string(variantomicron_severity_factor, 'omicron'),
                 
                 variantomicronBAX_severity_factor = 0.6*variantomicron_severity_factor,
                 variantomicronBAX_cross_protection_prob = 0.5,
                 variantomicronBAX_transmissibility = 0.85*variantomicron_transmissibility,
                 get_ifr_param_string(variantomicronBAX_severity_factor, 'omicronBAX')
                 )  
        
 plus_omicronBAX2_cross_042_transexcess_0_85_sev_060_omicron_cross_030_transexcess_0_5_severity_420_df =  intervention_tmp_df %>%
          mutate(intervention_id = sprintf("plus_omicronBAX2_cross_042_transexcess_0_85_sev_060_omicron_cross_030_transexcess_0_5_severity_420-vax_%d-mov_%s", vv, lockdown_list[ll]),
                 primary_cases_file = primary_cases_plus_omicron_lineage_file_in,
                 enable_vaccination = vv,
                 vaccination_phases_enable_discrete_timing = 1,
                 variantdelta_severity_factor = 0.025*variantdelta_severity_factor,
                 
                 variantomicron_transmissibility = 0.5*variantdelta_transmissibility, ## Trans Excess 0.9, 1.2, 1.5, 2
                 variantomicron_cross_protection_prob = 0.30, ## Escape 0.5, 0.8 (Cross 0.5, 0.2)
                 variantomicron_severity_factor = 4.2*variantdelta_severity_factor,
                 get_ifr_param_string(variantomicron_severity_factor, 'omicron'),
                 
                 variantomicronBAX_severity_factor = 0.6*variantomicron_severity_factor,
                 variantomicronBAX_cross_protection_prob = 0.42,
                 variantomicronBAX_transmissibility = 0.85*variantomicron_transmissibility,
                 get_ifr_param_string(variantomicronBAX_severity_factor, 'omicronBAX')
                 )  

 normal_omicronBAX2_cross_042_transexcess_0_85_sev_060_omicron_cross_030_transexcess_0_5_severity_420_df =  intervention_tmp_df %>%
          mutate(intervention_id = sprintf("normal_omicronBAX2_cross_042_transexcess_0_85_sev_060_omicron_cross_030_transexcess_0_5_severity_420-vax_%d-mov_%s", vv, lockdown_list[ll]),
                 primary_cases_file = primary_cases_normal_omicron_lineage_file_in,
                 enable_vaccination = vv,
                 vaccination_phases_enable_discrete_timing = 1,
                 variantdelta_severity_factor = 0.025*variantdelta_severity_factor,
                 
                 variantomicron_transmissibility = 0.5*variantdelta_transmissibility, ## Trans Excess 0.9, 1.2, 1.5, 2
                 variantomicron_cross_protection_prob = 0.30, ## Escape 0.5, 0.8 (Cross 0.5, 0.2)
                 variantomicron_severity_factor = 4.2*variantdelta_severity_factor,
                 get_ifr_param_string(variantomicron_severity_factor, 'omicron'),
                 
                 variantomicronBAX_severity_factor = 0.6*variantomicron_severity_factor,
                 variantomicronBAX_cross_protection_prob = 0.42,
                 variantomicronBAX_transmissibility = 0.85*variantomicron_transmissibility,
                 get_ifr_param_string(variantomicronBAX_severity_factor, 'omicronBAX')
                 ) 
        
 normal_omicronBAX2_cross_046_transexcess_0_85_sev_060_omicron_cross_030_transexcess_0_5_severity_420_df =  intervention_tmp_df %>%
          mutate(intervention_id = sprintf("normal_omicronBAX2_cross_046_transexcess_0_85_sev_060_omicron_cross_030_transexcess_0_5_severity_420-vax_%d-mov_%s", vv, lockdown_list[ll]),
                 primary_cases_file = primary_cases_normal_omicron_lineage_file_in,
                 enable_vaccination = vv,
                 vaccination_phases_enable_discrete_timing = 1,
                 variantdelta_severity_factor = 0.025*variantdelta_severity_factor,
                 
                 variantomicron_transmissibility = 0.5*variantdelta_transmissibility, ## Trans Excess 0.9, 1.2, 1.5, 2
                 variantomicron_cross_protection_prob = 0.30, ## Escape 0.5, 0.8 (Cross 0.5, 0.2)
                 variantomicron_severity_factor = 4.2*variantdelta_severity_factor,
                 get_ifr_param_string(variantomicron_severity_factor, 'omicron'),
                 
                 variantomicronBAX_severity_factor = 0.6*variantomicron_severity_factor,
                 variantomicronBAX_cross_protection_prob = 0.46,
                 variantomicronBAX_transmissibility = 0.85*variantomicron_transmissibility,
                 get_ifr_param_string(variantomicronBAX_severity_factor, 'omicronBAX')
                 ) 


        #intervention_default_df,
        #scalars_intervention
        scalars_intervention = bind_rows(normal_omicronBAX2_cross_050_transexcess_0_85_sev_060_omicron_cross_030_transexcess_0_5_severity_420_df,
                                         plus_omicronBAX2_cross_050_transexcess_0_85_sev_060_omicron_cross_030_transexcess_0_5_severity_420_df,
                                         plus_omicronBAX2_cross_042_transexcess_0_85_sev_060_omicron_cross_030_transexcess_0_5_severity_420_df,
                                         normal_omicronBAX2_cross_042_transexcess_0_85_sev_060_omicron_cross_030_transexcess_0_5_severity_420_df,
                                         normal_omicronBAX2_cross_046_transexcess_0_85_sev_060_omicron_cross_030_transexcess_0_5_severity_420_df
                                        )                                
    }
}

scalars_intervention = bind_cols(scalars_intervention, get_ifr_param_string(scalars_intervention$variantgamma_severity_factor, 'gamma'))
scalars_intervention = bind_cols(scalars_intervention, get_ifr_param_string(scalars_intervention$variantkappa_severity_factor, 'kappa'))
scalars_intervention = bind_cols(scalars_intervention, get_ifr_param_string(scalars_intervention$variantdelta_severity_factor, 'delta'))

##==============================================#
## Create parameters to sweep-----------------
##==============================================#
## 1. Create DF with parameters
enable_age_specific_susceptibility_min_in = 0
influenza_susceptibility_by_age_minage_in = 10
influenza_susceptibility_by_age_minvalue_in = 1

## if(kids_susceptibility_in < 1.0){
##     enable_age_specific_susceptibility_min_in = 1
## }

scalars = scalars_intervention %>%
    mutate(days = numDays,
           track_infection_events = track_infection_events_in,
           enable_age_specific_susceptibility = enable_age_specific_susceptibility_in,
           enable_age_specific_susceptibility_min = enable_age_specific_susceptibility_min_in,
           influenza_susceptibility_by_age_minage = influenza_susceptibility_by_age_minage_in,
           influenza_susceptibility_by_age_minvalue = influenza_susceptibility_by_age_minvalue_in,
           influenza_asymp_infectivity = asymp_infectivity_in,
           influenza_face_mask_transmission_efficacy = face_mask_transmission_efficacy_in,
           enable_nursing_homes_importations = enable_nursing_homes_importations_in,
           track_fatality_events = track_fatality_events_in,
           report_age_of_infection = report_age_of_infection_in,
           report_incidence_by_county = report_incidence_by_county_in,
           report_place_of_infection = report_place_of_infection_in,
           isolation_rate = isolation_rate,
           shelter_in_place_students = shelter_in_place_students,
           num_demes = num_demes_in,
           synthetic_population_id = synthetic_population_id,
           advance_seeding = advance_seeding,
           epidemic_offset = epidemic_offset,
           enable_shelter_in_place = enable_shelter_in_place,
           enable_shelter_in_place_timeseries = enable_shelter_in_place_timeseries,
           shelter_in_place_delay_mean = shelter_in_place_delay_mean,           
           enable_face_mask_timeseries = enable_face_mask_timeseries_in,
           min_age_face_masks = min_age_face_masks_in,
           
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
           
           variantomicron_asymp_infectivity = asymp_infectivity_in,
           variantomicron_susceptibility_by_age_minage = influenza_susceptibility_by_age_minage_in,
           variantomicron_susceptibility_by_age_minvalue = influenza_susceptibility_by_age_minvalue_in,
           variantomicron_face_mask_transmission_efficacy = face_mask_transmission_efficacy_in,
           
           variantomicronBAX_asymp_infectivity = asymp_infectivity_in,
           variantomicronBAX_susceptibility_by_age_minage = influenza_susceptibility_by_age_minage_in,
           variantomicronBAX_susceptibility_by_age_minvalue = influenza_susceptibility_by_age_minvalue_in,
           variantomicronBAX_face_mask_transmission_efficacy = face_mask_transmission_efficacy_in,
           start_date = start_date)


##===============================================##
## Write the parameters to files---------------
##===============================================##
#defaults_params = './input_files/params_covid.txt'
defaults_covid_params = './input_files/params_covid.txt'
defaults_alpha_params = './input_files/params_covid_alpha.txt'
defaults_gamma_params = './input_files/params_covid_gamma.txt'
defaults_kappa_params = './input_files/params_covid_kappa.txt'
defaults_delta_params = './input_files/params_covid_delta.txt'
defaults_omicron_params = './input_files/params_covid_omicron.txt'
defaults_omicronBAX_params = './input_files/params_covid_omicronBAX.txt'
defaults_vaccine_params = './input_files/params_covid_vaccine.txt'

if(variants_in >= 1){
    if(vaccination_in > 0){
        defaults_vaccine_params = sprintf('./input_files/params_covid_vaccine_%d.txt', vaccination_in)
    }
    defaults_params = './input_files/params_covid_combined.txt'
    system(sprintf("/bin/cat %s %s %s %s %s %s %s %s > %s",
                   defaults_covid_params,
                   defaults_alpha_params, 
                   defaults_gamma_params,
                   defaults_kappa_params,
                   defaults_delta_params, 
                   defaults_omicron_params,
                   defaults_omicronBAX_params,
                   defaults_vaccine_params,
                   defaults_params, 
                   intern = T))
}

basename_params = sprintf('covid_%d_params',state_code)
basename_jobs = sprintf('FRED_%d_projections_asymp',state_code)
write_fred_parameters(scalars, defaults_params, output.dir,basename.in=basename_params, fred_defaults = fred_defaults)

## print report scalars parameters file with IDs
report_scalars = dplyr::select(
                            scalars, influenza_transmissibility, influenza_asymp_infectivity,
                            influenza_face_mask_transmission_efficacy,
                            enable_age_specific_susceptibility_min,
                            influenza_susceptibility_by_age_minage,
                            influenza_susceptibility_by_age_minvalue,
                            shelter_in_place_compliance, start_date,
                            imports_factor, days, seed, 
                            primary_cases_file, school_closure_policy, school_student_teacher_ratio,
                            enable_face_mask_usage, enable_face_mask_timeseries, enable_community_contact_timeseries,
                            facemask_compliance, min_age_face_masks,
                            enable_school_reduced_capacity, school_reduced_capacity, school_reduced_capacity_day,
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
                            variantomicron_transmissibility,
                            variantomicron_transmissibility_factor,
                            variantomicron_cross_protection_prob, 
                            variantomicron_imports_factor,
                            variantomicron_severity_factor,
                            variantomicronBAX_transmissibility,
                            variantomicronBAX_transmissibility_factor,
                            variantomicronBAX_cross_protection_prob, 
                            variantomicronBAX_imports_factor,
                            variantomicronBAX_severity_factor,
                            school_closure_duration, school_closure_day, shelter_in_place_delay_mean,
                            school_contacts,school_contact_factor,
                            nursing_home_incidence_importations_factor,
                            classroom_contacts,
                            neighborhood_contacts, neighborhood_contact_factor,
                            workplace_contacts, office_contacts, workplace_contact_factor, household_contact_factor, household_contacts,
                            enable_holiday_contacts,
                            community_contact_rate_1,
                            holiday_contact_rate,
                            holiday_start,
                            holiday_end,
                            enable_vaccination,
                            neighborhood_same_age_bias,
                            intervention_id) %>%
    mutate(job_id = sprintf("%s_%d", basename_jobs, row_number()), run_id = row_number(),
           params_file = sprintf('%s_%d.txt',basename_params,row_number()),
           reps = reps_per_job,
           state_code = state_code,
           advance_seeding = advance_seeding,
           epidemic_offset = epidemic_offset
           )

write.csv(report_scalars, file.path(output.dir, 'FRED_parameters.csv'), row.names= F, quote = F)
##===============================================##
## Write Slurm job function---------------
##===============================================##
write_submission_array_slurm = function(experiment_supername_in,
                                  experiment_name_in,
                                  experiment_dir_in,
                                  params_base,
                                  job_base,
                                  reps, scalars, FUN, cores_in=1, walltime_in = "0:45:00",
                                  fred_home_dir_in="~/Coronavirus/FRED", fred_results_in="~/Coronavirus/FRED_RESULTS"){
    print('submit array')    
    jobname = sprintf("%s-%s",experiment_supername_in, experiment_name_in)
    tmp_cmd_file = sprintf('tmp_execute_cmd_%s.txt',jobname)
    FUN(scalars, tmp_cmd_file)    
    n = nrow(scalars)
    submission_template = "#!/bin/bash
#SBATCH --job-name=JOBNAME
#SBATCH --array=1-JOBSQUEUE
#SBATCH --nodes=1
#SBATCH --cpus-per-task=JOBCORES
#SBATCH --time=9:00:00                        # Time limit hrs:min:sec
#SBATCH --output=errors_job_%j.log            # Standard output and error log

#module load R/3.5.0
#cd $SLURM_WORKDIR

export FRED_HOME=FREDHOMESTR
export FRED_RESULTS=FREDRESULTSSTR
export PATH=${FRED_HOME}/bin:/shared/home/Azure-ASIS/anaconda3/bin:$PATH
export LD_LIBRARY_PATH=/shared/home/Azure-ASIS/packages/lib:/shared/home/Azure-ASIS/anaconda3/lib

file='TMPCMDFILE'
cmd=`head -n ${SLURM_ARRAY_TASK_ID} $file | tail -n 1`
cd EXPERIMENTDIR
eval $cmd
"    
    submission_str = submission_template %>%
        str_replace_all(pattern="JOBNAME", replacement = jobname) %>%
        str_replace_all(pattern="EXPERIMENTDIR", replacement = experiment_dir_in) %>%
        str_replace_all(pattern="FREDHOMESTR", replacement = fred_home_dir_in) %>%
        str_replace_all(pattern="FREDRESULTSSTR", replacement = fred_results_in) %>%
        str_replace_all(pattern="JOBWALLTIME", replacement = walltime_in) %>%
        str_replace_all(pattern="JOBSQUEUE", replacement = as.character(n)) %>%
        str_replace_all(pattern="PARAMSBASE", replacement = params_base) %>%
        str_replace_all(pattern="JOBBASE", replacement = job_base) %>%
        str_replace_all(pattern="REPS", replacement = as.character(reps)) %>%
        str_replace_all(pattern="TMPCMDFILE", replacement = tmp_cmd_file) %>%
        str_replace_all(pattern="JOBCORES", replacement = as.character(cores_in))
        
    submission_file = sprintf("run_files/%s-%s.sh",experiment_supername_in,experiment_name_in)
    file.connection = file(submission_file)
    write(submission_str,file.connection)
    close(file.connection)
    return(submission_file)
}

##===============================================##
## Submit job function---------------
##===============================================##
submit_jobs = function(experiment_supername_in,
                       experiment_name_in,
                       experiment_dir_in,
                       params_base,
                       job_base,
                       reps, scalars,
                       FUN,
                       cores_in = 1,
                       delete_files = F,
                       sub_array = TRUE,
                       submit_job = TRUE,
                       subsys = "SLURM", walltime_in = "0:45:00",
                       fred_home_dir_in="~/Coronavirus/FRED",
                       fred_results_in="~/Coronavirus/FRED_RESULTS"){
    if(sub_array == TRUE){
        if(subsys == "UGE"){
            submission_file = write_submission_array(
                experiment_supername_in = experiment_supername_in,
                experiment_name_in = experiment_name_in,
                experiment_dir_in = experiment_dir_in,
                params_base = params_base,
                job_base = job_base,
                reps = reps, scalars = scalars, FUN = FUN, cores_in = cores_in, walltime_in=walltime_in,
                fred_home_dir_in=fred_home_dir_in, fred_results_in=fred_results_in)
        }else if(subsys == "PBS"){
            submission_file = write_submission_array_pbs(
                experiment_supername_in = experiment_supername_in,
                experiment_name_in = experiment_name_in,
                experiment_dir_in = experiment_dir_in,
                params_base = params_base,
                job_base = job_base,
                reps = reps, scalars = scalars, FUN = FUN,
                cores_in = cores_in,walltime_in=walltime_in,
                fred_home_dir_in=fred_home_dir_in, fred_results_in=fred_results_in)
        }else if(subsys == "SLURM"){
            submission_file = write_submission_array_slurm(
                experiment_supername_in = experiment_supername_in,
                experiment_name_in = experiment_name_in,
                experiment_dir_in = experiment_dir_in,
                params_base = params_base,
                job_base = job_base,
                reps = reps, scalars = scalars, 
                FUN = FUN, 
                cores_in = cores_in, 
                walltime_in=walltime_in,
                fred_home_dir_in=fred_home_dir_in, 
                fred_results_in=fred_results_in)
        }
        if(submit_job){
            if(subsys == "UGE" | subsys == "UGE"){
                system(sprintf("qsub %s", submission_file))
            }else if(subsys == "SLURM"){
                system(sprintf("/usr/bin/sbatch %s", submission_file))
            }
        }
        if(delete_files == TRUE){
            unlink(submission_file)
        }   
    }
    
}

##===============================================##
## submit to CRC---------------
##===============================================##
## unlink(fred_results_dir,recursive = TRUE)
if(!dir.exists(fred_results_dir)){
    dir.create(fred_results_dir)
}

submit_jobs(experiment_supername_in = sprintf('FRED_LINEAGE_PROJ_asymp_%.2f_FM_%.2f_KSUS_%.2f_V%.0f_VAX%d', asymp_infectivity_in, face_mask_transmission_efficacy_in, kids_susceptibility_age_in, variants_in, vaccination_in),
            experiment_name_in = as.character(state_code),
            experiment_dir_in = output.dir,
            params_base = basename_params,
            job_base = basename_jobs,
            reps = reps_per_job,
            scalars = report_scalars,
            FUN = write_cmd_function,
            cores_in=2,
            walltime_in = "6:00:00",
            subsys="SLURM",
            submit_job = subm_jobs,
            fred_home_dir_in=fred_home, fred_results_in=fred_results_dir)
