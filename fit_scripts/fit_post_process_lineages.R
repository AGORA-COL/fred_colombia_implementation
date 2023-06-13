##===============================#
## Plot model fit for Bogotá
## Author: Guido España
## 2020
##===============================#
## Setup-------------
##===============================#
setwd('/shared/home/Azure-ASIS/FRED_Implementation/fit_scripts')

set.seed(123456)
library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(lubridate)

library(sf)
library(raster)
library(rgdal)
library(deSolve)

args = (commandArgs(TRUE))
projection_label = args[1]
##===============================#
## Functions------------
##===============================#
hosp_pre <- function(t, state, params){
    gam_h = params[['gamma']]
    h = state["h"]
    hh = state["hh"]
    x = state["x"]
    ch = ch_time[t+1]
    max_h = uci_available_time[t + 1] * 0.97
    cancel_t = 0
    if(h + ch - gam_h[t+1]*h > max_h){
        cancel_t = h + ch - max_h - gam_h[t+1]*h
    }
    ##print(sprintf("t %.0f max h: %d, H: %.2f, cancel: %.2f",t, max_h, h, cancel_t))
    dx = cancel_t
    dh = - gam_h[t+1]*h  + ch - cancel_t
    dhh = - gam_h[t+1]*hh  + ch
    return(list(c(dh, dhh, dx)))
}

error.bar <- function(x, y, upper, lower, length=0.01,...){
  arrows(x,upper, x, lower, angle=90, code=3, length=length, ...)
}

get_SSE <- function(data_in, model_in){
    return(sum((data_in - model_in)^2))
}

lowCI <- function(x){
    return(quantile(x, probs = 0.025, na.rm = T))
}
highCI <- function(x){
    return(quantile(x, probs = 0.975, na.rm = T))
}
medCI  <- function(x){
    return(quantile(x, probs = 0.5, na.rm = T))
}


get_loglikelihood_tests <- function(data_pos_in, data_neg_in, model_symp_in, model_pop_in, model_prev_inf_in,
                                    ratio_in, test_sensitivity = 0.777, test_specificity = 1.0){
    ## For some reason some positive cases are negative:
    ## Exclude those data points
    data_pos_in[which(data_pos_in < 0)] = NA
    data_neg_in[which(data_neg_in < 0)] = NA
    
    prevalence = model_symp_in / model_pop_in
    infected_prev = (model_prev_inf_in) / model_pop_in
    
    probC = prevalence / (prevalence + ratio_in * (1 - prevalence))
    probI = infected_prev*ratio_in / (prevalence + ratio_in * (1 - prevalence))
    probU = (1 - infected_prev - prevalence)*ratio_in / (prevalence + ratio_in * (1 - prevalence))

    prob = (probC + probI)*test_sensitivity + probU*(1-test_specificity)
    ll_out = sum(dbinom(data_pos_in, data_pos_in + data_neg_in, prob, log=T), na.rm = T)
    return(-ll_out)
}


get_pos_test_model <- function(model_symp_in, model_pop_in, model_prev_inf_in,
                               ratio_in, test_sensitivity = 0.777,
                               test_specificity = 1.0){
    ratio_in = ratio_in[1]
    prevalence = model_symp_in / model_pop_in
    infected_prev = (model_prev_inf_in) / model_pop_in
    
    probC = prevalence / (prevalence + ratio_in * (1 - prevalence))
    probI = infected_prev*ratio_in / (prevalence + ratio_in * (1 - prevalence))
    probU = (1 - infected_prev - prevalence)*ratio_in / (prevalence + ratio_in * (1 - prevalence))
    
    prob = (probC + probI)*test_sensitivity + probU*(1-test_specificity)
    return(prob)
}


fit_icu_residence <- function(model_in, data_in){                               
    data_in$Change = data_in$UCI_prevalence_bog - lag(data_in$UCI_prevalence_bog)
    data_in = data_in %>% left_join(model_in, by = 'Date')
    
    Rolling_data <- data_in %>%
        tq_mutate(   
            select     = Chosp_median,
            mutate_fun = rollapply, 
            width      = 7,
            align      = "right",
            FUN        = mean,
            na.rm      = TRUE,
            col_rename = "Admis_7"
        ) %>%
        tq_mutate(
            select     = UCI_prevalence_bog,
            mutate_fun = rollapply,
            width      = 7,
            align      = "right",
            FUN        = mean,
            na.rm      = TRUE,
            col_rename = "Ocupadas_7"
        ) %>%
        mutate(RTi = Ocupadas_7 / Admis_7)
    
    min_RTi = min(Rolling_data$RTi, na.rm = T)
    ResidenceTime <- Rolling_data %>%
        dplyr::select(Date, RTi) %>%
        replace_na(list(RTi = min_RTi))
    return(ResidenceTime)
}


##===============================#
## Data ------------
##===============================#
ss = 11001
outdir = '../FRED_run/output/SHORT_FORECAST'
outdir_fit = '../FRED_run/output/CALIBRATION'
asymp_inf = 1.0
fm_ef = 0.73
ksus = 10.0
var_in = 1
vax_in = 70
outdir_st = file.path(outdir, sprintf('FRED_%d_projections_asymp_%.2f_fm_%.2f_ksus_%.2f_var_%.0f_vax_%03d_mov_%s_out', ss, asymp_inf,fm_ef, ksus, var_in,vax_in,projection_label))
outdir_fit_st = file.path(outdir_fit, sprintf('FRED_%d_calibration_asymp_%.2f_fm_%.2f_ksus_%.2f_var_%.0f_vax_%03d_mov_out', ss, asymp_inf,fm_ef, ksus, var_in,vax_in))
interventions_df = read_csv('../input_files/interventions_Colombia.csv') 

print(outdir_st)

brk_ages = c(0, seq(from=5,by = 5, to = 80), 120)
brk_lbls = sprintf("ACF%d_%d", brk_ages[-length(brk_ages)], brk_ages[-1])
brk_data_lbls = sprintf("A%d_%d", brk_ages[-length(brk_ages)], brk_ages[-1])
base_age_df = data.frame(AgeGroup = brk_lbls, Deaths = 0, stringsAsFactors = F)
base_age_data = data.frame(AgeGroup = brk_data_lbls, stringsAsFactors = F)

localidad_list = read_csv('../data/Bogota_localidades_ID.csv') %>%
    filter(Localidad_ID != 20)


BOG_data = read_csv('../FRED_run/data/BOG_covid_death_data.csv') %>%
    mutate(MunCode = as.numeric(MunCode)) %>%
    filter(MunCode %in% interventions_df$State) %>%
    left_join(interventions_df, by = c("MunCode" = "State")) 

validation_data = read_csv(file.path('../FRED_run/data', 'COL_covid_death_data.csv')) %>%
    mutate(MunCode = as.numeric(MunCode)) %>%
    filter(MunCode %in% interventions_df$State) %>%
    left_join(interventions_df, by = c("MunCode" = "State")) 

fit_data = read_csv(file.path(outdir_fit_st, 'COL_covid_death_data.csv')) %>%
    mutate(MunCode = as.numeric(MunCode)) %>%
    filter(MunCode %in% interventions_df$State) %>%
    left_join(interventions_df, by = c("MunCode" = "State")) 

fit_date = as.Date("2021-10-12")
forecast_date = as.Date("2021-12-31")
brk_ages = c(0, seq(from=5,by = 5, to = 80), 120)
brk_lbls = sprintf("ACF%d_%d", brk_ages[-length(brk_ages)], brk_ages[-1])
brk_data_lbls = sprintf("A%d_%d", brk_ages[-length(brk_ages)], brk_ages[-1])
base_age_df = data.frame(AgeGroup = brk_lbls, Deaths = 0, stringsAsFactors = F)
base_age_data = data.frame(AgeGroup = brk_data_lbls, stringsAsFactors = F)

age_data = read_csv(file.path('../FRED_run/data','Age_BOG_covid_data.csv')) %>%
    filter(MunCode == 11001) %>%
    filter(Date <= fit_date) 

age_sero_data = read_csv('../epidata/Bogota_seroprevalence_study_age.csv')
uci_date_str = format(lubridate::today() - 2, "%d%m%Y")
uci_url = sprintf("https://datosabiertos.bogota.gov.co/dataset/adebac34-7974-4aa8-a31d-97b2b6796e4b/resource/ecbfe485-f4a4-40d0-bc2c-628e9e5d3786/download/osb_ocupacion_ucis_covid-19_%s.csv", uci_date_str)

# download.file(uci_url,
#               '../epidata/Bogota_UCI_ocupacion.csv')

uci_df = read.delim('../epidata/Bogota_UCI_ocupacion.csv', fileEncoding = 'latin1', sep = ";")

uci_df = read_delim('../epidata/Bogota_UCI_ocupacion.csv', locale = locale(encoding = "windows-1252"),
                    delim = ";", col_types = lapply(1:ncol(uci_df),function(x)col_character())) %>%
    rename(UCI_prevalence = "Camas UCI ocupadas Covid-19",
           UCI_available = "Total camas UCI COVID 19 reportadas por IPS") %>%
    mutate(UCI_prevalence = as.numeric(str_replace(UCI_prevalence, '\\.','')),
           UCI_available = as.numeric(str_replace(UCI_available, '\\.',''))) %>%
    mutate(Date = as.Date(parse_date_time(Fecha, orders = c("dmy", "dmy HMS")))) %>%
    filter(!is.na(Date))

imports_uci = data.frame(Prop = c(1,2,2,3,5,6,12,10,10,11,13,12)/100, stringsAsFactors = F)
imports_uci$Date = seq(from = as.Date('2021-04-01'), by = 'day', length.out = nrow(imports_uci))

uci_incidence_file = '../input_files/BOG_UCI_timeseries.csv'

uci_incidence = read_csv(uci_incidence_file)
uci_incidence = uci_incidence[1:(nrow(uci_incidence) - 1),]

jpeg('../figures/report_figure_model_fit_UCI_data_demand.jpeg', width=6.5,height=8, units="in", res = 300)
plot(uci_incidence$Date, uci_incidence$ICU_admissions, type = 'o', col = 'red', lwd = 1, ylab = 'Ingresos UCI',xlab = 'Fecha', xlim = c(as.Date('2021-01-01'), max(uci_incidence$Date)), ylim = c(0,320))
lines(uci_incidence$Date, uci_incidence$ICU_admissions_demand, col = 'orange')
dev.off()

uci_df = uci_df %>% left_join(imports_uci, by = 'Date') %>%
    replace_na(list(Prop = 0)) %>%
    mutate(UCI_prevalence_bog = UCI_prevalence * (1 - Prop))

## Read fitted hosp parameters
icu_duration = read_csv('../output/FRED_fit_ICU_residence.csv')

##===============================#
## Process output-------------
##===============================#
data_out = file.path(outdir_st,'fred_output.csv')

params_out = file.path(outdir_st, 'FRED_parameters_out.csv')
params_sweep_df = read_csv(params_out)
fred_sweep_df = read_csv(data_out) %>% 
    right_join(params_sweep_df, by = c("job_id" = "job_id")) %>%
    mutate(Date = as.Date("2020-01-01") + Day)    
fred_sweep_df[is.na(fred_sweep_df)] <- 0 

fred_sweep_df = fred_sweep_df %>% group_by(seed,state_code) %>% mutate(CumCF = cumsum(CF_mean), RR_max = max(RR_mean, na.rm = T)) %>% ungroup()

param_names_fit = c('workplace_contact_factor', 'neighborhood_contact_factor', 'community_contact_rate_1', 'shelter_in_place_compliance', 'influenza_transmissibility',
                    'facemask_compliance',
                    'variantalpha_transmissibility_factor',
                    'variantalpha_cross_protection_prob',
                    'variantalpha_imports_factor',
                    'variantgamma_imports_factor', 
                    'variantkappa_introduction_day',
                    'variantdelta_imports_factor',
                    'variantgamma_transmissibility_factor',
                    'variantgamma_cross_protection_prob',
                    'variantkappa_transmissibility_factor', 
                    'variantkappa_cross_protection_prob',
                    'variantdelta_transmissibility_factor', 
                    'variantdelta_cross_protection_prob',
                    'variantomicron_transmissibility_factor', 
                    'variantomicron_cross_protection_prob',
                    'neighborhood_same_age_bias')

(summary(params_sweep_df[,param_names_fit]))

##===============================#
## process  Model fit --------
##===============================#
## First, summarize all data
                                                                                         
state_fred = fred_sweep_df %>%
    group_by(job_id) %>%
    mutate(AR_mean = 100 * (cumsum(C_mean + C_1_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean + C_6_mean - reinfections)) / N_mean[1],
           AR_0_mean = 100 * (cumsum(C_mean)) / N_mean[1],
           AR_1_mean = 100 * (cumsum(C_1_mean)) / N_mean[1],
           AR_2_mean = 100 * (cumsum(C_2_mean)) / N_mean[1],
           AR_3_mean = 100 * (cumsum(C_3_mean)) / N_mean[1],
           AR_4_mean = 100 * (cumsum(C_4_mean)) / N_mean[1],
           AR_5_mean = 100 * (cumsum(C_5_mean)) / N_mean[1],
           AR_6_mean = 100 * (cumsum(C_6_mean)) / N_mean[1],
           VAX_1 = V1A_0_4_mean + V1A_5_9_mean + V1A_10_14_mean + V1A_15_19_mean + V1A_20_24_mean + V1A_25_29_mean + V1A_30_34_mean + V1A_35_39_mean + V1A_40_44_mean + V1A_45_49_mean + V1A_50_54_mean + V1A_55_59_mean  + V1A_60_64_mean + V1A_65_69_mean  + V1A_70_74_mean + V1A_75_79_mean  + V1A_80_120_mean,
           VAX_2 = V2A_0_4_mean + V2A_5_9_mean + V2A_10_14_mean + V2A_15_19_mean + V2A_20_24_mean + V2A_25_29_mean + V2A_30_34_mean + V2A_35_39_mean + V2A_40_44_mean + V2A_45_49_mean + V2A_50_54_mean + V2A_55_59_mean  + V2A_60_64_mean + V2A_65_69_mean  + V2A_70_74_mean + V2A_75_79_mean  + V2A_80_120_mean,
           VAX_3 = V3A_0_4_mean + V3A_5_9_mean + V3A_10_14_mean + V3A_15_19_mean + V3A_20_24_mean + V3A_25_29_mean + V3A_30_34_mean + V3A_35_39_mean + V3A_40_44_mean + V3A_45_49_mean + V3A_50_54_mean + V3A_55_59_mean  + V3A_60_64_mean + V3A_65_69_mean  + V3A_70_74_mean + V3A_75_79_mean  + V3A_80_120_mean
          ) %>%
    ungroup()


tmp_fred_all = state_fred %>%
    group_by(intervention_id, Day, Date, start_date, school_closure_day) %>%
    summarize(VAX_1_median = quantile(VAX_1, probs = c(0.5), na.rm = T),
              VAX_1_low = quantile(VAX_1, probs = c(0.025), na.rm = T),
              VAX_1_high = quantile(VAX_1, probs = c(0.975), na.rm = T),
              
              VAX_2_median = quantile(VAX_2, probs = c(0.5), na.rm = T),
              VAX_2_low = quantile(VAX_2, probs = c(0.025), na.rm = T),
              VAX_2_high = quantile(VAX_2, probs = c(0.975), na.rm = T),
              
              VAX_3_median = quantile(VAX_3, probs = c(0.5), na.rm = T),
              VAX_3_low = quantile(VAX_3, probs = c(0.025), na.rm = T),
              VAX_3_high = quantile(VAX_3, probs = c(0.975), na.rm = T),
        
              CF_median = quantile(CF_mean + CF_1_mean + CF_2_mean + CF_3_mean + CF_4_mean + CF_5_mean + CF_6_mean, probs = c(0.5), na.rm = T),
              CF_low = quantile(CF_mean + CF_1_mean + CF_2_mean + CF_3_mean + CF_4_mean + CF_5_mean + CF_6_mean, probs = c(0.025), na.rm = T),
              CF_high = quantile(CF_mean + CF_1_mean + CF_2_mean + CF_3_mean +  CF_4_mean + CF_5_mean + CF_6_mean, probs = c(0.975), na.rm = T),
              Cs_median = quantile(Cs_mean + Cs_1_mean + Cs_2_mean + Cs_3_mean + Cs_4_mean + Cs_5_mean + Cs_6_mean, probs = c(0.5), na.rm = T),
              Cs_low = quantile(Cs_mean + Cs_1_mean + Cs_2_mean + Cs_3_mean + Cs_4_mean + Cs_5_mean + Cs_6_mean, probs = c(0.025), na.rm = T),
              Cs_high = quantile(Cs_mean + Cs_1_mean + Cs_2_mean + Cs_3_mean + Cs_4_mean + Cs_5_mean + Cs_6_mean, probs = c(0.975), na.rm = T),
              reinfections_median = quantile(reinfections / (C_1_mean + C_mean+ C_2_mean + C_3_mean + C_4_mean + C_5_mean + C_6_mean), probs = c(0.5), na.rm = T),
              reinfections_low = quantile(reinfections / (C_1_mean + C_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean +C_6_mean), probs = c(0.025), na.rm = T),
              reinfections_high = quantile(reinfections / (C_1_mean + C_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean + C_6_mean), probs = c(0.975), na.rm = T),
              C_median = quantile(C_mean + C_1_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean + C_6_mean, probs = c(0.5), na.rm = T),
              C_low = quantile(C_mean + C_1_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean + C_6_mean, probs = c(0.025), na.rm = T),
              C_high = quantile(C_mean + C_1_mean+ C_2_mean + C_3_mean + C_4_mean + C_5_mean + C_6_mean, probs = c(0.975), na.rm = T),
              
              C_0_median = quantile(C_mean, probs = c(0.5), na.rm = T),
              C_0_low = quantile(C_mean, probs = c(0.025), na.rm = T),
              C_0_high = quantile(C_mean, probs = c(0.975), na.rm = T),
              
              C_1_median = quantile(C_1_mean, probs = c(0.5), na.rm = T),
              C_1_low = quantile(C_1_mean, probs = c(0.025), na.rm = T),
              C_1_high = quantile(C_1_mean, probs = c(0.975), na.rm = T),
              
              C_2_median = quantile(C_2_mean, probs = c(0.5), na.rm = T),
              C_2_low = quantile(C_2_mean, probs = c(0.025), na.rm = T),
              C_2_high = quantile(C_2_mean, probs = c(0.975), na.rm = T),
              
              C_3_median = quantile(C_3_mean, probs = c(0.5), na.rm = T),
              C_3_low = quantile(C_3_mean, probs = c(0.025), na.rm = T),
              C_3_high = quantile(C_3_mean, probs = c(0.975), na.rm = T),
              
              C_4_median = quantile(C_4_mean, probs = c(0.5), na.rm = T),
              C_4_low = quantile(C_4_mean, probs = c(0.025), na.rm = T),
              C_4_high = quantile(C_4_mean, probs = c(0.975), na.rm = T),
              
              C_5_median = quantile(C_5_mean, probs = c(0.5), na.rm = T),
              C_5_low = quantile(C_5_mean, probs = c(0.025), na.rm = T),
              C_5_high = quantile(C_5_mean, probs = c(0.975), na.rm = T),
              
              C_6_median = quantile(C_6_mean, probs = c(0.5), na.rm = T),
              C_6_low = quantile(C_6_mean, probs = c(0.025), na.rm = T),
              C_6_high = quantile(C_6_mean, probs = c(0.975), na.rm = T),
              
              CF_0_median = quantile(CF_mean, probs = c(0.5), na.rm = T),
              CF_0_low = quantile(CF_mean, probs = c(0.025), na.rm = T),
              CF_0_high = quantile(CF_mean, probs = c(0.975), na.rm = T),
              
              CF_1_median = quantile(CF_1_mean, probs = c(0.5), na.rm = T),
              CF_1_low = quantile(CF_1_mean, probs = c(0.025), na.rm = T),
              CF_1_high = quantile(CF_1_mean, probs = c(0.975), na.rm = T),
              
              CF_2_median = quantile(CF_2_mean, probs = c(0.5), na.rm = T),
              CF_2_low = quantile(CF_2_mean, probs = c(0.025), na.rm = T),
              CF_2_high = quantile(CF_2_mean, probs = c(0.975), na.rm = T),
              
              CF_3_median = quantile(CF_3_mean, probs = c(0.5), na.rm = T),
              CF_3_low = quantile(CF_3_mean, probs = c(0.025), na.rm = T),
              CF_3_high = quantile(CF_3_mean, probs = c(0.975), na.rm = T),
              
              CF_4_median = quantile(CF_4_mean, probs = c(0.5), na.rm = T),
              CF_4_low = quantile(CF_4_mean, probs = c(0.025), na.rm = T),
              CF_4_high = quantile(CF_4_mean, probs = c(0.975), na.rm = T),
              
              CF_5_median = quantile(CF_5_mean, probs = c(0.5), na.rm = T),
              CF_5_low = quantile(CF_5_mean, probs = c(0.025), na.rm = T),
              CF_5_high = quantile(CF_5_mean, probs = c(0.975), na.rm = T),
              
              CF_6_median = quantile(CF_6_mean, probs = c(0.5), na.rm = T),
              CF_6_low = quantile(CF_6_mean, probs = c(0.025), na.rm = T),
              CF_6_high = quantile(CF_6_mean, probs = c(0.975), na.rm = T),
              
              V_median = quantile(V_mean, probs = c(0.5), na.rm = T),
              V_low = quantile(V_mean, probs = c(0.025), na.rm = T),
              V_high = quantile(V_mean, probs = c(0.975), na.rm = T),
              Vtd_median = quantile(Vtd_mean, probs = c(0.5), na.rm = T),
              Vtd_low = quantile(Vtd_mean, probs = c(0.025), na.rm = T),
              Vtd_high = quantile(Vtd_mean, probs = c(0.975), na.rm = T),
              
              DomOriginal_median = quantile(C_mean / (C_1_mean + C_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean + C_6_mean), probs = c(0.5), na.rm = T),
              DomOriginal_low = quantile(C_mean / (C_1_mean + C_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean + C_6_mean), probs = c(0.025), na.rm = T),
              DomOriginal_high = quantile(C_mean / (C_1_mean + C_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean + C_6_mean), probs = c(0.975), na.rm = T),
              
              DomAlpha_median = quantile(C_1_mean / (C_1_mean + C_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean + C_6_mean), probs = c(0.5), na.rm = T),
              DomAlpha_low = quantile(C_1_mean / (C_1_mean + C_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean + C_6_mean), probs = c(0.025), na.rm = T),
              DomAlpha_high = quantile(C_1_mean /(C_1_mean + C_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean + C_6_mean), probs = c(0.975), na.rm = T),
              
              DomGamma_median = quantile(C_2_mean / (C_1_mean + C_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean + C_6_mean), probs = c(0.5), na.rm = T),
              DomGamma_low = quantile(C_2_mean / (C_1_mean + C_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean + C_6_mean), probs = c(0.025), na.rm = T),
              DomGamma_high = quantile(C_2_mean /(C_1_mean + C_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean + C_6_mean), probs = c(0.975), na.rm = T),
              
              DomKappa_median = quantile(C_3_mean / (C_1_mean + C_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean + C_6_mean), probs = c(0.5), na.rm = T),
              DomKappa_low = quantile(C_3_mean / (C_1_mean + C_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean + C_6_mean), probs = c(0.025), na.rm = T),
              DomKappa_high = quantile(C_3_mean /(C_1_mean + C_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean + C_6_mean), probs = c(0.975), na.rm = T),
              
              DomDelta_median = quantile(C_4_mean / (C_1_mean + C_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean + C_6_mean), probs = c(0.5), na.rm = T),
              DomDelta_low = quantile(C_4_mean / (C_1_mean + C_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean + C_6_mean), probs = c(0.025), na.rm = T),
              DomDelta_high = quantile(C_4_mean /(C_1_mean + C_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean + C_6_mean), probs = c(0.975), na.rm = T),
              
              DomOmicron_median = quantile(C_5_mean / (C_1_mean + C_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean + C_6_mean), probs = c(0.5), na.rm = T),
              DomOmicron_low = quantile(C_5_mean / (C_1_mean + C_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean + C_6_mean), probs = c(0.025), na.rm = T),
              DomOmicron_high = quantile(C_5_mean /(C_1_mean + C_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean + C_6_mean), probs = c(0.975), na.rm = T),
              
              DomOmicronBAX_median = quantile(C_6_mean / (C_1_mean + C_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean + C_6_mean), probs = c(0.5), na.rm = T),
              DomOmicronBAX_low = quantile(C_6_mean / (C_1_mean + C_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean + C_6_mean), probs = c(0.025), na.rm = T),
              DomOmicronBAX_high = quantile(C_6_mean /(C_1_mean + C_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean + C_6_mean), probs = c(0.975), na.rm = T),
              
              Phosp_median = 0,
              Phosp_low = 0,
              Phosp_high = 0,
              Chosp_median = quantile(Chosp_mean + Chosp_1_mean + Chosp_2_mean +  Chosp_3_mean  + Chosp_4_mean + Chosp_5_mean + Chosp_6_mean, probs = c(0.5), na.rm = T),
              Chosp_low = quantile(Chosp_mean + Chosp_1_mean + Chosp_2_mean + Chosp_3_mean + Chosp_4_mean + Chosp_5_mean + Chosp_6_mean, probs = c(0.025), na.rm = T),
              Chosp_high = quantile(Chosp_mean + Chosp_1_mean + Chosp_2_mean + Chosp_3_mean + Chosp_4_mean + Chosp_5_mean + Chosp_6_mean, probs = c(0.975), na.rm = T),
              PrevInf_median = quantile(PrevInf_mean, probs = c(0.5), na.rm = T),
              PrevInf_low = quantile(PrevInf_mean, probs = c(0.025), na.rm = T),
              PrevInf_high = quantile(PrevInf_mean, probs = c(0.975), na.rm = T),
              AR_median = medCI(AR_mean)/100,
              AR_low = lowCI(AR_mean)/100,
              AR_high = highCI(AR_mean)/100,
              Sch_median = mean(Sch_mean),
              Wrk_median = mean(Wrk_mean),
              Nbr_median = mean(Nbr_mean),
              H_median = medCI(H_mean),
              N_median = medCI(N_mean),
              Nursing_home_CF_median = medCI(Nursing_Home_CF_mean),
              N_sheltering_median = 1 -  medCI(N_sheltering_mean/N_mean[1]),
              N_sheltering_low = 1-  lowCI(N_sheltering_mean/N_mean[1]),
              N_sheltering_high = 1 -  highCI(N_sheltering_mean/N_mean[1]),
              
              AR_0_low = quantile(AR_0_mean, probs = c(0.025), na.rm = T),
              AR_0_high = quantile(AR_0_mean, probs = c(0.975), na.rm = T),
              AR_0_median = quantile(AR_0_mean, probs = c(0.5), na.rm = T),
              
              AR_1_low = quantile(AR_1_mean, probs = c(0.025), na.rm = T),
              AR_1_high = quantile(AR_1_mean, probs = c(0.975), na.rm = T),
              AR_1_median = quantile(AR_1_mean, probs = c(0.5), na.rm = T),
              
              AR_2_low = quantile(AR_2_mean, probs = c(0.025), na.rm = T),
              AR_2_high = quantile(AR_2_mean, probs = c(0.975), na.rm = T),
              AR_2_median = quantile(AR_2_mean, probs = c(0.5), na.rm = T),
              
              AR_3_low = quantile(AR_3_mean, probs = c(0.025), na.rm = T),
              AR_3_high = quantile(AR_3_mean, probs = c(0.975), na.rm = T),
              AR_3_median = quantile(AR_3_mean, probs = c(0.5), na.rm = T),
              
              AR_4_low = quantile(AR_4_mean, probs = c(0.025), na.rm = T),
              AR_4_high = quantile(AR_4_mean, probs = c(0.975), na.rm = T),
              AR_4_median = quantile(AR_4_mean, probs = c(0.5), na.rm = T),
              
              AR_5_low = quantile(AR_5_mean, probs = c(0.025), na.rm = T),
              AR_5_high = quantile(AR_5_mean, probs = c(0.975), na.rm = T),
              AR_5_median = quantile(AR_5_mean, probs = c(0.5), na.rm = T),
              
              AR_6_low = quantile(AR_6_mean, probs = c(0.025), na.rm = T),
              AR_6_high = quantile(AR_6_mean, probs = c(0.975), na.rm = T),
              AR_6_median = quantile(AR_6_mean, probs = c(0.5), na.rm = T),
              
              RR_median = medCI(RR_mean) * (1 - (DomAlpha_median + DomGamma_median + DomKappa_median + DomDelta_median + DomOmicron_median + DomOmicronBAX_median)),
              RR_1_median = medCI(RR_1_mean) * DomAlpha_median,
              RR_2_median = medCI(RR_2_mean) * DomGamma_median,
              RR_3_median = medCI(RR_3_mean) * DomKappa_median,
              RR_4_median = medCI(RR_4_mean) * DomDelta_median,   
              RR_5_median = medCI(RR_5_mean) * DomOmicron_median, 
              RR_6_median = medCI(RR_6_mean) * DomOmicronBAX_median, 
              shelter_in_place_delay_mean = mean(shelter_in_place_delay_mean))%>%                               
    ungroup() %>%
    mutate(RR_total = RR_median + RR_1_median + RR_2_median + RR_3_median + RR_4_median + RR_5_median + RR_6_median)

## adjust uci occupancy -----------

fred_adjusted = tibble()
for(in_id in unique(tmp_fred_all$intervention_id)){
    tmp_df = tmp_fred_all %>% filter(intervention_id == in_id) %>%
        left_join(icu_duration, by = 'Date') %>%
        left_join(uci_df[,c('Date', 'UCI_available')], by = 'Date') %>%
        replace_na(list(ResidenceTimeICU = 11))
    tmp_df$UCI_available[is.na(tmp_df$UCI_available)] = max(tmp_df$UCI_available, na.rm = T)
    ## replace_na(list(ResidenceTimeICU = icu_duration$ResidenceTimeICU[nrow(icu_duration)]))
    
    for(mm in c('median','high','low')){
        var_cancelled = sprintf('UCI_Cancelled_%s',mm)
        var_inci = sprintf('Chosp_%s',mm)
        var_pre = sprintf('Phosp_%s',mm)
        var_pre_dem = sprintf('PhospDem_%s',mm)
        ch_time = tmp_df %>% pull(var_inci)
        uci_available_time = tmp_df %>% pull(UCI_available)
        ode.res = ode(y = c(h = 0, hh = 0, x = 0), times = seq(from = 0, to = nrow(tmp_df) - 1, by = 1), func = hosp_pre, parms = list(gamma = 1/tmp_df$ResidenceTimeICU))
        tmp_df[[var_pre]] = ode.res[,'h']
        tmp_df[[var_pre_dem]] = ode.res[,'hh']
        tmp_df[[var_cancelled]] = c(0,diff(ode.res[,'x']))
    }
    fred_adjusted = bind_rows(fred_adjusted, tmp_df)
}
tmp_fred_all = fred_adjusted

tmp_fred = tmp_fred_all %>% filter(intervention_id == "default")

tmp_fred_vax = tmp_fred_all %>% filter(intervention_id == "default_vax")

tmp_fred_age = fred_sweep_df %>% drop_na() %>%
    filter(intervention_id == "high_contacts") %>%
    dplyr::select(job_id, Day, Date,school_closure_day, start_date,  matches('^ARA[0-9]+'))%>%
    gather(key = AgeGroup, value = Infections, -c('job_id', 'Day', 'Date', 'school_closure_day', 'start_date')) %>%
    mutate(AgeGroup = str_replace(AgeGroup, "AR","")) %>%
    group_by(Day, Date, start_date, school_closure_day, AgeGroup) %>%
    summarize(C_median = median(Infections))%>%
    ungroup()

variants_list = c('Alpha', 'Gamma', 'Kappa', 'Delta')
variants_names_list = c('Alpha', 'Gamma', 'B.1.621', 'Delta')                                                    
write_csv(tmp_fred_all, sprintf('../output/fred_output_model_fit_%s.csv', projection_label))