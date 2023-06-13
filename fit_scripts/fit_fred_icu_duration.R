library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(lubridate)

library(deSolve)
library(splines)
library(fda)

##===============================#
## Functions------------
##===============================#
lowCI <- function(x){
    return(quantile(x, probs = 0.025, na.rm = T))
}
highCI <- function(x){
    return(quantile(x, probs = 0.975, na.rm = T))
}
medCI  <- function(x){
    return(quantile(x, probs = 0.5, na.rm = T))
}


hosp_pre <- function(t, state, params){
    gam_h = params[['gamma']]
    h = state["h"]
    ch = ch_time[t+1]
    dh = - gam_h[t+1]*h + ch
    return(list(dh))
}

llik <- function(params){
    params_in = params
    fs.bs = fd(coef = params_in, basisobj = bs)
    gamma_t = 1 / predict(fs.bs,1:nrow(data_in))
    params.ode = list(gamma = gamma_t)
    ode.res = ode(y = init_state, times = tt, func = hosp_pre, parms = params.ode, method = "ode23")
    fit_d = ode.res[,2]
    lik_out = -sum(dpois(data_in$UCI_prevalence, lambda = fit_d + 0.01, log = T))
    return(lik_out)
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
vax_in = 100
outdir_st = file.path(outdir, sprintf('FRED_%d_postcalibration_asymp_%.2f_fm_%.2f_ksus_%.2f_var_%.0f_vax_%03d_mov_out', ss, asymp_inf,fm_ef, ksus, var_in, vax_in))
outdir_fit_st = file.path(outdir_fit, sprintf('FRED_%d_calibration_asymp_%.2f_fm_%.2f_ksus_%.2f_var_%.0f_vax_%03d_mov_out', ss, asymp_inf,fm_ef, ksus, var_in, vax_in))
interventions_df = read_csv('../FRED_run/input_files/interventions_Colombia.csv') 


brk_ages = c(0, seq(from=10,by = 10, to = 80), 120)
brk_lbls = sprintf("ACF%d_%d", brk_ages[-length(brk_ages)], brk_ages[-1])
brk_data_lbls = sprintf("A%d_%d", brk_ages[-length(brk_ages)], brk_ages[-1])
base_age_df = data.frame(AgeGroup = brk_lbls, Deaths = 0, stringsAsFactors = F)
base_age_data = data.frame(AgeGroup = brk_data_lbls, stringsAsFactors = F)

localidad_list = read_csv('../input_files/Bogota_localidades_ID.csv') %>%
    filter(Localidad_ID != 20)


BOG_data = read_csv('../FRED_run/BOG_covid_death_data.csv') %>%
    mutate(MunCode = as.numeric(MunCode)) %>%
    filter(MunCode %in% interventions_df$State) %>%
    left_join(interventions_df, by = c("MunCode" = "State")) 

validation_data = read_csv(file.path('../FRED_run/input_files', 'COL_covid_death_data.csv')) %>%
    mutate(MunCode = as.numeric(MunCode)) %>%
    filter(MunCode %in% interventions_df$State) %>%
    left_join(interventions_df, by = c("MunCode" = "State")) 

fit_data = read_csv(file.path(outdir_fit_st, 'COL_covid_death_data.csv')) %>%
    mutate(MunCode = as.numeric(MunCode)) %>%
    filter(MunCode %in% interventions_df$State) %>%
    left_join(interventions_df, by = c("MunCode" = "State")) 

fit_date = as.Date("2021-07-31")
brk_ages = c(0, seq(from=10,by = 10, to = 80), 120)
brk_lbls = sprintf("ACF%d_%d", brk_ages[-length(brk_ages)], brk_ages[-1])
brk_data_lbls = sprintf("A%d_%d", brk_ages[-length(brk_ages)], brk_ages[-1])
base_age_df = data.frame(AgeGroup = brk_lbls, Deaths = 0, stringsAsFactors = F)
base_age_data = data.frame(AgeGroup = brk_data_lbls, stringsAsFactors = F)

age_data = read_csv(file.path('../FRED_run','Age_BOG_covid_data.csv')) %>%
    mutate(MunCode = 11001) %>%
    filter(Date <= fit_date, MunCode %in% interventions_df$State) 

age_sero_data = read_csv('../epidata/Bogota_seroprevalence_study_age.csv')

hosp_data = readxl::read_xlsx('../epidata/2021-03-26_Incidencia_Hosp.xlsx') %>%
    mutate(Date = as.Date(FECHA_DE_HOSPITALIZACION))

rt_url = "https://saludata.saludcapital.gov.co/osb/wp-content/uploads/medios/Numero-de-Reproduccion-Efectivo-R(t)%20todos%20los%20casos.csv"
download.file(rt_url,'../epidata/Bogota_Reproduction_Number.csv')

data_rt = read.delim('../epidata/Bogota_Reproduction_Number.csv', fileEncoding  = "latin1" , skip = 5, sep = ";") %>%
    mutate(DateStart = parse_date_time(Fecha_Inicio_Ventana, orders = c("dmy", "dmy HMS")))%>%
    mutate(Date = parse_date_time(Fecha_Fin_Ventana, orders = c("dmy", "dmy HMS")))%>%
    mutate(R_mean = `Mean.R.`, R_std = `Std.R.`)%>%
    mutate(R_mean = as.numeric(str_replace(R_mean, ",", "\\."))) %>%
    filter(grepl(".*Bogot",Localidad))

uci_date_str = format(lubridate::today() - 1, "%d%m%Y")

uci_url = sprintf("https://datosabiertos.bogota.gov.co/dataset/adebac34-7974-4aa8-a31d-97b2b6796e4b/resource/6657cbbc-8277-4caf-8823-3350f253ae64/download/osb_ocupacion_ucis_covid-19_%s.csv", uci_date_str)

download.file(uci_url,
              '../epidata/Bogota_UCI_ocupacion.csv')

uci_df = read.delim('../epidata/Bogota_UCI_ocupacion.csv', fileEncoding = 'latin1', sep = ";")

uci_df = read_delim('../epidata/Bogota_UCI_ocupacion.csv', locale = locale(encoding = "windows-1252"),
                    delim = ";", col_types = lapply(1:ncol(uci_df),function(x)col_character())) %>%
    rename(UCI_prevalence = "Camas UCI ocupadas Covid-19") %>%
    mutate(UCI_prevalence = as.numeric(str_replace(UCI_prevalence, '\\.',''))) %>%
    mutate(Date = as.Date(parse_date_time(Fecha, orders = c("dmy", "dmy HMS")))) %>%
    filter(!is.na(Date))

uci_incidence_file = '../input_files/BOG_UCI_timeseries.csv'

uci_incidence = read_csv(uci_incidence_file)
uci_incidence = uci_incidence[1:(nrow(uci_incidence) - 1),]

##===============================#
## Process output-------------
##===============================#
data_out = file.path(outdir_st,'fred_output.csv')

params_out = file.path(outdir_st, 'FRED_parameters_out.csv')
params_sweep_df = read_csv(params_out)
fred_sweep_df = read_csv(data_out) %>% 
    right_join(params_sweep_df, by = c("job_id" = "job_id")) %>%
    mutate(Date = as.Date("2020-01-01") + Day)    


fred_sweep_df = fred_sweep_df %>% group_by(seed,state_code) %>% mutate(CumCF = cumsum(CF_mean), RR_max = max(RR_mean, na.rm = T)) %>% ungroup()

param_names_fit = c('workplace_contact_factor', 'neighborhood_contact_factor', 'community_contact_rate_1', 'shelter_in_place_compliance', 'influenza_transmissibility', 'facemask_compliance')

(summary(params_sweep_df[,param_names_fit]))

##===============================#
## process  Model fit --------
##===============================#
## First, summarize all data
state_fred = fred_sweep_df %>% drop_na() %>%
    group_by(job_id) %>%
    mutate(AR_mean = 100 * (cumsum(C_mean) + cumsum(C_1_mean)) / N_mean[1]) %>%
    ungroup()

tmp_fred_all = state_fred %>%
    group_by(intervention_id, Day, Date, start_date, school_closure_day) %>%
    summarize(CF_median = quantile(CF_mean + CF_1_mean , probs = c(0.5), na.rm = T),
              CF_low = quantile(CF_mean + CF_1_mean, probs = c(0.025), na.rm = T),
              CF_high = quantile(CF_mean + CF_1_mean, probs = c(0.975), na.rm = T),
              Cs_median = quantile(Cs_mean + Cs_1_mean, probs = c(0.5), na.rm = T),
              Cs_low = quantile(Cs_mean + Cs_1_mean, probs = c(0.025), na.rm = T),
              Cs_high = quantile(Cs_mean + Cs_1_mean, probs = c(0.975), na.rm = T),
              reinfections_median = quantile(reinfections / (C_1_mean + C_mean), probs = c(0.5), na.rm = T),
              reinfections_low = quantile(reinfections / (C_1_mean + C_mean), probs = c(0.025), na.rm = T),
              reinfections_high = quantile(reinfections / (C_1_mean + C_mean), probs = c(0.975), na.rm = T),
              C_median = quantile(C_mean + C_1_mean, probs = c(0.5), na.rm = T),
              C_low = quantile(C_mean + C_1_mean, probs = c(0.025), na.rm = T),
              C_high = quantile(C_mean + C_1_mean, probs = c(0.975), na.rm = T),
              C_1_median = quantile(C_1_mean, probs = c(0.5), na.rm = T),
              C_1_low = quantile(C_1_mean, probs = c(0.025), na.rm = T),
              C_1_high = quantile(C_1_mean, probs = c(0.975), na.rm = T),
              DomVariant_median = quantile(C_1_mean / (C_1_mean + C_mean), probs = c(0.5), na.rm = T),
              DomVariant_low = quantile(C_1_mean / (C_mean + C_1_mean), probs = c(0.025), na.rm = T),
              DomVariant_high = quantile(C_1_mean / (C_mean + C_1_mean), probs = c(0.975), na.rm = T),
              A5_15_median = quantile(A5_15, probs = c(0.5), na.rm = T) / 1073709,
              A35_65_median = quantile(A15_25 + A25_35 + A35_45 + A45_55 + A55_65, probs = c(0.5), na.rm = T) / 2432587,
              Phosp_median = quantile(Phosp_mean + Phosp_1_mean, probs = c(0.5), na.rm = T),
              Phosp_low = quantile(Phosp_mean + Phosp_1_mean, probs = c(0.025), na.rm = T),
              Phosp_high = quantile(Phosp_mean + Phosp_1_mean, probs = c(0.975), na.rm = T),
              Chosp_median = quantile(Chosp_mean + Chosp_1_mean, probs = c(0.5), na.rm = T),
              Chosp_low = quantile(Chosp_mean + Chosp_1_mean, probs = c(0.025), na.rm = T),
              Chosp_high = quantile(Chosp_mean + Chosp_1_mean, probs = c(0.975), na.rm = T),
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
              RR_median = medCI(RR_mean),
              RR_1_median = medCI(RR_1_mean),              
              shelter_in_place_delay_mean = mean(shelter_in_place_delay_mean))%>%
    ungroup()

tmp_fred = tmp_fred_all %>% filter(intervention_id == "default")
tmp_fred$RR_total = (tmp_fred$RR_1_median * tmp_fred$DomVariant_median + tmp_fred$RR_median * (1 - tmp_fred$DomVariant_median))
# tmp_fred_mov = tmp_fred_all %>% filter(intervention_id == "high_mov") 
# tmp_fred_mov$RR_total = (tmp_fred_mov$RR_1_median * tmp_fred_mov$DomVariant_median + tmp_fred_mov$RR_median * (1 - tmp_fred_mov$DomVariant_median))
# tmp_fred_contacts = tmp_fred_all %>% filter(intervention_id == "high_contacts") 
# tmp_fred_contacts$RR_total = (tmp_fred_contacts$RR_1_median * tmp_fred_contacts$DomVariant_median + tmp_fred_contacts$RR_median * (1 - tmp_fred_contacts$DomVariant_median))


tmp_fred_vax = tmp_fred_all %>% filter(intervention_id == "default_vax")
tmp_fred_vax$RR_total = (tmp_fred_vax$RR_1_median * tmp_fred_vax$DomVariant_median + tmp_fred_vax$RR_median * (1 - tmp_fred_vax$DomVariant_median))
# tmp_fred_mov_vax = tmp_fred_all %>% filter(intervention_id == "high_mov_vax") 
# tmp_fred_mov_vax$RR_total = (tmp_fred_mov_vax$RR_1_median * tmp_fred_mov_vax$DomVariant_median + tmp_fred_mov_vax$RR_median * (1 - tmp_fred_mov_vax$DomVariant_median))
# tmp_fred_contacts_vax = tmp_fred_all %>% filter(intervention_id == "high_contacts_vax") 
# tmp_fred_contacts_vax$RR_total = (tmp_fred_contacts_vax$RR_1_median * tmp_fred_contacts_vax$DomVariant_median + tmp_fred_contacts_vax$RR_median * (1 - tmp_fred_contacts_vax$DomVariant_median))


##===============================#
## Fit occupancy model------------
##===============================#
## Data
uci_incidence$Chosp_median = uci_incidence$ICU_admissions

data_in = left_join(uci_df, tmp_fred_mov, by = 'Date')
##data_in = left_join(uci_df, tmp_fred_mov_vax, by = 'Date')
##data_in = left_join(uci_df, uci_incidence, by = 'Date')

ch_time = data_in$Chosp_median

tt = seq(from=0,to = nrow(data_in)-1, by = 1)

init_state = c(h = 0)
nodes_bs = round(nrow(data_in)/15)
bs = fda::create.bspline.basis(rangeval = c(1,nrow(data_in)), nbasis = nodes_bs, norder = 3)
llik(rep(10, nodes_bs))

init_nodes = rep(11,nodes_bs)


par_fit = optim(par = init_nodes, fn = llik, control = list(maxit = 5e3), method = "L-BFGS-B", lower = init_nodes - 5, upper = init_nodes + 20)

fit_bs = par_fit$par
fs.bs = fd(coef = fit_bs, basisobj = bs)
gamma_t = 1 / predict(fs.bs,1:nrow(data_in))
params.ode = list(gamma = gamma_t)
ode.res = ode(y = init_state, times = tt, func = hosp_pre, parms = params.ode, method = "ode23")
fit_d = ode.res[,2]

plot(data_in$Date, data_in$UCI_prevalence, type = 'l', ylim = c(0, 3000))
lines(data_in$Date, fit_d, col = 'green')
data_in$ResidenceTimeICU = 1 / gamma_t

##======================================#
## Save params ---------------
##======================================#
write_csv(data_in %>% dplyr::select(Date, ResidenceTimeICU), '../output/FRED_fit_ICU_residence.csv')

