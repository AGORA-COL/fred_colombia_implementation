##===============================#
## Plot model fit for Bogotá
## 2020
####################################################
## local variables
####################################################
repo_name = 'fred_colombia_implementation'
AGORA_path = '/zine/HPC02S1/ex-dveloza/AGORA/apps'

####################################################
## FRED set cluster enviromental variables
####################################################
setwd(sprintf('%s/%s/fred_run_experiments/post_process/model-reports', AGORA_path, repo_name))

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
args = (commandArgs(TRUE))
calibration_label = 'production'
projection_label = calibration_label
ss = 27
outdir = '../../../fred_run_experiments/output/CALIBRATION'
asymp_inf = 1.0
fm_ef = 0.73
ksus = 10.0
stage = 'calibration'
var_in = 1
vax_in = 70
report_label = 'CHOCO/Grupo_parametros_1'

if(length(args) >= 1){
    calibration_label = args[1]
    if(length(args) >= 2){
        ss = as.numeric(args[2])
        if(length(args) >= 3){
            outdir = args[3]
            if(length(args) >= 4){
                asymp_inf = as.numeric(args[4])
                if(length(args) >= 5){
                    fm_ef = as.numeric(args[5])
                    if(length(args) >= 6){
                        ksus = as.numeric(args[6])
                        if(length(args) >= 7){
                            stage = args[7]
                            if(length(args) >= 8){
                                var_in = as.numeric(args[8])
                                if(length(args) >= 9){
                                    vax_in = as.numeric(args[9])
                                    if(length(args) >= 10){
                                        report_label = args[10]
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

simul_label = sprintf('FRED_%d_%s_asymp_%.2f_fm_%.2f_ksus_%.2f_var_%.0f_vax_%03d_mov_%s_out', ss, stage, asymp_inf,fm_ef, ksus, var_in, vax_in, calibration_label)
outdir_st = file.path(outdir, simul_label)
interventions_df = read_csv('../../../fred_input_files/interventions/interventions_Colombia.csv') 

brk_ages = c(0, seq(from=10,by = 10, to = 80), 120)
brk_lbls = sprintf("ACF%d_%d", brk_ages[-length(brk_ages)], brk_ages[-1])
brk_data_lbls = sprintf("A%d_%d", brk_ages[-length(brk_ages)], brk_ages[-1])
base_age_df = data.frame(AgeGroup = brk_lbls, Deaths = 0, stringsAsFactors = F)
base_age_data = data.frame(AgeGroup = brk_data_lbls, stringsAsFactors = F)

fit_date = as.Date("2021-07-31")
brk_ages = c(0, seq(from=10,by = 10, to = 80), 120)
brk_lbls = sprintf("ACF%d_%d", brk_ages[-length(brk_ages)], brk_ages[-1])
brk_data_lbls = sprintf("A%d_%d", brk_ages[-length(brk_ages)], brk_ages[-1])
base_age_df = data.frame(AgeGroup = brk_lbls, Deaths = 0, stringsAsFactors = F)
base_age_data = data.frame(AgeGroup = brk_data_lbls, stringsAsFactors = F)

uci_df_dept <- read.csv('../../../process_inputs/covid_uci_beds/output_files/camas_uci_departamentos_colombia.csv') %>% 
              dplyr::filter(capacity == "Cuidado Intensivo Adulto") %>%
              mutate(Date = as.Date(date), dept_code = as.numeric(gsub(".*?(\\d+).*", "\\1", department))) %>%
              dplyr::filter(dept_code == ss) %>%
              mutate(UCI_prevalence = covid_confirmed + covid_suspected) %>%
              dplyr::select(Date, UCI_prevalence, total_beds)
              
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
if(stage == 'calibration'){
    fred_sweep_df$intervention_id = 'calibration'
}

# List of expected base column names
base_columns <- c("C", "CF", "Cs", "Chosp", "RR")

# List of variants to apply to each base column
variants <- paste0("_", 1:6)

# Generate the expected column names with variants and '_mean' suffix
column_variants <- c(sapply(base_columns, function(col) {
    paste0(col, variants, "_mean")
}), "reinfections")

# Initialize missing columns with zero
fred_sweep_df <- as_tibble(fred_sweep_df)  # Ensure fred_sweep_df is a tibble
for (col in column_variants) {
  if (!col %in% names(fred_sweep_df)) {
    fred_sweep_df[[col]] <- 0
  }
}

state_fred = fred_sweep_df %>% drop_na() %>%
    group_by(job_id) %>%
    mutate(AR_mean = 100 * (cumsum(C_mean) + cumsum(C_1_mean)) / N_mean[1]) %>%
    ungroup()

tmp_fred_all = state_fred %>%
    group_by(intervention_id, Day, Date, start_date, school_closure_day) %>%
    summarize(CF_median = quantile(CF_mean + CF_1_mean + CF_2_mean + CF_3_mean + CF_4_mean + CF_5_mean + CF_6_mean, probs = c(0.5), na.rm = T),
              CF_low = quantile(CF_mean + CF_1_mean + CF_2_mean + CF_3_mean + CF_4_mean + CF_5_mean + CF_6_mean, probs = c(0.025), na.rm = T),
              CF_high = quantile(CF_mean + CF_1_mean + CF_2_mean + CF_3_mean + CF_4_mean + CF_5_mean + CF_6_mean, probs = c(0.975), na.rm = T),
              Cs_median = quantile(Cs_mean + Cs_1_mean + Cs_2_mean + Cs_3_mean + Cs_4_mean + Cs_5_mean + Cs_6_mean, probs = c(0.5), na.rm = T),
              Cs_low = quantile(Cs_mean + Cs_1_mean + Cs_2_mean + Cs_3_mean + Cs_4_mean + Cs_5_mean + Cs_6_mean, probs = c(0.025), na.rm = T),
              Cs_high = quantile(Cs_mean + Cs_1_mean + Cs_2_mean + Cs_3_mean + Cs_4_mean + Cs_5_mean + Cs_6_mean, probs = c(0.975), na.rm = T),
              reinfections_median = quantile(reinfections / (C_1_mean + C_mean + Cs_2_mean + Cs_3_mean + Cs_4_mean + Cs_5_mean + Cs_6_mean), probs = c(0.5), na.rm = T),
              reinfections_low = quantile(reinfections / (C_1_mean + C_mean + Cs_2_mean + Cs_3_mean + Cs_4_mean + Cs_5_mean + Cs_6_mean), probs = c(0.025), na.rm = T),
              reinfections_high = quantile(reinfections / (C_1_mean + C_mean + Cs_2_mean + Cs_3_mean + Cs_4_mean + Cs_5_mean + Cs_6_mean), probs = c(0.975), na.rm = T),
              C_median = quantile(C_mean + C_1_mean + Cs_2_mean + Cs_3_mean + Cs_4_mean + Cs_5_mean + Cs_6_mean, probs = c(0.5), na.rm = T),
              C_low = quantile(C_mean + C_1_mean + Cs_2_mean + Cs_3_mean + Cs_4_mean + Cs_5_mean + Cs_6_mean, probs = c(0.025), na.rm = T),
              C_high = quantile(C_mean + C_1_mean + Cs_2_mean + Cs_3_mean + Cs_4_mean + Cs_5_mean + Cs_6_mean, probs = c(0.975), na.rm = T),
              C_1_median = quantile(C_1_mean + Cs_2_mean + Cs_3_mean + Cs_4_mean + Cs_5_mean + Cs_6_mean, probs = c(0.5), na.rm = T),
              C_1_low = quantile(C_1_mean + Cs_2_mean + Cs_3_mean + Cs_4_mean + Cs_5_mean + Cs_6_mean, probs = c(0.025), na.rm = T),
              C_1_high = quantile(C_1_mean + Cs_2_mean + Cs_3_mean + Cs_4_mean + Cs_5_mean + Cs_6_mean, probs = c(0.975), na.rm = T),
              DomVariant_median = quantile(C_1_mean / (C_1_mean + C_mean + Cs_2_mean + Cs_3_mean + Cs_4_mean + Cs_5_mean + Cs_6_mean), probs = c(0.5), na.rm = T),
              DomVariant_low = quantile(C_1_mean / (C_mean + C_1_mean + Cs_2_mean + Cs_3_mean + Cs_4_mean + Cs_5_mean + Cs_6_mean), probs = c(0.025), na.rm = T),
              DomVariant_high = quantile(C_1_mean / (C_mean + C_1_mean + Cs_2_mean + Cs_3_mean + Cs_4_mean + Cs_5_mean + Cs_6_mean), probs = c(0.975), na.rm = T),
              #Phosp_median = quantile(Phosp_mean + Phosp_1_mean, probs = c(0.5), na.rm = T),
              #Phosp_low = quantile(Phosp_mean + Phosp_1_mean, probs = c(0.025), na.rm = T),
              #Phosp_high = quantile(Phosp_mean + Phosp_1_mean, probs = c(0.975), na.rm = T),
              Chosp_median = quantile(Chosp_mean + Chosp_1_mean + Chosp_2_mean + Chosp_3_mean + Chosp_4_mean + Chosp_5_mean + Chosp_6_mean, probs = c(0.5), na.rm = T),
              Chosp_low = quantile(Chosp_mean + Chosp_1_mean + Chosp_2_mean + Chosp_3_mean + Chosp_4_mean + Chosp_5_mean + Chosp_6_mean, probs = c(0.025), na.rm = T),
              Chosp_high = quantile(Chosp_mean + Chosp_1_mean + Chosp_2_mean + Chosp_3_mean + Chosp_4_mean + Chosp_5_mean + Chosp_6_mean, probs = c(0.975), na.rm = T),
              #PrevInf_median = quantile(PrevInf_mean, probs = c(0.5), na.rm = T),
              #PrevInf_low = quantile(PrevInf_mean, probs = c(0.025), na.rm = T),
              #PrevInf_high = quantile(PrevInf_mean, probs = c(0.975), na.rm = T),
              AR_median = medCI(AR_mean)/100,
              AR_low = lowCI(AR_mean)/100,
              AR_high = highCI(AR_mean)/100,
              #Sch_median = mean(Sch_mean),
            #   Wrk_median = mean(Wrk_mean),
            #   Nbr_median = mean(Nbr_mean),
            #   H_median = medCI(H_mean),
            #   N_median = medCI(N_mean),
            #   Nursing_home_CF_median = medCI(Nursing_Home_CF_mean),
            #   N_sheltering_median = 1 -  medCI(N_sheltering_mean/N_mean[1]),
            #   N_sheltering_low = 1-  lowCI(N_sheltering_mean/N_mean[1]),
            #   N_sheltering_high = 1 -  highCI(N_sheltering_mean/N_mean[1]),
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
names(uci_df_dept)

data_in = left_join(tmp_fred_all, uci_df_dept, by = 'Date') %>%
  mutate(across(everything(), ~replace_na(., 0)))
##data_in = left_join(uci_df, tmp_fred_mov, by = 'Date')
##data_in = left_join(uci_df, tmp_fred_mov_vax, by = 'Date')
##data_in = left_join(uci_df, uci_incidence, by = 'Date')

ch_time = data_in$Chosp_median

tt = seq(from=0,to = nrow(data_in)-1, by = 1)

init_state = c(h = 0)
nodes_bs = round(nrow(data_in)/50) #100
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
write.csv(data_in %>% dplyr::select(Date, ResidenceTimeICU), paste0("./output/", report_label, '/fitted_data/FRED_fit_ICU_residence.csv'))
