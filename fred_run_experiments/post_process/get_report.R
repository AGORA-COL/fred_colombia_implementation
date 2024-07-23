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

set.seed(123456)
library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(lubridate)
library(ggplot2)
library(sf)
library(raster)
library(rgdal)
library(deSolve)
library(rlang)
library(fda)
library(tidyr)
library(zoo)

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


##===============================#
## Data ------------
##===============================#
args = (commandArgs(TRUE))
calibration_label = 'production'
stage = 'calibration'
projection_label = calibration_label
ss = 27
outdir = '../../../fred_run_experiments/output/CALIBRATION'
asymp_inf = 1.0
fm_ef = 0.73
ksus = 10.0
var_in = 0
vax_in = 70
report_label = 'CHOCO/Grupo_parametros_0'

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



simul_label = sprintf('FRED_%d_%s_asymp_%.2f_fm_%.2f_ksus_%.2f_var_%.0f_vax_%03d_mov_%s_out', ss, stage, asymp_inf,fm_ef, ksus, var_in,vax_in, calibration_label)
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
fred_sweep_df_ = read_csv(data_out) %>% 
    right_join(params_sweep_df, by = c("job_id" = "job_id")) %>%
    mutate(Date = as.Date("2020-01-01") + Day)    
fred_sweep_df_[is.na(fred_sweep_df_)] <- 0 


params_limits = file.path(outdir_st, 'FRED_parameters_limits.csv')
params_limits_df = read_csv(params_limits)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
params_sweep_ll = params_sweep_df %>%
    filter(state_code == ss) %>%
    mutate(LL_total = LL_deaths)

params_sweep_ll <- params_sweep_ll %>% replace_na(list(LL_deaths = 0))
max_LL <- max(params_sweep_ll$LL_deaths)
particles_w <- 1
num_selected_jobs <- 100  # Change this to the desired number

while (TRUE) {
  print(particles_w)
  prob_array <- exp(-(params_sweep_ll$LL_deaths) * particles_w  / max_LL)
  
  if (sum(prob_array) > 0) {
    indx_sampled <- sample.int(n = length(prob_array), size = 3000, prob = prob_array, replace = TRUE)
    n_sampled <- length(unique(indx_sampled))
    
    if (n_sampled <= num_selected_jobs) {
      break
    }
  }
  particles_w <- particles_w + 100
}

selected_jobs <- params_sweep_ll[indx_sampled, ] %>% pull(job_id)
unique(selected_jobs)

fred_sweep_df = fred_sweep_df_ %>% dplyr::filter(job_id %in% selected_jobs)
fred_sweep_df = fred_sweep_df %>% group_by(seed,state_code) %>% mutate(CumCF = cumsum(CF_mean), RR_max = max(RR_mean, na.rm = T)) %>% ungroup()

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

###========================================###
## Refine selection
###========================================###
fred_sweep_df$CF <- fred_sweep_df$CF_mean + 
                    fred_sweep_df$CF_1_mean + 
                    fred_sweep_df$CF_2_mean + 
                    fred_sweep_df$CF_3_mean + 
                    fred_sweep_df$CF_4_mean + 
                    fred_sweep_df$CF_5_mean + 
                    fred_sweep_df$CF_6_mean

################################################################################################################################################
# Reshape the data to wide format
df_wide <- fred_sweep_df %>% 
            dplyr::select(Day, CF, job_id) %>%
            pivot_wider(names_from = job_id, values_from = CF) %>%
            arrange(Day)

# Convert to a matrix for functional data analysis
death_matrix <- as.matrix(df_wide[-1])  # Assuming first column is 'Day' after arrange

# Time points (assuming these are evenly spaced and known)
times <- unique(fred_sweep_df$Day)
dates <- unique(fred_sweep_df$Date)

# Create a basis object for smoothing
nbasis <- 30  # Adjust based on your data
basis_obj <- create.bspline.basis(rangeval = range(times), nbasis = nbasis)

# Convert the matrix to a functional data object
fd <- Data2fd(argvals = times, y = death_matrix, basisobj = basis_obj)

### Smooth curves
smoothed_deaths <- eval.fd(times, fd)
fd_df <- data.frame(Time = times, Deaths = smoothed_deaths)
fd_df_long <- pivot_longer(fd_df, 
                           cols = starts_with("Deaths"), 
                           names_to = "Calibration", 
                           values_to = "Deaths") %>%
              mutate(Calibration = str_replace(Calibration, "Deaths.", ""))

fd_df_long$Day <- fd_df_long$Time + 1
fd_df_long$Date <- as.Date('2020-01-01') + fd_df_long$Time

### Boxplot info
boxplot_output <- boxplot.fd(smoothed_deaths)

outliers_data <- data.frame(
  Calibration = names(boxplot_output$depth)[boxplot_output$outpoint],
  Depth = boxplot_output$depth[boxplot_output$outpoint],
  Outlier = TRUE
)

clean_curves_df <- fd_df_long %>%
  filter(!(Calibration %in% outliers_data$Calibration)) %>%
  mutate(Deaths = ifelse(Deaths < 0, 0, Deaths))

outliers_curves_df <- fd_df_long %>%
  filter((Calibration %in% outliers_data$Calibration)) %>%
  mutate(Deaths = ifelse(Deaths < 0, 0, Deaths))

###========================================###
## Refine selection
###========================================###
fred_sweep_df <- fred_sweep_df %>% dplyr::filter(job_id %in% clean_curves_df$Calibration)
params_sweep_ll_selected = params_sweep_ll %>% dplyr::filter(job_id %in% clean_curves_df$Calibration)
write_csv(params_sweep_ll_selected, paste0("./output/", report_label, '/parameters/FRED_parameters_out_selected.csv'))

###========================================###
## Process selected
###========================================###                                         
state_fred = fred_sweep_df %>%
    group_by(job_id) %>%
    mutate(AR_mean = 100 * (cumsum(C_mean + C_1_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean - reinfections)) / N_mean[1],
           AR_0_mean = 100 * (cumsum(C_mean)) / N_mean[1],
           AR_1_mean = 100 * (cumsum(C_1_mean)) / N_mean[1],
           AR_2_mean = 100 * (cumsum(C_2_mean)) / N_mean[1],
           AR_3_mean = 100 * (cumsum(C_3_mean)) / N_mean[1],
           AR_4_mean = 100 * (cumsum(C_4_mean)) / N_mean[1],
           AR_5_mean = 100 * (cumsum(C_5_mean)) / N_mean[1]) %>%
    ungroup()


tmp_fred_all = state_fred %>%
    group_by(intervention_id, Day, Date, start_date, school_closure_day) %>%
    summarize(CF_median = quantile(CF_mean + CF_1_mean + CF_2_mean + CF_3_mean + CF_4_mean + CF_5_mean, probs = c(0.5), na.rm = T),
              CF_low = quantile(CF_mean + CF_1_mean + CF_2_mean + CF_3_mean + CF_4_mean + CF_5_mean, probs = c(0.0), na.rm = T),
              CF_high = quantile(CF_mean + CF_1_mean + CF_2_mean + CF_3_mean +  CF_4_mean + CF_5_mean, probs = c(1.0), na.rm = T),
              
              Cs_median = quantile(Cs_mean + Cs_1_mean + Cs_2_mean + Cs_3_mean + Cs_4_mean + Cs_5_mean, probs = c(0.5), na.rm = T),
              Cs_low = quantile(Cs_mean + Cs_1_mean + Cs_2_mean + Cs_3_mean + Cs_4_mean + Cs_5_mean, probs = c(0.0), na.rm = T),
              Cs_high = quantile(Cs_mean + Cs_1_mean + Cs_2_mean + Cs_3_mean + Cs_4_mean + Cs_5_mean, probs = c(1.0), na.rm = T),
              
              reinfections_median = quantile(reinfections / (C_1_mean + C_mean+ C_2_mean + C_3_mean + C_4_mean + C_5_mean), probs = c(0.5), na.rm = T),
              reinfections_low = quantile(reinfections / (C_1_mean + C_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean), probs = c(0.0), na.rm = T),
              reinfections_high = quantile(reinfections / (C_1_mean + C_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean), probs = c(1.0), na.rm = T),
              C_median = quantile(C_mean + C_1_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean, probs = c(0.5), na.rm = T),
              C_low = quantile(C_mean + C_1_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean, probs = c(0.0), na.rm = T),
              C_high = quantile(C_mean + C_1_mean+ C_2_mean + C_3_mean + C_4_mean + C_5_mean, probs = c(1.0), na.rm = T),
              
              C_0_median = quantile(C_mean, probs = c(0.5), na.rm = T),
              C_0_low = quantile(C_mean, probs = c(0.0), na.rm = T),
              C_0_high = quantile(C_mean, probs = c(1.0), na.rm = T),
              
              C_1_median = quantile(C_1_mean, probs = c(0.5), na.rm = T),
              C_1_low = quantile(C_1_mean, probs = c(0.0), na.rm = T),
              C_1_high = quantile(C_1_mean, probs = c(1.0), na.rm = T),
              
              C_2_median = quantile(C_2_mean, probs = c(0.5), na.rm = T),
              C_2_low = quantile(C_2_mean, probs = c(0.0), na.rm = T),
              C_2_high = quantile(C_2_mean, probs = c(1.0), na.rm = T),
              
              C_3_median = quantile(C_3_mean, probs = c(0.5), na.rm = T),
              C_3_low = quantile(C_3_mean, probs = c(0.0), na.rm = T),
              C_3_high = quantile(C_3_mean, probs = c(1.0), na.rm = T),
              
              C_4_median = quantile(C_4_mean, probs = c(0.5), na.rm = T),
              C_4_low = quantile(C_4_mean, probs = c(0.0), na.rm = T),
              C_4_high = quantile(C_4_mean, probs = c(1.0), na.rm = T),
              
              C_5_median = quantile(C_5_mean, probs = c(0.5), na.rm = T),
              C_5_low = quantile(C_5_mean, probs = c(0.0), na.rm = T),
              C_5_high = quantile(C_5_mean, probs = c(1.0), na.rm = T),
              
              CF_0_median = quantile(CF_mean, probs = c(0.5), na.rm = T),
              CF_0_low = quantile(CF_mean, probs = c(0.0), na.rm = T),
              CF_0_high = quantile(CF_mean, probs = c(1.0), na.rm = T),
              
              CF_1_median = quantile(CF_1_mean, probs = c(0.5), na.rm = T),
              CF_1_low = quantile(CF_1_mean, probs = c(0.0), na.rm = T),
              CF_1_high = quantile(CF_1_mean, probs = c(1.0), na.rm = T),
              
              CF_2_median = quantile(CF_2_mean, probs = c(0.5), na.rm = T),
              CF_2_low = quantile(CF_2_mean, probs = c(0.0), na.rm = T),
              CF_2_high = quantile(CF_2_mean, probs = c(1.0), na.rm = T),
              
              CF_3_median = quantile(CF_3_mean, probs = c(0.5), na.rm = T),
              CF_3_low = quantile(CF_3_mean, probs = c(0.0), na.rm = T),
              CF_3_high = quantile(CF_3_mean, probs = c(1.0), na.rm = T),
              
              CF_4_median = quantile(CF_4_mean, probs = c(0.5), na.rm = T),
              CF_4_low = quantile(CF_4_mean, probs = c(0.0), na.rm = T),
              CF_4_high = quantile(CF_4_mean, probs = c(1.0), na.rm = T),
              
              CF_5_median = quantile(CF_5_mean, probs = c(0.5), na.rm = T),
              CF_5_low = quantile(CF_5_mean, probs = c(0.0), na.rm = T),
              CF_5_high = quantile(CF_5_mean, probs = c(1.0), na.rm = T),

              ACF0_5_median = quantile(ACF0_5, probs = c(0.5), na.rm = TRUE),
              ACF0_5_low = quantile(ACF0_5, probs = c(0.0), na.rm = TRUE),
              ACF0_5_high = quantile(ACF0_5, probs = c(1.0), na.rm = TRUE),
              
              ACF5_10_median = quantile(ACF5_10, probs = c(0.5), na.rm = TRUE),
              ACF5_10_low = quantile(ACF5_10, probs = c(0.0), na.rm = TRUE),
              ACF5_10_high = quantile(ACF5_10, probs = c(1.0), na.rm = TRUE),
              
              ACF10_15_median = quantile(ACF10_15, probs = c(0.5), na.rm = TRUE),
              ACF10_15_low = quantile(ACF10_15, probs = c(0.0), na.rm = TRUE),
              ACF10_15_high = quantile(ACF10_15, probs = c(1.0), na.rm = TRUE),
              
              ACF15_20_median = quantile(ACF15_20, probs = c(0.5), na.rm = TRUE),
              ACF15_20_low = quantile(ACF15_20, probs = c(0.0), na.rm = TRUE),
              ACF15_20_high = quantile(ACF15_20, probs = c(1.0), na.rm = TRUE),
              
              ACF20_25_median = quantile(ACF20_25, probs = c(0.5), na.rm = TRUE),
              ACF20_25_low = quantile(ACF20_25, probs = c(0.0), na.rm = TRUE),
              ACF20_25_high = quantile(ACF20_25, probs = c(1.0), na.rm = TRUE),
              
              ACF25_30_median = quantile(ACF25_30, probs = c(0.5), na.rm = TRUE),
              ACF25_30_low = quantile(ACF25_30, probs = c(0.0), na.rm = TRUE),
              ACF25_30_high = quantile(ACF25_30, probs = c(1.0), na.rm = TRUE),
              
              ACF30_35_median = quantile(ACF30_35, probs = c(0.5), na.rm = TRUE),
              ACF30_35_low = quantile(ACF30_35, probs = c(0.0), na.rm = TRUE),
              ACF30_35_high = quantile(ACF30_35, probs = c(1.0), na.rm = TRUE),
              
              ACF35_40_median = quantile(ACF35_40, probs = c(0.5), na.rm = TRUE),
              ACF35_40_low = quantile(ACF35_40, probs = c(0.0), na.rm = TRUE),
              ACF35_40_high = quantile(ACF35_40, probs = c(1.0), na.rm = TRUE),
              
              ACF40_45_median = quantile(ACF40_45, probs = c(0.5), na.rm = TRUE),
              ACF40_45_low = quantile(ACF40_45, probs = c(0.0), na.rm = TRUE),
              ACF40_45_high = quantile(ACF40_45, probs = c(1.0), na.rm = TRUE),
              
              ACF45_50_median = quantile(ACF45_50, probs = c(0.5), na.rm = TRUE),
              ACF45_50_low = quantile(ACF45_50, probs = c(0.0), na.rm = TRUE),
              ACF45_50_high = quantile(ACF45_50, probs = c(1.0), na.rm = TRUE),
              
              ACF50_55_median = quantile(ACF50_55, probs = c(0.5), na.rm = TRUE),
              ACF50_55_low = quantile(ACF50_55, probs = c(0.0), na.rm = TRUE),
              ACF50_55_high = quantile(ACF50_55, probs = c(1.0), na.rm = TRUE),
              
              ACF55_60_median = quantile(ACF55_60, probs = c(0.5), na.rm = TRUE),
              ACF55_60_low = quantile(ACF55_60, probs = c(0.0), na.rm = TRUE),
              ACF55_60_high = quantile(ACF55_60, probs = c(1.0), na.rm = TRUE),
              
              ACF60_65_median = quantile(ACF60_65, probs = c(0.5), na.rm = TRUE),
              ACF60_65_low = quantile(ACF60_65, probs = c(0.0), na.rm = TRUE),
              ACF60_65_high = quantile(ACF60_65, probs = c(1.0), na.rm = TRUE),
              
              ACF65_70_median = quantile(ACF65_70, probs = c(0.5), na.rm = TRUE),
              ACF65_70_low = quantile(ACF65_70, probs = c(0.0), na.rm = TRUE),
              ACF65_70_high = quantile(ACF65_70, probs = c(1.0), na.rm = TRUE),
              
              ACF70_75_median = quantile(ACF70_75, probs = c(0.5), na.rm = TRUE),
              ACF70_75_low = quantile(ACF70_75, probs = c(0.0), na.rm = TRUE),
              ACF70_75_high = quantile(ACF70_75, probs = c(1.0), na.rm = TRUE),
              
              ACF75_80_median = quantile(ACF75_80, probs = c(0.5), na.rm = TRUE),
              ACF75_80_low = quantile(ACF75_80, probs = c(0.0), na.rm = TRUE),
              ACF75_80_high = quantile(ACF75_80, probs = c(1.0), na.rm = TRUE),
              
              ACF80_120_median = quantile(ACF80_120, probs = c(0.5), na.rm = TRUE),
              ACF80_120_low = quantile(ACF80_120, probs = c(0.0), na.rm = TRUE),
              ACF80_120_high = quantile(ACF80_120, probs = c(1.0), na.rm = TRUE),
              
              #V_median = quantile(V_mean, probs = c(0.5), na.rm = T),
              # V_low = quantile(V_mean, probs = c(0.0), na.rm = T),
              # V_high = quantile(V_mean, probs = c(1.0), na.rm = T),
              # Vtd_median = quantile(Vtd_mean, probs = c(0.5), na.rm = T),
              # Vtd_low = quantile(Vtd_mean, probs = c(0.0), na.rm = T),
              # Vtd_high = quantile(Vtd_mean, probs = c(1.0), na.rm = T),
              
              DomOriginal_median = quantile(C_mean / (C_1_mean + C_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean), probs = c(0.5), na.rm = T),
              DomOriginal_low = quantile(C_mean / (C_1_mean + C_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean), probs = c(0.0), na.rm = T),
              DomOriginal_high = quantile(C_mean / (C_1_mean + C_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean), probs = c(1.0), na.rm = T),
              DomAlpha_median = quantile(C_1_mean / (C_1_mean + C_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean), probs = c(0.5), na.rm = T),
              DomAlpha_low = quantile(C_1_mean / (C_1_mean + C_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean), probs = c(0.0), na.rm = T),
              DomAlpha_high = quantile(C_1_mean /(C_1_mean + C_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean), probs = c(1.0), na.rm = T),
              DomGamma_median = quantile(C_2_mean / (C_1_mean + C_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean), probs = c(0.5), na.rm = T),
              DomGamma_low = quantile(C_2_mean / (C_1_mean + C_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean), probs = c(0.0), na.rm = T),
              DomGamma_high = quantile(C_2_mean /(C_1_mean + C_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean), probs = c(1.0), na.rm = T),
              DomKappa_median = quantile(C_3_mean / (C_1_mean + C_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean), probs = c(0.5), na.rm = T),
              DomKappa_low = quantile(C_3_mean / (C_1_mean + C_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean), probs = c(0.0), na.rm = T),
              DomKappa_high = quantile(C_3_mean /(C_1_mean + C_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean), probs = c(1.0), na.rm = T),
              DomDelta_median = quantile(C_4_mean / (C_1_mean + C_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean), probs = c(0.5), na.rm = T),
              DomDelta_low = quantile(C_4_mean / (C_1_mean + C_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean), probs = c(0.0), na.rm = T),
              DomDelta_high = quantile(C_4_mean /(C_1_mean + C_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean), probs = c(1.0), na.rm = T),
              DomOmicron_median = quantile(C_5_mean / (C_1_mean + C_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean), probs = c(0.5), na.rm = T),
              DomOmicron_low = quantile(C_5_mean / (C_1_mean + C_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean), probs = c(0.0), na.rm = T),
              DomOmicron_high = quantile(C_5_mean /(C_1_mean + C_mean + C_2_mean + C_3_mean + C_4_mean + C_5_mean), probs = c(1.0), na.rm = T),
              Phosp_median = 0,
              Phosp_low = 0,
              Phosp_high = 0,

              Chosp_median = quantile(Chosp_mean + Chosp_1_mean + Chosp_2_mean +  Chosp_3_mean  + Chosp_4_mean + Chosp_5_mean, probs = c(0.5), na.rm = T),
              Chosp_low = quantile(Chosp_mean + Chosp_1_mean + Chosp_2_mean + Chosp_3_mean + Chosp_4_mean + Chosp_5_mean, probs = c(0.0), na.rm = T),
              Chosp_high = quantile(Chosp_mean + Chosp_1_mean + Chosp_2_mean + Chosp_3_mean + Chosp_4_mean + Chosp_5_mean, probs = c(1.0), na.rm = T),
              #PrevInf_median = quantile(PrevInf_mean, probs = c(0.5), na.rm = T),
              #PrevInf_low = quantile(PrevInf_mean, probs = c(0.0), na.rm = T),
              #PrevInf_high = quantile(PrevInf_mean, probs = c(1.0), na.rm = T),
              AR_median = medCI(AR_mean)/100,
              AR_low = lowCI(AR_mean)/100,
              AR_high = highCI(AR_mean)/100,
              #Sch_median = mean(Sch_mean),
              #Wrk_median = mean(Wrk_mean),
              #Nbr_median = mean(Nbr_mean),
              # H_median = medCI(H_mean),
              # N_median = medCI(N_mean),
              # Nursing_home_CF_median = medCI(Nursing_Home_CF_mean),
              # N_sheltering_median = 1 -  medCI(N_sheltering_mean/N_mean[1]),
              # N_sheltering_low = 1-  lowCI(N_sheltering_mean/N_mean[1]),
              # N_sheltering_high = 1 -  highCI(N_sheltering_mean/N_mean[1]),
              
              AR_0_low = quantile(AR_0_mean, probs = c(0.0), na.rm = T),
              AR_0_high = quantile(AR_0_mean, probs = c(1.0), na.rm = T),
              AR_0_median = quantile(AR_0_mean, probs = c(0.5), na.rm = T),
              
              AR_1_low = quantile(AR_1_mean, probs = c(0.0), na.rm = T),
              AR_1_high = quantile(AR_1_mean, probs = c(1.0), na.rm = T),
              AR_1_median = quantile(AR_1_mean, probs = c(0.5), na.rm = T),
              
              AR_2_low = quantile(AR_2_mean, probs = c(0.0), na.rm = T),
              AR_2_high = quantile(AR_2_mean, probs = c(1.0), na.rm = T),
              AR_2_median = quantile(AR_2_mean, probs = c(0.5), na.rm = T),
              
              AR_3_low = quantile(AR_3_mean, probs = c(0.0), na.rm = T),
              AR_3_high = quantile(AR_3_mean, probs = c(1.0), na.rm = T),
              AR_3_median = quantile(AR_3_mean, probs = c(0.5), na.rm = T),
              
              AR_4_low = quantile(AR_4_mean, probs = c(0.0), na.rm = T),
              AR_4_high = quantile(AR_4_mean, probs = c(1.0), na.rm = T),
              AR_4_median = quantile(AR_4_mean, probs = c(0.5), na.rm = T),
              
              AR_5_low = quantile(AR_5_mean, probs = c(0.0), na.rm = T),
              AR_5_high = quantile(AR_5_mean, probs = c(1.0), na.rm = T),
              AR_5_median = quantile(AR_5_mean, probs = c(0.5), na.rm = T),
              
              RR_median = medCI(RR_mean) * (1 - (DomAlpha_median + DomGamma_median + DomKappa_median + DomDelta_median + DomOmicron_median)),
              RR_1_median = medCI(RR_1_mean) * DomAlpha_median,
              RR_2_median = medCI(RR_2_mean) * DomGamma_median,
              RR_3_median = medCI(RR_3_mean) * DomKappa_median,
              RR_4_median = medCI(RR_4_mean) * DomDelta_median,   
              RR_5_median = medCI(RR_5_mean) * DomOmicron_median,   
              shelter_in_place_delay_mean = mean(shelter_in_place_delay_mean))%>%                               
    ungroup() %>%
    mutate(RR_total = RR_median + RR_1_median + RR_2_median + RR_3_median + RR_4_median + RR_5_median)


## adjust uci occupancy -----------

print(colnames(tmp_fred_all))
## Read fitted hosp parameters
#icu_duration = read_csv(sprintf('../../../fit_scripts/fitted_data/fred_fit_model_%s_var_%s_%s.csv', ss, var_in, calibration_label))
icu_duration = read_csv(paste0("./output/", report_label, '/fitted_data/FRED_fit_ICU_residence.csv'))

uci_df_dept = uci_df_dept %>% rename('UCI_available' = 'total_beds')#, 'UCI_prevalence' = 'prevalence'

fred_adjusted = tibble()
for(in_id in unique(tmp_fred_all$intervention_id)){
    tmp_df = tmp_fred_all %>% filter(intervention_id == in_id) %>%
        left_join(icu_duration, by = 'Date') %>%
        left_join(uci_df_dept[,c('Date', 'UCI_available')], by = 'Date') %>%
        replace_na(list(ResidenceTimeICU = 11))
    tmp_df$UCI_available[is.na(tmp_df$UCI_available)] = max(tmp_df$UCI_available, na.rm = T)
    ## replace_na(list(ResidenceTimeICU = uci_df_dept$ResidenceTimeICU[nrow(uci_df_dept)]))

    print(colnames(tmp_df))
    
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
                                                 
write_csv(tmp_fred_all, paste0("./output/", report_label, sprintf('/fitted_data/fred_output_model_%s_var_%s_%s.csv', ss, var_in, projection_label)))
