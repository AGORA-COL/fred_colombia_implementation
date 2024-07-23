####################################################
## local variables
####################################################
repo_name = 'fred_colombia_implementation'
AGORA_path = '/zine/HPC02S1/ex-dveloza/AGORA/apps'

####################################################
## FRED set cluster enviromental variables
####################################################
setwd(sprintf('%s/%s/fred_run_experiments/post_process', AGORA_path, repo_name))

#Load packages
library(tidyverse)
library(dplyr)
library(ggplot2)
library(EpiEstim)


args = (commandArgs(TRUE))
calibration_label = 'production'
projection_label = calibration_label
departament_code = 27
outdir = '../../fred_run_experiments/output/CALIBRATION'
asymp_inf = 1.0
fm_ef = 0.73
ksus = 10.0
stage = 'calibration'
var_in = 0
vax_in = 70

if(length(args) >= 1){
    calibration_label = args[1]
    if(length(args) >= 2){
        departament_code = as.numeric(args[2])
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
                                }
                            }
                        }                        
                    }
                }
            }
        }
    }
}
##===============================#
## Setup output-------------
##===============================#
geoinfo_file = '../../fred_input_files/geoinfo_municipios_colombia.csv'
df_geo_info = read.csv(geoinfo_file) %>% 
              dplyr::select(COD_DEPTO, NOM_DEPART) %>% 
              distinct() %>%
              mutate(NOM_DEPART = str_replace_all(NOM_DEPART, "[ÁÄÂÀÃÅ]", "A")) %>%
              mutate(NOM_DEPART = str_replace_all(NOM_DEPART, "[ÉËÊÈ]", "E")) %>%
              mutate(NOM_DEPART = str_replace_all(NOM_DEPART, "[ÍÏÎÌ]", "I")) %>%
              mutate(NOM_DEPART = str_replace_all(NOM_DEPART, "[ÓÖÔÒÕ]", "O")) %>%
              mutate(NOM_DEPART = str_replace_all(NOM_DEPART, "[ÚÜÛÙ]", "U")) %>%
              mutate(NOM_DEPART = str_replace_all(NOM_DEPART, "[Ñ]", "N")) %>%
              mutate(NOM_DEPART = str_replace_all(NOM_DEPART, "[ .,]", "_"))

##===============================#
## Setup output-------------
##===============================#
dept_name <- df_geo_info %>% filter(COD_DEPTO == departament_code) %>% pull(NOM_DEPART)
report_label = sprintf('%s_%s/grupo_parametros_%s', departament_code, dept_name, var_in)
outdir_report <- file.path("model-reports/input_files", report_label)

if(!dir.exists(outdir_report)){
    dir.create(outdir_report, recursive = TRUE, showWarnings = FALSE)
}

deaths_file            <- '../../fred_input_files/covid_data/COL_covid_death_data.csv'
age_group_deaths_files <- '../../fred_input_files/covid_data/Age_COL_covid_data.csv'


df_age_deaths <- read.csv(age_group_deaths_files) %>% 
                  mutate(DeptCode = substr(MunCode, 1, nchar(MunCode) - 3)) %>%
                  group_by(Date, DeptCode, AgeGroup) %>%
                  summarise(Deaths = sum(Deaths)) %>% ungroup()

df_age_dept_data <- df_age_deaths %>% filter(DeptCode == departament_code)

df_deaths <- read.csv(deaths_file) %>% 
            mutate(DeptCode = substr(MunCode, 1, nchar(MunCode) - 3)) %>%
            group_by(Date, DeptCode) %>%
            summarise(Cases = sum(Cases), Deaths = sum(Deaths)) %>% ungroup()

df_dept_data <- df_deaths %>% filter(DeptCode == departament_code)

#Load data
daily_cases <- df_dept_data %>% dplyr::select(Date, Cases) %>% rename("date" = "Date", "cases" = "Cases")
daily_cases$date <- as.Date(daily_cases$date)

all_dates <- data.frame(
  date = seq.Date(
    min(daily_cases$date, na.rm = TRUE),
    max(daily_cases$date, na.rm = TRUE),
    by = "day",
    origin = "1970-01-01" # Add the origin parameter
  )
)

#Parameters for COVID-19
rt_window <- 14
incubation_period <- 5
mean_si <- 5.2
std_si <- 0.3

#`compute_rt()` uses `EpiEstim::estimate_R()` to compute $R_t$ from a given incidence data frame
compute_rt <- function(
    df_incidence,
    method = "parametric_si",
    mean_si, # Mean serial interval
    std_si, # Standard deviation of the serial interval
    rt_window, # Time window length
    incubation_period # Incubation period
) {
  t_start <- seq(incubation_period, nrow(df_incidence) - rt_window)
  t_end <- t_start + rt_window
  
  rt_data <- estimate_R(
    df_incidence,
    method = method,
    config = make_config(
      list(
        mean_si = mean_si,
        std_si = std_si,
        t_start = t_start,
        t_end = t_end
      )
    )
  )
  
  df_rt <- rt_data$R
  df_rt$window_start <- min(df_incidence$infection) + df_rt$t_start
  df_rt$window_end <- min(df_incidence$infection) + df_rt$t_end
  
  return(df_rt)
}

daily_cases$date <- as.Date(daily_cases$date) 
daily_cases$infection <- daily_cases$date - incubation_period

daily_cases$infection <- as.Date(daily_cases$infection) 
daily_cases$infection <- ymd(daily_cases$infection)


department_rt <- daily_cases %>% 
          rename(I = cases) %>% 
          left_join(all_dates, by = "date") %>% 
          replace(is.na(df), 0) %>% 
    compute_rt(
    rt_window = rt_window,
    mean_si = mean_si,
    std_si = std_si,
    incubation_period = incubation_period
  ) %>% 
  select(c('window_start', 'window_end', 'Mean(R)', 'Quantile.0.05(R)', 'Quantile.0.975(R)')) %>%
  rename(
    rt_mean = `Mean(R)`,
    rt_mean_lower = `Quantile.0.05(R)`,
    rt_mean_upper = `Quantile.0.975(R)`
  )


write.csv(department_rt, sprintf('%s/rt_reported.csv', outdir_report))
write.csv(df_dept_data, sprintf('%s/%s_deaths.csv', outdir_report, departament_code))
write.csv(df_age_dept_data, sprintf('%s/%s_age_groups_deaths.csv', outdir_report, departament_code))