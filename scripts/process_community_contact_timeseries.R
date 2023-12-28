setwd('/zine/HPC02S1/ex-dveloza/AGORA/apps/fred_colombia_implementation/scripts')
##===============================#
## Process community timeseries
## Author: Guido España
## 2020
##===============================#
## Setup-------------
##===============================#
library(dplyr)
library(tidyverse)
library(lubridate)
library(mgcv)
library(stats)
##===============================#
## Inputs-------------
##===============================#
start_date = '2020-01-01'

##===============================#
## Read data-------------
##===============================#
survey_df = readxl::read_xlsx('../data/20210412_Pregunta_36_tracking.xlsx', skip =4)
colnames(survey_df)[ncol(survey_df)] = 'Date'
survey_df$Date = as.Date(lubridate::parse_date_time(survey_df$Date, locale='es_ES', orders = c("dmy")))
survey_df$Contact_unknown = ifelse(survey_df$P36_E == "SI", 1, 0)
survey_df$Contact_family = ifelse(survey_df$P36_A == "SI", 1, 0)
survey_df$Contact_friends = ifelse(survey_df$P36_B == "SI", 1, 0)

contacts_df = survey_df %>%
    group_by(Date) %>%
    summarize(Contact_unknown = sum(Contact_unknown, na.rm = T) / n(),
              Contact_family = sum(Contact_family, na.rm = T) / n(),
              Contact_friends = sum(Contact_friends, na.rm = T) / n()
              ) %>%
    ungroup() %>%
    mutate(Contacts_total = Contact_unknown + Contact_family + Contact_friends) %>%
    mutate(Date = Date - 10) %>%
    mutate(day = as.numeric(Date - as.Date(start_date)))


mod_gam0 = gam(Contacts_total ~ s(day), data=contacts_df)
mod_splin = smooth.spline(contacts_df$day, contacts_df$Contacts_total)

contacts_model = data.frame(Date = as.Date(seq(from = contacts_df$Date[1], to = contacts_df$Date[nrow(contacts_df)], by = 'day')), stringsAsFactors = F)

contacts_model$day = as.numeric(contacts_model$Date - as.Date(start_date))

contacts_model$Contacts_total = predict(mod_splin, contacts_model$day)$y
contacts_model$Contacts_Prop = contacts_model$Contacts_total / contacts_model$Contacts_total[1]
write_csv(contacts_model, '../experiments_colombia/input_files/community_survey_11001.csv')
