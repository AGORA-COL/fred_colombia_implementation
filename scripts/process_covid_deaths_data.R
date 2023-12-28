setwd('/zine/HPC02S1/ex-dveloza/AGORA/apps/fred_colombia_implementation/scripts')
##===============================#
## Author: Diego Veloza
## 21.12.2023
##===============================#
## Setup-------------
##===============================#
set.seed(123456)
library(dplyr)
library(tidyverse)


colombia_cases_file = '../data/Covid_Colombia_cases_20231207.csv'
df_raw_data = read.csv(colombia_cases_file)

df_ = df_raw_data %>% 
    dplyr::filter(Estado == 'Fallecido') %>%
    mutate(MunCode = Código.DIVIPOLA.municipio, 
           Date = as.Date(Fecha.de.muerte),
           AgeUnits = Unidad.de.medida.de.edad,
           Age = Edad) %>% 
    dplyr::select(MunCode, Date, AgeUnits, Age)


df_$Age <- ifelse(df_$AgeUnits == 2, df_$Age / 12, 
                 ifelse(df_$AgeUnits == 3, df_$Age / 365, df_$Age))

age_breaks <- c(0, seq(from=5,by = 5, to = 80), 120)
age_labels <- paste0("ACF", head(age_breaks, -1), "_", tail(age_breaks, -1))

df_$AgeGroup <- cut(df_$Age, breaks = age_breaks, labels = age_labels, right = FALSE)

df_count <- df_ %>%
  group_by(MunCode, Date, AgeGroup) %>%
  summarize(Deaths = n(), .groups = 'drop')

write.csv(df_count, '../input_files/Age_COL_covid_data.csv')
