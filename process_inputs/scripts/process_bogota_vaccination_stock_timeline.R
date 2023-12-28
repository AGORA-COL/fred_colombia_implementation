##===============================#
## Process vaccination rollout
## Author: Guido España
## 2021
##===============================#
## Setup-------------
##===============================#
library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(lubridate)

##===============================#
## Inputs-------------
##===============================#
start_date = '2020-01-01'
forecast_date = '2021-12-31'
prop_bogota = 7.5/50

##===============================#
## Read data-------------
##===============================#
## 1. Download vaccine stock for each vaccine
download.file("https://docs.google.com/spreadsheets/d/1N4V-qqp6q8exuYAZgZ79GT9GazsSCHUK/export?format=xlsx",
              "./input_files/COL_Vaccines_Administered.xlsx")


download.file("https://docs.google.com/spreadsheets/d/1z2KYfMvDMLHb3f1xQMDHM5Q9ll_vIwe764XBBQF7P2E/export?format=xlsx",
              "./input_files/COL_Vaccines_Stock.xlsx")

## 2. Assign a number to each type of vaccine
vaccine_params = read_csv('../input_files/COVID19_Vaccine_Parameters_Colombia.csv')%>%
    mutate(Vaccine = tolower(Vaccine)) %>%
    dplyr::select(Vaccine, ID)

## TODO: How to deal with COVAX?
vaccine_stock_df = readxl::read_xlsx('./input_files/COL_Vaccines_Stock.xlsx', skip = 9, sheet = 'Originales') %>%
    rename(StockType = 'Operación', TotalVaccines = 'Bogotá', Vaccine = 'TipoVacuna') %>%
    mutate(Date = as.Date(Fecha)) %>%
    dplyr::select(Date, Vaccine, StockType,TotalVaccines) %>%
    mutate(Vaccine = ifelse(Vaccine == "Pfizer-Covax", "Pfizer", Vaccine)) %>%
    mutate(Vaccine = tolower(Vaccine)) %>%
    drop_na() %>%
    dplyr::filter(StockType == 1, TotalVaccines > 0) %>%
    left_join(vaccine_params, by = "Vaccine")     



## 3. Process daily vaccination numbers
vaccine_daily_df = readxl::read_xlsx("./input_files/COL_Vaccines_Administered.xlsx",sheet = "Dosis aplicadas") %>%
    mutate(Date = as.Date(as.numeric(FECHAS), origin="1899-12-30")) %>%
    filter(!is.na(Date)) %>%
    dplyr::select(-FECHAS) %>%
    gather(key = Department, value = VaccinesApplied, -Date) %>%
    filter(Department == "Bogotá")

##TODO: Project vaccination into the future

write_csv(vaccine_daily_df, './input_files/11001_vaccine_capacity_timeseries.csv')
write_csv(vaccine_stock_df, './input_files/11001_vaccine_stock_timeseries.csv')

##===============================#
## plot figure-------------
##===============================#
col_palette = brewer.pal(6, 'Set1')

jpeg('./fred_bog_vaccination_stock.jpeg', width = 9, height = 4, units = 'in', res = 300)
layout(mat = matrix(1:2,nrow = 1))
barplot(height = vaccine_stock_df$TotalVaccines, col = col_palette[vaccine_stock_df$ID + 1],
        lwd = 2, ylim = c(0, max(vaccine_stock_df$TotalVaccines)*1.1), xlab = '', ylab = 'Vacunas', main = 'Llegada de vacunas a Bogotá',names.arg = vaccine_stock_df$Date, las = 2, cex.axis = 0.5, cex.names = 0.5)

legend('topleft',legend = vaccine_params$Vaccine, fill = col_palette[vaccine_params$ID + 1],col = col_palette[vaccine_params$ID + 1], ncol = 2, cex = 0.5)

plot(vaccine_daily_df$Date, vaccine_daily_df$VaccinesApplied, col = 'black', lwd = 1, type = 'b',ylab = 'Número de dosis aplicadas', xlab = '', xaxt = 'n', main = 'Dosis aplicadas diariamente', xaxs = 'i', yaxs = 'i')

date_ticks =  seq(from=min(vaccine_daily_df$Date),to=max(vaccine_daily_df$Date), by = '1 month')
date_labels = (date_ticks)
axis.Date(1,vaccine_daily_df$Date, at = date_ticks, labels = date_labels, cex.axis = 0.8, las = 1)

grid(nx = length(date_ticks)+1, ny = NULL)
dev.off()


##===============================#
## process data-------------
##===============================#
## output_file = '../misc_scripts/vaccine_stock_bogota_file.txt'

## vaccine_stock_lines = sprintf("%.0f %.0f %.0f %.0f",vacc_df$day, vacc_df$day,vacc_df$bog_doses,0)

## fileConn<-file(output_file)
## writeLines(vaccine_stock_lines, fileConn)
## close(fileConn)


## output_file = '../misc_scripts/vaccine_capacity_bogota_file.txt'

## ## vaccine_capacity_lines = c(sprintf("%.0f %.0f",bog_df$day, bog_df$doses),
## ##                         sprintf("%.0f %.0f",bog_df$day[nrow(bog_df)]+1, 0))

## vaccine_capacity_lines = c(sprintf("%.0f %.0f",daily_bog_df$day, daily_bog_df$doses))

## fileConn<-file(output_file)
## writeLines(vaccine_capacity_lines, fileConn)
## close(fileConn)
