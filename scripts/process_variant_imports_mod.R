setwd('/zine/HPC02S1/ex-dveloza/AGORA/apps/fred_colombia_implementation/scripts')
##===============================#
## Plot model fit for Bogotá
## Author: Guido España
## 2020
##===============================#
## Setup-------------
##===============================#
set.seed(123456)
library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(lubridate)
library(rjson)
library(mgcv)
library(MASS)

date_format = '20210818'

args = (commandArgs(TRUE))

if(length(args) >= 1){
    date_format = args[1]
}

##===============================#
##Functions-------------
##===============================#
parse_inner_json <- function(x_in, xx){
    cluster_names = unlist(x_in[['cluster_names']])
    country_data = x_in[['distributions']][[xx]]

    country_df = do.call(rbind,lapply(country_data[['distribution']], FUN = function(x){data.frame(variants = names(unlist(x[['cluster_counts']])), counts = as.numeric(unlist(x[['cluster_counts']])), week = x[['week']],total_counts = x[['total_sequences']], stringsAsFactors = F)}))
    country_df$country = country_data$country

    return(country_df)
}

##===============================#
## Read data-------------
##===============================#
data_url = "https://raw.githubusercontent.com/hodcroftlab/covariants/master/web/data/perCountryData.json"
download.file(data_url, '../data/coronavirus_variants_per_country.json')

covariants_list = rjson::fromJSON(file = '../data/coronavirus_variants_per_country.json',
                                  simplify = F)

covariants_regions = covariants_list[['regions']]


region_1_df = do.call(rbind, lapply(X = 1:length(covariants_regions[[1]][['distributions']]), FUN = function(x){parse_inner_json(covariants_regions[[1]], x)}))
region_1_df$dominance = region_1_df$counts / region_1_df$total_counts


uk_df = region_1_df %>% filter(country == 'United Kingdom')
total_variants = unique(uk_df$variants)
uk_df$week = as.Date(uk_df$week)
vv = "20I (Alpha, V1)"

tmp_df = filter(uk_df, variants == vv)
##plot(tmp_df$week, tmp_df$dominance, type = 'o')
max_date = as.Date(date_format, format = "%Y%m%d")

## Get a minimum and maximum week
## For each country, fit a gam and extrapolate dominance to the timeseries
total_variants_df = data.frame(stringsAsFactors = F)
for(cc in unique(region_1_df$country)){
    tmp_df = filter(region_1_df, country == cc)
    tmp_df$day = as.numeric(as.Date(tmp_df$week) - min(as.Date(tmp_df$week) ))    
    country_variants = unique(tmp_df$variants)
    for(vv in country_variants){
        tmp_v_df = filter(tmp_df, variants == vv)
        if(nrow(tmp_v_df) > 10){
            var_df = data.frame(Date = seq(from = min(as.Date(tmp_df$week)), to = max_date, by = '1 day'), variant = vv, country = cc, stringsAsFactors = F)
            var_df$day = as.numeric(var_df$Date - min(var_df$Date))
            mod_gam0 = gam(dominance ~ s(day), data=tmp_v_df, family=gaussian(link='identity'))        
            var_df$dominance = predict(mod_gam0, newdata = var_df)
            total_variants_df = bind_rows(total_variants_df, var_df)
        }
    }
}
total_variants_df$dominance[total_variants_df$dominance < 0] = 0

##==================================#
## Read import data--------
##==================================#
country_names_list = list(
    "ALEMANIA" = "Germany",
    "ARGENTINA" = "Argentina",
    "ARUBA" = "Aruba",
    "BRASIL" = "Brazil",
    "CANADÁ" = "Canada",
    "CHILE" = "Chile",
    "COSTA RICA" = "Costa Rica",
    "ECUADOR" = "Ecuador",
    "EL SALVADOR" = "El Salvador",
    "EGIPTO" = "Egypt",
    "ESPAÑA" = "Spain",
    "ESTADOS UNIDOS DE AMÉRICA" = "USA",
    "FRANCIA" = "France",
    "GUATEMALA" = "Guatemala",
    "ITALIA" = "Italy",
    "JAMAICA" = "Jamaica",
    "MÉXICO" = "Mexico",
    "PANAMA" = "Panama",
    "PANAMÁ" = "Panama",
    "PERU" = "Peru",
    "PERÚ" = "Peru",
    "PUERTO RICO" = "USA",
    "REINO UNIDO DE GRAN BRETAÑA E IRLANDA DEL NORTE" = "United Kingdom",
    "REPÚBLICA DOMINICANA" = "Dominican Republic",
    "TURQUÍA" = "Turkey",
    "VENEZUELA" = "Venezuela")

country_names_df = data.frame(Country = unlist(country_names_list), CountryName = names(country_names_list))

covid_file = sprintf('../epidata/Covid_Colombia_cases_%s.csv',date_format)
reports_covid_df = read_csv(covid_file) %>%
    rename(DateDiagnosis = "Fecha de diagnóstico",
           CountryName = "Nombre del país",
           DateReport = "fecha reporte web",
           SymptomOnset = "Fecha de inicio de síntomas",
           MunCode = "Código DIVIPOLA municipio",
           CaseType = "Tipo de contagio",
           DateDeath = "Fecha de muerte") %>%
    dplyr::select(MunCode, CaseType, DateDeath,  DateReport,DateDiagnosis, SymptomOnset, CountryName) %>%
    mutate(DateDeath = str_replace(DateDeath, "-   -", "")) %>%
    filter(DateDiagnosis != as.Date('1899-12-31')) %>%
    mutate(DateDiagnosis = str_replace_all(DateDiagnosis, "0+:.*", ""),
           DateReport = str_replace_all(DateReport, "0+:.*", ""),
           DateDeath = str_replace_all(DateDeath, "0+:.*", "")) %>%
    mutate(DateDiagnosis = parse_date_time(DateDiagnosis, orders = c("dmy", "dmy HMS")))%>%
    mutate(SymptomOnset = parse_date_time(SymptomOnset, orders = c("dmy", "dmy HMS")))%>%
    mutate(DateDeath = parse_date_time(DateDeath, orders = c("dmy", "dmy HMS")))%>%
    mutate(DateDiagnosis = lubridate::ymd(DateDiagnosis))%>%
    mutate(DateReport = parse_date_time(DateReport, orders = c("dmy", "dmy HMS"))) %>%
    mutate(DateReport = lubridate::ymd(DateReport)) %>%
    mutate(SymptomOnset = lubridate::ymd(SymptomOnset))%>%
    mutate(DateDeath = lubridate::ymd(DateDeath)) 


##==================================#
## Combine with import data--------
##==================================#
alternative_imports_df = read_csv('../experiments_colombia/input_files/11001_imports_alternative.csv') %>%
    group_by(Date) %>%
    summarize(Imports = mean(Imports)) %>%
    ungroup() %>%
    mutate(Replicate = 1)

imports_df = reports_covid_df %>% filter(CaseType == "Importado") %>%
    filter(MunCode == 11001) %>%
    dplyr::select(DateReport, CountryName) %>%    
    left_join(country_names_df, by = 'CountryName') %>%
    drop_na() %>%
    dplyr::select(DateReport, Country)

total_imports = imports_df %>% group_by(DateReport) %>% summarize(Imports = n()) %>% ungroup() %>%
    rename(Date = DateReport) %>%
    mutate(Date = Date - 3,
           Imports = Imports * 3.3) %>%
    bind_rows(alternative_imports_df) %>%
    group_by(Date) %>%
    summarize(Imports = max(Imports)) %>%
    ungroup()

total_imports$Day = as.numeric(total_imports$Date - min(total_imports$Date))

mod_gam_imp = gam(Imports ~ s(Day), data=total_imports, family=gaussian(link='identity'))
total_imports_pr = data.frame(Date = seq(from = min(total_imports$Date), to = max_date, by = '1 day'),
                              stringsAsFactors = F)

total_imports_pr$Day = as.numeric(total_imports_pr$Date - min(total_imports_pr$Date))
total_imports_pr$Imports = predict(mod_gam_imp, newdata = total_imports_pr)

total_imports = total_imports_pr

country_imports = imports_df %>% group_by(Country) %>% summarize(Imports = n()) %>%
    filter(Imports > 10) %>%
    mutate(ImportsProp = Imports / sum(Imports))

total_variants = unique(region_1_df$variants)

variants_imports = total_variants_df %>% left_join(country_imports, by = c('country' = 'Country')) %>%
    mutate(PropVariantImports = dominance * ImportsProp) %>%
    group_by(variant, Date) %>%
    summarize(PropVariantImports = sum(PropVariantImports, na.rm = T)) %>%
    ungroup() %>%
    left_join(total_imports, by = c('Date' = 'Date')) %>%
    replace_na(list(Imports = 0)) %>%
    mutate(TotalVariantImports = Imports * PropVariantImports)

write_csv(variants_imports, '../experiments_colombia/input_files/COL_variant_imports.csv')

##==================================#
## Plot imports by variant--------
##==================================#
jpeg('../figures/process_figure_variants_imports_colombia.jpeg', width=8.5,height=6, units="in", res = 300)
date_ticks =  seq(from=min(variants_imports$Date),to=max(variants_imports$Date), by = '1 month')
date_labels = sprintf("%s-%s", as.character(month(date_ticks,label = T)), as.character(day(date_ticks)))

voc_list = c("20I (Alpha, V1)", "20J (Gamma, V3)", "21A (Delta)", "21B (Kappa)",
             "20H (Beta, V2)", "21C (Epsilon)", "21D (Eta)", "21F (Iota)")

col_variants = brewer.pal(n = length(voc_list), name = "Set3")
layout(matrix(1:4, ncol = 2))
pplt = par("plt")
adjx = (0.1 - pplt[1]) / (pplt[2] - pplt[1])


barplot(country_imports$ImportsProp, names.arg = country_imports$Country, las = 2)
mtext("Proportion of imports in Bogota", side = 2, line = 2.5, cex = 0.6)
mtext(LETTERS[1], side = 3, line = 0.5, outer = F, cex = 1.0, adj = adjx)

plot(total_imports$Date, total_imports$Imports, lwd = 2, type = 'l',
     col = 'black', ylab = "Imports in Bogota", xlab = "Date", xlim = c(min(variants_imports$Date),max(variants_imports$Date)), xaxt = 'n')
axis.Date(1,variants_imports$Date, at = date_ticks, labels = date_labels, cex.axis = 0.8, las = 2)
mtext(LETTERS[2], side = 3, line = 0.5, outer = F, cex = 1.0, adj = adjx)

for(vv in 1:length(voc_list)){
    tmp_var = filter(variants_imports, variant == voc_list[vv])
    if(vv == 1){
        plot(tmp_var$Date, tmp_var$TotalVariantImports*1, type = 'l', col = col_variants[vv], lwd = 1.5,
             xaxt = 'n', ylab = "Total Imports from Variant", xlab = "Date")
        axis.Date(1,tmp_var$Date, at = date_ticks, labels = date_labels, cex.axis = 0.8, las = 2)
    }else{
        lines(tmp_var$Date, tmp_var$TotalVariantImports*1, type = 'l', col = col_variants[vv], lwd = 1.5)
    }    
}
legend("topleft",legend = voc_list,col = col_variants, lwd = 2, ncol = 2, cex = 0.5)
mtext(LETTERS[3], side = 3, line = 0.5, outer = F, cex = 1.0, adj = adjx)

for(vv in 1:length(voc_list)){
    tmp_var = filter(variants_imports, variant == voc_list[vv])
    if(vv == 1){
        plot(tmp_var$Date, tmp_var$PropVariantImports, type = 'l', col = col_variants[vv], lwd = 1.5,
             xaxt = 'n', ylab = "Proportion of Imports from Variant", xlab = "Date")
        axis.Date(1,tmp_var$Date, at = date_ticks, labels = date_labels, cex.axis = 0.8, las = 2)
    }else{
        lines(tmp_var$Date, tmp_var$PropVariantImports, type = 'l', col = col_variants[vv], lwd = 1.5)
    }    
}
mtext(LETTERS[4], side = 3, line = 0.5, outer = F, cex = 1.0, adj = adjx)

dev.off()
