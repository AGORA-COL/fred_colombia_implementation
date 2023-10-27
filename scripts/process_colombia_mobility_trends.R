setwd('/zine/HPC02S1/ex-dveloza/AGORA/apps/fred_colombia_implementation/scripts')
##===============================================================#
## Get mobility data from Grandata
## Author: Diego Veloza
## Date: 2023/10/26
##===============================================================#
## Setup-------------
##===============================================================#
library(httr)
library(jsonlite)
library(lubridate)
library(dplyr)
library(tidyverse)
library(rjson)
library(sf)
library(raster)
library(rgdal)
library(maptools)
library(osmdata)
library(RColorBrewer)
library(osmplotr)
library(RCurl)
library(data.table)
options(digits = 22,scipen = 999)

library(stringi)
library(stringr)
library(scales)


metadata_file <- "/zine/HPC02S1/ex-dveloza/AGORA/apps/synthetic_populations/data/param_files/colombia_municipios_metadata.json"
json_data = rjson::fromJSON(
                                file = metadata_file,
                                simplify = F)


# Initialize an empty data frame to store the results
colombia_municipalities_info <- data.frame(
  department_name = character(),
  department_code = numeric(),
  mun_name = character(),
  divipola_code = numeric(),
  stringsAsFactors = FALSE
)

# Iterate through the departments
for (department_name in names(json_data$colombia)) {
  
  department_data <- json_data$colombia[[department_name]]
  
  # Iterate through the municipalities within each department
  for (mun_name in names(department_data)) {
    
    mun_data <- department_data[[mun_name]]
    
    department_code <- mun_data$department_code[[1]]
    divipola_code <- mun_data$divipola_code[[1]]
    
    # Create a new row with the extracted data
    new_row <- data.frame(
      department_name = department_name,
      department_code = department_code,
      mun_name = mun_name,
      divipola_code = divipola_code,
      stringsAsFactors = FALSE
    )
    
    # Append the new row to the result data frame
    colombia_municipalities_info <- rbind(colombia_municipalities_info, new_row)
  }
}

colombia_municipalities_info <- colombia_municipalities_info %>%
                                mutate(department_name = if_else(department_name == "ARCHIPIÉLAGO DE SAN ANDRÉS PROVIDENCIA Y",
                                                                            "SAN ANDRÉS Y PROVIDENCIA", department_name)) %>%
                                mutate(mun_name = stri_trans_general(mun_name, "Latin-ASCII"), department_name = stri_trans_general(department_name, "Latin-ASCII"))

##===============================================================#
## Read shp data-------------
##===============================================================#
## Maybe only download shp data if it's not downloaded yet')
mun_shp = rgdal::readOGR('/zine/HPC02S1/ex-dveloza/AGORA/apps/synthetic_populations/data/raw_data/geodata/Colombia_shp/Municipios.shp')
population_mun = read_csv('/zine/HPC02S1/ex-dveloza/AGORA/apps/synthetic_populations/data/processed_data/popdata/colombia_population_data_municp.csv')

##===============================================================#
## Read data-------------
##===============================================================#
date_seq = seq(from = as.Date('2020-03-02'), to = as.Date('2021-01-17'), by = 1)

file_facebook_data <- '/zine/HPC02S1/ex-dveloza/AGORA/apps/fred_colombia_implementation/input_files/facebook_open_data/movement-range-2022-05-22.txt'
file_gadm <- '/zine/HPC02S1/ex-dveloza/AGORA/apps/fred_colombia_implementation/input_files/gadm_data/gadm41_COL_2.shp'

facebook_data <- read.csv(file_facebook_data, sep = '\t') %>% 
                                filter(country == 'COL') %>% 
                                mutate(polygon_name = toupper(polygon_name))

facebook_data$num_dpto <- as.numeric(str_extract(facebook_data$polygon_id, "(?<=\\.)\\d+"))
facebook_data$num_mun  <- as.numeric(str_extract(facebook_data$polygon_id, "(?<=\\.)\\d+(?=_1)"))

facebook_data <- facebook_data %>%
  mutate(
    num_dpto = case_when(
      num_dpto >= 5 ~ num_dpto + 1,
      TRUE ~ num_dpto
    ),
    num_dpto = case_when(
      polygon_name == "SANTAFÉ DE BOGOTÁ" ~ 5,
      TRUE ~ num_dpto
    ),
    num_mun = case_when(
      polygon_name == "SANTAFÉ DE BOGOTÁ" ~ 1,
      TRUE ~ num_mun
    ),
    polygon_name = case_when(
      polygon_name == "SANTAFÉ DE BOGOTÁ" ~ "BOGOTÁ D.C.",
      TRUE ~ polygon_name
    )
  )


gadm_shp = rgdal::readOGR(file_gadm)

gadm_dpto_data <- gadm_shp@data %>% dplyr::select(GID_2, NAME_1, NAME_2) %>%
                mutate(num_dpto_gadm = as.numeric(str_extract(GID_2, "(?<=\\.)\\d+"))) %>%
                #mutate(num_mun = as.numeric(str_extract(GID_2, "(?<=\\.)\\d+(?=_2)"))) %>%
                rename("true_polygon_id" = "GID_2", "name_dpto" = "NAME_1", "name_mun" = "NAME_2") %>%
                dplyr::select(num_dpto_gadm, name_dpto) %>%
                distinct(num_dpto_gadm, name_dpto) %>%
                mutate(name_dpto = toupper(name_dpto), name_dpto = stri_trans_general(name_dpto, "Latin-ASCII"))


colombia_municipalities_info <- left_join(colombia_municipalities_info, gadm_dpto_data, 
                                            by = c("department_name" = "name_dpto")) 

replacement_values <- c(
  "PURISIMA"                    = "PURISIMA DE LA CONCEPCION",
  "SAN BERNARDINO DE SAHAGUN"   = "SAHAGUN",
  "SANTA CRUZ DE LORICA"        = "LORICA",
  "PUERTO INIRIDA"              = "INIRIDA",
  "SANTA MARTA (DIST. ESP.)"    = "SANTA MARTA",
  "DON MATIAS"                  = "DONMATIAS",
  "LA UNION DE SUCRE"           = "LA UNION",
  "SAN VICENTE"                 = "SAN VICENTE FERRER",
  "SAN LUIS DE CUBARRAL"        = "CUBARRAL",
  "VISTA HERMOSA"               = "VISTAHERMOSA",
  "SAN JOSE DE CUCUTA"          = "CUCUTA",
  "SAN MIGUEL DE MOCOA"         = "MOCOA",
  "SANTIAGO DE CALI"            = "CALI",
  "SAN JUAN DE PASTO"           = "PASTO",
  "SINCE"                       = "SINCELEJO",
  "TOLU"                        = "SANTIAGO DE TOLU",
  "ARMERO"                      = "ARMERO GUAYABAL",
  "TUMACO"                      = "SAN ANDRES DE TUMACO"
)

facebook_data <- facebook_data %>%
                mutate(polygon_name = stri_trans_general(polygon_name, "Latin-ASCII")) %>%
                mutate(polygon_name = case_when(
                    num_dpto == 2 & polygon_name == "BOLIVAR" ~ "CIUDAD BOLIVAR",
                    TRUE ~ recode(polygon_name, !!!replacement_values)
                ))


facebook_data_divipola_code <- left_join(facebook_data, colombia_municipalities_info, by = c("polygon_name" = "mun_name", "num_dpto" = "num_dpto_gadm"))
# write_csv(facebook_data_divipola_code, '../input_files/COL_mobility_trends_facebook.csv')

##===============================================================#
## Combine shape files-------------
##===============================================================#

population_mun <- population_mun %>% filter(Code != "00", Code != Zone, Year == '2018', Gender == 'Total')
dept_list <- unique(population_mun$Code)
for(dept_code in dept_list){
    cat(crayon::red(sprintf("\nCurrent department : %s\n", dept_code)))

    dept_population <- population_mun %>% filter(Code == dept_code)
    total_population <- sum(dept_population$Pop)
    mun_list <- unique(dept_population$Zone)

    dept_all_days_list <- list()

    # Loop through each unique date in the facebook_data_divipola_code dataframe
    unique_dates <- unique(facebook_data_divipola_code$ds)
    for(date in unique_dates){
        # Get the mobility data for the current department and date
        dept_mobility <- facebook_data_divipola_code %>% 
                    filter(department_code == as.numeric(dept_code) & ds == date)
        
        # Identify the municipalities with and without mobility data
        mun_with_data <- unique(dept_mobility$divipola_code)
        mun_without_data <- setdiff(as.numeric(mun_list), as.numeric(mun_with_data))
        
        # If there are municipalities without data, estimate their mobility trend
        if(length(mun_without_data) > 0){
            # Sum up the trends of municipalities with data, weighted by their population
            weighted_sum <- sum(dept_mobility$all_day_bing_tiles_visited_relative_change * 
                                    dept_population %>% filter(as.numeric(Zone) %in% mun_with_data) %>% pull(Pop))
            
            # Get the total population of municipalities with data
            total_pop_with_data <- sum(dept_population %>% filter(as.numeric(Zone) %in% mun_with_data) %>% pull(Pop))
            
            # Calculate the average trend for the department
            avg_trend <- weighted_sum / total_pop_with_data
            
            # Create a data frame for the municipalities without data
            mun_without_data_df <- dept_population %>% 
                filter(as.numeric(Zone) %in% mun_without_data) %>%
                mutate(
                department_code = as.numeric(dept_code),
                divipola_code = as.numeric(Zone),
                ds = date,
                all_day_bing_tiles_visited_relative_change = avg_trend * (Pop / total_pop_with_data)
                ) %>%
                dplyr::select(ds, divipola_code, department_code, all_day_bing_tiles_visited_relative_change)
        }
        dept_all_days_list[[length(dept_all_days_list) + 1]] = bind_rows(dept_mobility %>% 
                                                                        dplyr::select(ds, 
                                                                                      divipola_code, 
                                                                                      department_code,
                                                                                      all_day_bing_tiles_visited_relative_change), mun_without_data_df)
    }

    dept_mobility <- bind_rows(dept_all_days_list)
    write_csv(dept_mobility, sprintf('/zine/HPC02S1/ex-dveloza/AGORA/apps/fred_colombia_implementation/input_files/%s_facebook_mobility_trends.csv', dept_code))
}


##===============================================================#
## Plot trends by locality-------------
##===============================================================#
my_colors = brewer.pal(6,"RdYlBu")
my_colors = colorRampPalette(my_colors)(13)
mov_brk <- cut(as.numeric(dept_mobility$all_day_bing_tiles_visited_relative_change), breaks = c(-1,-0.5,-0.2,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,2.0), include.lowest = T, right = F)
my_colors = my_colors[as.numeric(mov_brk)]

dept_mobility$ds <- as.Date(dept_mobility$ds)

# Normalizing the values
max_values <- dept_mobility %>%
  group_by(divipola_code) %>%
  summarize(max_value = max(all_day_bing_tiles_visited_relative_change, na.rm = TRUE))

dept_mobility <- dept_mobility %>%
  left_join(max_values, by = "divipola_code") %>%
  mutate(normalized_values = all_day_bing_tiles_visited_relative_change / max_value)

# Using ggplot2 for the visualization
p <- ggplot(dept_mobility, aes(x = ds, y = normalized_values, color = as.factor(divipola_code))) +
  geom_line() +
  scale_color_manual(values = my_colors) +
  labs(title = "Normalized Mobility Trend by Region", 
       x = "Date", 
       y = "Normalized Mobility Trend", 
       color = "Region Code") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(filename = "mobility_plot.png", plot = p, width = 10, height = 6, dpi = 300)


##===============================================================#
## Plot mobility in shape files-------------
##===============================================================#

# For simplicity, let's take the latest mobility data for each region
latest_mobility <- dept_mobility %>%
  group_by(divipola_code) %>%
  arrange(desc(ds)) %>%
  slice(1)

mun_shp_sf <- st_as_sf(mun_shp)
mun_shp_sf$ID_ESPACIA <- as.numeric(mun_shp_sf$ID_ESPACIA)

# Merge shapefile with the mobility data
merged_data <- left_join(mun_shp_sf[as.numeric(mun_shp@data$ID_ESPACIA) %in% latest_mobility$divipola_code, ], latest_mobility, by = c("ID_ESPACIA" = "divipola_code"))

# Plotting the map
q = ggplot(data = merged_data) +
  geom_sf(aes(fill = normalized_values), color = "white") + # Assuming normalized_values is what you want to plot
  scale_fill_viridis_c() + # You can use any other color scale
  labs(title = "Latest Normalized Mobility Trends by Region",
       fill = "Mobility") +
  theme_minimal()

ggsave(filename = "mobility_plot.png", plot = q, width = 10, height = 6, dpi = 300)

##===============================================================#
## Just for fun, make a video-------------
##===============================================================#
# Directory to save temporary images
img.dir <- "tmpimg"

# Clean up old images
if(file.exists(img.dir)){
    system(paste('rm -rf ', img.dir,sep = ''))
}
system(paste('mkdir ', img.dir,sep = ''))

# Assuming that 'date_seq' is the sorted unique set of dates in your dataset
date_seq <- sort(unique(merged_data_all_dates$ds))

# Calculate the overall min and max values for normalized_values
overall_min <- min(merged_data_all_dates$normalized_values, na.rm = TRUE)
overall_max <- max(merged_data_all_dates$normalized_values, na.rm = TRUE)

for(dd in 1:length(date_seq)){
    cat(crayon::red(sprintf("\r\b%s", date_seq[dd])))
    filename <- file.path(img.dir, sprintf("mobility_trends_%06d.png", dd - 1))
    png(filename, width = 1024, height = 800)
    
    tmp_df <- merged_data_all_dates %>% filter(ds == date_seq[dd])
    
    # Adapted ggplot code with fixed color scale
    q <- ggplot(data = tmp_df) +
        geom_sf(aes(fill = normalized_values), color = "white") + 
        scale_fill_viridis_c(limits = c(overall_min, overall_max)) +  # Fixed scale limits
        labs(title = sprintf("Mobility Trends: %s", date_seq[dd]),
             fill = "Mobility") +
        theme_minimal()
    
    print(q)
    dev.off()
}


# Generate the video
video_name <- "Mobility_trends_animation"

output_movie <- sprintf('%s.mp4', video_name)

ffmpeg_str <- sprintf(
    "ffmpeg -r 3 -i %s/%s -c:v libx264 -r 30 -pix_fmt yuv420p %s -y ", 
    img.dir,  'mobility_trends_%06d.png',   output_movie )

system(ffmpeg_str)
