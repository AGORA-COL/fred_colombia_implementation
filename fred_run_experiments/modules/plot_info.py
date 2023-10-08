import pandas as pd
import matplotlib.pyplot as plt
import warnings
import matplotlib.colors as mcolors
import matplotlib.transforms as mtransforms
import datetime as dt
import subprocess
import matplotlib as mpl
import utils as ut
import os
import unidecode
mpl.rcParams['axes.linewidth'] = 1.75

colores = []
for name, color in mcolors.TABLEAU_COLORS.items():
    colores.append(name[4:])
warnings.filterwarnings('ignore')

######################
### Simul settings ###
######################
forecast_day = dt.datetime(2022,12,1)
initial_day = dt.datetime(2020,6,1)
projection_label = 'omicron_lineages'

folder_name = 'FRED_11001_projections_asymp_1.00_fm_0.73_ksus_10.00_var_1_vax_070_mov_'+projection_label+'_out'

## Succesfull simulation parameters
df_params = pd.read_csv('../../output/SHORT_FORECAST/'+folder_name+'/FRED_parameters_out.csv')

######################
#### Update info #####
######################
## Deaths data
ut.FRED_muertes().to_csv('../../data/COL_covid_death_data.csv')

## Rt
if os.path.exists('../../data/Numero-de-Reproduccion-Efectivo-R(t) todos los casos.csv'):
    os.remove('../../data/Numero-de-Reproduccion-Efectivo-R(t) todos los casos.csv')
    
file_url = '"https://saludata.saludcapital.gov.co/osb/wp-content/uploads/medios/Numero-de-Reproduccion-Efectivo-R(t)%20todos%20los%20casos.csv"'
subprocess.run("wget "+file_url+" -P ../../data/", shell=True, check=True)

path = '../../data/Numero-de-Reproduccion-Efectivo-R(t) todos los casos.csv'
with open(path, encoding='latin1') as f:
    contents = f.readlines()

column_names = unidecode.unidecode(contents[0].upper())
for r in ((' ;', ';'), (' ;', ';'),(' ','_')):
    column_names = column_names.replace(*r)

os.remove('../../data/Numero-de-Reproduccion-Efectivo-R(t) todos los casos.csv')
file = open('../../data/Numero-de-Reproduccion-Efectivo-R(t) todos los casos.csv', 'a')
for i in range(5, len(contents)):
    file.writelines(unidecode.unidecode(contents[i].replace(',','.')))
file.close()

###############################
##### Process simuls info #####
###############################
## Process simuls info
subprocess.run("Rscript ~/FRED_Implementation/fit_scripts/fit_post_process.R "+projection_label, shell=True, check=True)

#####################
##### Load info #####
#####################
## Simul post process
file = f'../../../output/fred_output_model_fit_{projection_label}.csv'
df_fit = pd.read_csv(file)

## Dominance data
#file = '../../input_files/Bogota_Covid_Variants_Dominance_IC.csv'
file = '../../input_files/var_prop.csv'
df_dom_data = pd.read_csv(file)

## Rt
file_data = '../../data/Numero-de-Reproduccion-Efectivo-R(t) todos los casos.csv'
data_Rt = pd.read_csv(file_data, delimiter = ';', encoding='latin')

## Muertes diarias
file_data = '../../data/COL_covid_death_data.csv'
df_data = pd.read_csv(file_data, encoding = 'latin')

## Simulations
file = '../../output/SHORT_FORECAST/'+folder_name+'/fred_output.csv'
df = pd.read_csv(file)

for scenario in df_params['intervention_id'].unique():
    ########################
    ##### Data setting #####
    ########################
    title = scenario
    df_scenario = df_params[df_params['intervention_id'] == scenario]

    ########################
    ##### Plot setting #####
    ########################
    plt.rcParams.update({'font.size': 15})
    fig1, ax = plt.subplots(2,2, tight_layout=True)
    fig1.set_size_inches(w=20, h=12)
    fig1.suptitle(title)
    P = []
    legend_label = []

    ##################
    ##### Deaths #####
    ##################
    death_data = df_data[df_data['MunCode'] == 11001]['Deaths'].to_numpy()
    death_days = (pd.to_datetime(df_data[df_data['MunCode'] == 11001]['Date']))
    ax[0,0].scatter(death_days, death_data, s = 10, c = 'black', label = 'Muertes reportadas', zorder = 2)
    ax[0,0].legend()

    start_day = dt.datetime(2020,1,1)
#     for simul_id in df_scenario['job_id'].to_numpy():
#         mask = df['job_id'] == simul_id
#         max_df = df[mask]

#         days = max_df['Day']
#         date_days = pd.np.array([start_day + dt.timedelta(days = i) for i in days])
#         deaths_0 = max_df['CF_mean'].to_numpy()
#         deaths_1 = max_df['CF_1_mean'].to_numpy()
#         deaths_2 = max_df['CF_2_mean'].to_numpy()
#         deaths_3 = max_df['CF_3_mean'].to_numpy()
#         deaths_4 = max_df['CF_4_mean'].to_numpy()
#         deaths_5 = max_df['CF_5_mean'].to_numpy()

#         ### Simulation
#         ax[0,0].plot(date_days, deaths_0, alpha = 1, zorder = 1, color = colores[0], linewidth = 2)
#         ax[0,0].plot(date_days, deaths_1, alpha = 1, zorder = 1, color = colores[1], linewidth = 2)
#         ax[0,0].plot(date_days, deaths_2, alpha = 1, zorder = 1, color = colores[2], linewidth = 2)
#         ax[0,0].plot(date_days, deaths_3, alpha = 1, zorder = 1, color = colores[3], linewidth = 2)
#         ax[0,0].plot(date_days, deaths_4, alpha = 1, zorder = 1, color = colores[4], linewidth = 2)
#         ax[0,0].plot(date_days, deaths_5, alpha = 1, zorder = 1, color = colores[5], linewidth = 2)
    mask = df_fit['intervention_id'] == scenario


    CF_fit_low = df_fit[mask]['CF_low'].to_numpy()
    CF_fit_high = df_fit[mask]['CF_high'].to_numpy()
    CF_fit_median = df_fit[mask]['CF_median'].to_numpy()
    
    CF_0_fit_low = df_fit[mask]['CF_0_low'].to_numpy()
    CF_0_fit_high = df_fit[mask]['CF_0_high'].to_numpy()
    CF_0_fit_median = df_fit[mask]['CF_0_median'].to_numpy()

    CF_1_fit_low = df_fit[mask]['CF_1_low'].to_numpy()
    CF_1_fit_high = df_fit[mask]['CF_1_high'].to_numpy()
    CF_1_fit_median = df_fit[mask]['CF_1_median'].to_numpy()

    CF_2_fit_low = df_fit[mask]['CF_2_low'].to_numpy()
    CF_2_fit_high = df_fit[mask]['CF_2_high'].to_numpy()
    CF_2_fit_median = df_fit[mask]['CF_2_median'].to_numpy()

    CF_3_fit_low = df_fit[mask]['CF_3_low'].to_numpy()
    CF_3_fit_high = df_fit[mask]['CF_3_high'].to_numpy()
    CF_3_fit_median = df_fit[mask]['CF_3_median'].to_numpy()

    CF_4_fit_low = df_fit[mask]['CF_4_low'].to_numpy()
    CF_4_fit_high = df_fit[mask]['CF_4_high'].to_numpy()
    CF_4_fit_median = df_fit[mask]['CF_4_median'].to_numpy()

    CF_5_fit_low = df_fit[mask]['CF_5_low'].to_numpy()
    CF_5_fit_high = df_fit[mask]['CF_5_high'].to_numpy()
    CF_5_fit_median = df_fit[mask]['CF_5_median'].to_numpy()

    date_CF_fit = pd.to_datetime(df_fit[mask]['Date'].to_numpy(), format='%Y-%m-%d', errors='coerce')
    p1 = ax[0,0].plot(date_CF_fit, CF_fit_median,  alpha = 1, zorder = 2, color = 'grey', linewidth = 2)
    ax[0,0].fill_between(date_CF_fit, CF_fit_low, CF_fit_high, where=CF_fit_high >= CF_fit_low, facecolor='grey', interpolate=True, alpha = 0.1)

    p2 = ax[0,0].plot(date_CF_fit, CF_0_fit_median,  alpha = 1, zorder = 2, color = colores[0], linewidth = 2)
    ax[0,0].fill_between(date_CF_fit, CF_0_fit_low, CF_0_fit_high, where=CF_0_fit_high >= CF_0_fit_low, facecolor=colores[0], interpolate=True, alpha = 0.1)

    p3 = ax[0,0].plot(date_CF_fit, CF_1_fit_median,  alpha = 1, zorder = 2, color = colores[1], linewidth = 2)
    ax[0,0].fill_between(date_CF_fit, CF_1_fit_low, CF_1_fit_high, where=CF_1_fit_high >= CF_1_fit_low, facecolor=colores[1], interpolate=True, alpha = 0.1)

    p4 = ax[0,0].plot(date_CF_fit, CF_2_fit_median,  alpha = 1, zorder = 2, color = colores[2], linewidth = 2)
    ax[0,0].fill_between(date_CF_fit, CF_2_fit_low, CF_2_fit_high, where=CF_2_fit_high >= CF_2_fit_low, facecolor=colores[2], interpolate=True, alpha = 0.1)

    p5 = ax[0,0].plot(date_CF_fit, CF_3_fit_median,  alpha = 1, zorder = 2, color = colores[3], linewidth = 2)
    ax[0,0].fill_between(date_CF_fit, CF_3_fit_low, CF_3_fit_high, where=CF_3_fit_high >= CF_3_fit_low, facecolor=colores[3], interpolate=True, alpha = 0.1)

    p6 = ax[0,0].plot(date_CF_fit, CF_4_fit_median,  alpha = 1, zorder = 2, color = colores[4], linewidth = 2)
    ax[0,0].fill_between(date_CF_fit, CF_4_fit_low, CF_4_fit_high, where=CF_4_fit_high >= CF_4_fit_low, facecolor=colores[4], interpolate=True, alpha = 0.1)

    p7 = ax[0,0].plot(date_CF_fit, CF_5_fit_median,  alpha = 1, zorder = 2, color = colores[5], linewidth = 2)
    ax[0,0].fill_between(date_CF_fit, CF_5_fit_low, CF_5_fit_high, where=CF_5_fit_high >= CF_5_fit_low, facecolor=colores[5], interpolate=True, alpha = 0.1)

    ## Subplot details
    ax[0,0].set_xticks(pd.date_range(start_day, periods=70, freq='15d'))
    plt.setp(ax[0,0].get_xticklabels(), rotation=90, fontsize=10)
    ax[0,0].set_ylim(0,400)
    ax[0,0].set_xlim(initial_day, forecast_day)
    ax[0,0].set_ylabel('Muertes')
    ax[0,0].grid()

    #####################
    ##### Dominance #####
    #####################
    mask = df_dom_data['var'] == 'Alpha'
    date_alpha = pd.to_datetime(df_dom_data[mask]['date'].to_numpy(), format='%Y-%m-%d', errors='coerce')
    data_alpha = df_dom_data[mask]['PointEst'].to_numpy()

    mask = df_dom_data['var'] == 'Gamma'
    date_gamma = pd.to_datetime(df_dom_data[mask]['date'].to_numpy(), format='%Y-%m-%d', errors='coerce')
    data_gamma = df_dom_data[mask]['PointEst'].to_numpy()

    mask = df_dom_data['var'] == 'Mu'
    date_mu = pd.to_datetime(df_dom_data[mask]['date'].to_numpy(), format='%Y-%m-%d', errors='coerce')
    data_mu = df_dom_data[mask]['PointEst'].to_numpy()

    mask = df_dom_data['var'] == 'Delta'
    date_delta = pd.to_datetime(df_dom_data[mask]['date'].to_numpy(), format='%Y-%m-%d', errors='coerce')
    data_delta = df_dom_data[mask]['PointEst'].to_numpy()

    mask = df_dom_data['var'] == 'Omicron'
    date_omicron = pd.to_datetime(df_dom_data[mask]['date'].to_numpy(), format='%Y-%m-%d', errors='coerce')
    data_omicron = df_dom_data[mask]['PointEst'].to_numpy()

    mask = df_dom_data['var'] == 'BA.4 Omicron'
    date_BA4_omicron = pd.to_datetime(df_dom_data[mask]['date'].to_numpy(), format='%Y-%m-%d', errors='coerce')
    data_BA4_omicron = df_dom_data[mask]['PointEst'].to_numpy()

    ## Plot dominance data
    ax[0,1].scatter(date_alpha, data_alpha, color =  colores[1])
    ax[0,1].scatter(date_gamma, data_gamma, color =  colores[2])
    ax[0,1].scatter(date_mu, data_mu, color =  colores[3])
    ax[0,1].scatter(date_delta, data_delta, color =  colores[4])
    ax[0,1].scatter(date_omicron, data_omicron, color =  colores[5])
    ax[0,1].scatter(date_BA4_omicron, data_BA4_omicron, color =  colores[6])

    mask = df_fit['intervention_id'] == scenario
    original_fit_low = df_fit[mask]['DomOriginal_low'].to_numpy()
    original_fit_high = df_fit[mask]['DomOriginal_high'].to_numpy() 
    original_fit_median = df_fit[mask]['DomOriginal_median'].to_numpy()

    alpha_fit_low = df_fit[mask]['DomAlpha_low'].to_numpy()
    alpha_fit_high = df_fit[mask]['DomAlpha_high'].to_numpy() 
    alpha_fit_median = df_fit[mask]['DomAlpha_median'].to_numpy()

    gamma_fit_low = df_fit[mask]['DomGamma_low'].to_numpy()
    gamma_fit_high = df_fit[mask]['DomGamma_high'].to_numpy() 
    gamma_fit_median = df_fit[mask]['DomGamma_median'].to_numpy()

    mu_fit_low = df_fit[mask]['DomKappa_low'].to_numpy()
    mu_fit_high = df_fit[mask]['DomKappa_high'].to_numpy() 
    mu_fit_median = df_fit[mask]['DomKappa_median'].to_numpy()

    delta_fit_low = df_fit[mask]['DomDelta_low'].to_numpy()
    delta_fit_high = df_fit[mask]['DomDelta_high'].to_numpy() 
    delta_fit_median = df_fit[mask]['DomDelta_median'].to_numpy()

    omicron_fit_low = df_fit[mask]['DomOmicron_low'].to_numpy()
    omicron_fit_high = df_fit[mask]['DomOmicron_high'].to_numpy() 
    omicron_fit_median = df_fit[mask]['DomOmicron_median'].to_numpy()

    date_original_fit = pd.to_datetime(df_fit[mask]['Date'].to_numpy(), format='%Y-%m-%d', errors='coerce')
    p1 = ax[0,1].plot(date_original_fit, original_fit_median,  alpha = 1, zorder = 2, color = colores[0], linewidth = 3)
    ax[0,1].fill_between(date_original_fit, original_fit_low, original_fit_high, where=original_fit_high >= original_fit_low, facecolor=colores[0], interpolate=True, alpha = 0.1)

    date_alpha_fit = pd.to_datetime(df_fit[mask]['Date'].to_numpy(), format='%Y-%m-%d', errors='coerce')
    p1 = ax[0,1].plot(date_alpha_fit, alpha_fit_median,  alpha = 1, zorder = 2, color = colores[1], linewidth = 3)
    ax[0,1].fill_between(date_alpha_fit, alpha_fit_low, alpha_fit_high, where=alpha_fit_high >= alpha_fit_low, facecolor=colores[1], interpolate=True, alpha = 0.1)

    date_gamma_fit = pd.to_datetime(df_fit[mask]['Date'].to_numpy(), format='%Y-%m-%d', errors='coerce')
    p2 = ax[0,1].plot(date_gamma_fit, gamma_fit_median,  alpha = 1, zorder = 2, color = colores[2], linewidth = 3)
    ax[0,1].fill_between(date_gamma_fit, gamma_fit_low, gamma_fit_high, where=gamma_fit_high >= gamma_fit_low, facecolor=colores[2], interpolate=True, alpha = 0.1)

    date_mu_fit = pd.to_datetime(df_fit[mask]['Date'].to_numpy(), format='%Y-%m-%d', errors='coerce')
    p3 = ax[0,1].plot(date_mu_fit, mu_fit_median,  alpha = 1, zorder = 2, color = colores[3], linewidth = 3)
    ax[0,1].fill_between(date_mu_fit, mu_fit_low, mu_fit_high, where=mu_fit_high >= mu_fit_low, facecolor=colores[3], interpolate=True, alpha = 0.1)

    date_delta_fit = pd.to_datetime(df_fit[mask]['Date'].to_numpy(), format='%Y-%m-%d', errors='coerce')
    p4 = ax[0,1].plot(date_delta_fit, delta_fit_median,  alpha = 1, zorder = 2, color = colores[4], linewidth = 3)
    ax[0,1].fill_between(date_delta_fit, delta_fit_low, delta_fit_high, where=delta_fit_high >= delta_fit_low, facecolor=colores[4], interpolate=True, alpha = 0.1)

    date_omicron_fit = pd.to_datetime(df_fit[mask]['Date'].to_numpy(), format='%Y-%m-%d', errors='coerce')
    p5 = ax[0,1].plot(date_omicron_fit, omicron_fit_median,  alpha = 1, zorder = 2, color = colores[5], linewidth = 3)
    ax[0,1].fill_between(date_omicron_fit, omicron_fit_low, omicron_fit_high, where=omicron_fit_high >= omicron_fit_low, facecolor=colores[5], interpolate=True, alpha = 0.1)

    ### Subplot details
    ax[0,1].set_xticks(pd.date_range(start_day, periods=70, freq='15d'))
    plt.setp(ax[0,1].get_xticklabels(), rotation=90, fontsize=10)
    ax[0,1].set_xlim(initial_day, forecast_day)
    ax[0,1].yaxis.set_ticks(pd.np.arange(0, 1.1, 0.1))
    ax[0,1].set_ylim(0, 1.05)
    ax[0,1].set_ylabel('Proporción variantes')
    ax[0,1].grid()

    #######################
    ######### AR ##########
    #######################
    AR_0_fit_low = df_fit[mask]['AR_0_low'].to_numpy()
    AR_0_fit_high = df_fit[mask]['AR_0_high'].to_numpy()
    AR_0_fit_median = df_fit[mask]['AR_0_median'].to_numpy()

    AR_1_fit_low = df_fit[mask]['AR_1_low'].to_numpy()
    AR_1_fit_high = df_fit[mask]['AR_1_high'].to_numpy()
    AR_1_fit_median = df_fit[mask]['AR_1_median'].to_numpy()

    AR_2_fit_low = df_fit[mask]['AR_2_low'].to_numpy()
    AR_2_fit_high = df_fit[mask]['AR_2_high'].to_numpy()
    AR_2_fit_median = df_fit[mask]['AR_2_median'].to_numpy()

    AR_3_fit_low = df_fit[mask]['AR_3_low'].to_numpy()
    AR_3_fit_high = df_fit[mask]['AR_3_high'].to_numpy()
    AR_3_fit_median = df_fit[mask]['AR_3_median'].to_numpy()

    AR_4_fit_low = df_fit[mask]['AR_4_low'].to_numpy()
    AR_4_fit_high = df_fit[mask]['AR_4_high'].to_numpy()
    AR_4_fit_median = df_fit[mask]['AR_4_median'].to_numpy()

    AR_5_fit_low = df_fit[mask]['AR_5_low'].to_numpy()
    AR_5_fit_high = df_fit[mask]['AR_5_high'].to_numpy()
    AR_5_fit_median = df_fit[mask]['AR_5_median'].to_numpy()

    date_AR_fit = pd.to_datetime(df_fit[mask]['Date'].to_numpy(), format='%Y-%m-%d', errors='coerce')
    p2 = ax[1,1].plot(date_AR_fit, AR_0_fit_median,  alpha = 1, zorder = 2, color = colores[0], linewidth = 2)
    ax[1,1].fill_between(date_AR_fit, AR_0_fit_low, AR_0_fit_high, where=AR_0_fit_high >= AR_0_fit_low, facecolor=colores[0], interpolate=True, alpha = 0.1)

    p3 = ax[1,1].plot(date_AR_fit, AR_1_fit_median,  alpha = 1, zorder = 2, color = colores[1], linewidth = 2)
    ax[1,1].fill_between(date_AR_fit, AR_1_fit_low, AR_1_fit_high, where=AR_1_fit_high >= AR_1_fit_low, facecolor=colores[1], interpolate=True, alpha = 0.1)

    p4 = ax[1,1].plot(date_AR_fit, AR_2_fit_median,  alpha = 1, zorder = 2, color = colores[2], linewidth = 2)
    ax[1,1].fill_between(date_AR_fit, AR_2_fit_low, AR_2_fit_high, where=AR_2_fit_high >= AR_2_fit_low, facecolor=colores[2], interpolate=True, alpha = 0.1)

    p5 = ax[1,1].plot(date_AR_fit, AR_3_fit_median,  alpha = 1, zorder = 2, color = colores[3], linewidth = 2)
    ax[1,1].fill_between(date_AR_fit, AR_3_fit_low, AR_3_fit_high, where=AR_3_fit_high >= AR_3_fit_low, facecolor=colores[3], interpolate=True, alpha = 0.1)

    p6 = ax[1,1].plot(date_AR_fit, AR_4_fit_median,  alpha = 1, zorder = 2, color = colores[4], linewidth = 2)
    ax[1,1].fill_between(date_AR_fit, AR_4_fit_low, AR_4_fit_high, where=AR_4_fit_high >= AR_4_fit_low, facecolor=colores[4], interpolate=True, alpha = 0.1)

    p7 = ax[1,1].plot(date_AR_fit, AR_5_fit_median,  alpha = 1, zorder = 2, color = colores[5], linewidth = 2)
    ax[1,1].fill_between(date_AR_fit, AR_5_fit_low, AR_5_fit_high, where=AR_5_fit_high >= AR_5_fit_low, facecolor=colores[5], interpolate=True, alpha = 0.1)

    ### Subplot details
    ax[1,1].set_xticks(pd.date_range(start_day, periods=70, freq='15d'))
    plt.setp(ax[1,1].get_xticklabels(), rotation=90, fontsize=10)
    ax[1,1].set_xlim(initial_day, forecast_day)
    ax[1,1].set_ylabel('Tasa de ataque')
    ax[1,1].grid() 

    ##################
    ######  Rt  ######
    ##################
    mask_Rt = data_Rt['Localidad'] == 'Bogota D.C.'
    dias_Rt = pd.to_datetime(data_Rt[mask_Rt]['Fecha_Fin_Ventana'], format='%d/%m/%Y', errors='coerce').to_numpy()
    datos_Rt = data_Rt[mask_Rt]['Mean(R)']

    ## Rt reportado
    pX = ax[1,0].plot(dias_Rt, datos_Rt, marker='o', label='Rt reportado', markeredgecolor='slategray', markerfacecolor='none', alpha = 0.9, zorder = 1, color = 'slategray', linewidth = 2, markersize = 5)
    #P.append((pX[0],))
    #legend_label.append('R(t) reportado')

    Rt_total_median = df_fit[mask]['RR_total'].to_numpy()
    Rt_median = df_fit[mask]['RR_median'].to_numpy()
    Rt_1_median = df_fit[mask]['RR_1_median'].to_numpy()
    Rt_2_median = df_fit[mask]['RR_2_median'].to_numpy()
    Rt_3_median = df_fit[mask]['RR_3_median'].to_numpy()
    Rt_4_median = df_fit[mask]['RR_4_median'].to_numpy()
    Rt_5_median = df_fit[mask]['RR_5_median'].to_numpy()

    date_Rt_fit = pd.to_datetime(df_fit[mask]['Date'].to_numpy(), format='%Y-%m-%d', errors='coerce')

    p1 = ax[1,0].plot(date_Rt_fit, Rt_total_median, label = 'Rt modelo',  alpha = 0.8, zorder = 2, color = 'black', linewidth = 2)
    ax[1,0].legend()

    p2 = ax[1,0].plot(date_Rt_fit, Rt_median,  alpha = 0.4, zorder = 1, color = colores[0], linewidth = 2)
    p3 = ax[1,0].plot(date_Rt_fit, Rt_1_median,  alpha = 0.4, zorder = 1, color = colores[1], linewidth = 2)
    p4 = ax[1,0].plot(date_Rt_fit, Rt_2_median,  alpha = 0.4, zorder = 1, color = colores[2], linewidth = 2)
    p5 = ax[1,0].plot(date_Rt_fit, Rt_3_median,  alpha = 0.4, zorder = 1, color = colores[3], linewidth = 2)
    p6 = ax[1,0].plot(date_Rt_fit, Rt_4_median,  alpha = 0.4, zorder = 1, color = colores[4], linewidth = 2)
    p7 = ax[1,0].plot(date_Rt_fit, Rt_5_median,  alpha = 0.4, zorder = 1, color = colores[5], linewidth = 2)

    ### Subplot details
    ax[1,0].set_xticks(pd.date_range(start_day, periods=70, freq='15d'))
    plt.setp(ax[1,0].get_xticklabels(), rotation=90, fontsize=10)
    ax[1,0].set_xlim(initial_day, forecast_day)
    ax[1,0].set_ylim(0, 3.5)
    ax[1,0].set_ylabel('Rt')
    ax[1,0].grid() 

    ######################
    #### Plot details ####
    ######################
    p1 = ax[0,0].plot(pd.np.NaN, pd.np.NaN, color = colores[0], linewidth = 6)
    p2 = ax[0,0].plot(pd.np.NaN, pd.np.NaN, color = colores[1], linewidth = 6)
    p3 = ax[0,0].plot(pd.np.NaN, pd.np.NaN, color = colores[2], linewidth = 6)
    p4 = ax[0,0].plot(pd.np.NaN, pd.np.NaN, color = colores[3], linewidth = 6)
    p5 = ax[0,0].plot(pd.np.NaN, pd.np.NaN, color = colores[4], linewidth = 6)
    p6 = ax[0,0].plot(pd.np.NaN, pd.np.NaN, color = colores[5], linewidth = 6)

    P.append((p1[0],))
    P.append((p2[0],))
    P.append((p3[0],))
    P.append((p4[0],))
    P.append((p5[0],))
    P.append((p6[0],))

    legend_label.append('Original')
    legend_label.append('Alpha')
    legend_label.append('Gamma')
    legend_label.append('Mu')
    legend_label.append('Delta')
    legend_label.append('Omicron')

    plt.yticks(fontsize=14)
    fig1.legend(P, legend_label, loc='upper center', bbox_to_anchor=(0.5, 1.05),
                fancybox=True, shadow=True, ncol=7, fontsize = 13.5)

    plt.savefig('../../figures/'+projection_label+'-'+scenario+'.png', bbox_inches='tight', pad_inches=0.03)