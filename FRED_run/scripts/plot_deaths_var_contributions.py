import pandas as pd
import matplotlib.pyplot as plt
import warnings
import matplotlib.colors as mcolors
import matplotlib.transforms as mtransforms
import datetime as dt

TEST = 'omicron_eq_delta_no_cross_no_vax_out'

colores = []
for name, color in mcolors.TABLEAU_COLORS.items():
    colores.append(name[4:])
warnings.filterwarnings('ignore')

### Muertes diarias
file_data = '../data/COL_covid_death_data.csv'
df_data = pd.read_csv(file_data, encoding = 'latin')

folder_name = 'FRED_11001_projections_asymp_1.00_fm_0.73_ksus_10.00_var_1_vax_070_mov_test_'+TEST

## Simulations
file = '../output/SHORT_FORECAST/'+folder_name+'/fred_output.csv'
df = pd.read_csv(file)

### Succesfull simulation parameters
file_data = '../output/SHORT_FORECAST/'+folder_name+'/FRED_parameters_out.csv'
df_params = pd.read_csv(file_data)

## PLOT
plt.rcParams.update({'font.size': 10})
fig, ax = plt.subplots()
fig.set_size_inches(w=15, h=7)

star_day = dt.datetime(2020,1,1)
death_data = df_data[df_data['MunCode'] == 11001]['Deaths'].to_numpy()

#############################################################
death_days = (pd.to_datetime(df_data[df_data['MunCode'] == 11001]['Date']))
##########################################################

no_sce = 0
sce_vec = ['omicron-eq-delta-vax_1-mov_base']

scenario = sce_vec[no_sce]
mask = df_params['intervention_id'] == scenario
df_scenario = df_params[mask]

for simul_id in df_scenario['job_id'].to_numpy():
    mask = df['job_id'] == simul_id
    max_df = df[mask]

    days = max_df['Day']
    star_day = dt.datetime(2020,1,1)
    date_days = [star_day + dt.timedelta(days = i) for i in days]

    deaths_0 = max_df['CF_mean'].to_numpy()
    deaths_1 = max_df['CF_1_mean'].to_numpy()
    deaths_2 = max_df['CF_2_mean'].to_numpy()
    deaths_3 = max_df['CF_3_mean'].to_numpy()
    deaths_4 = max_df['CF_4_mean'].to_numpy()
    deaths_5 = max_df['CF_5_mean'].to_numpy()
    
    date_days = pd.np.array(date_days)
        
    ### Simulation
    ax.plot(date_days, deaths_0, alpha = 1, zorder = 1, color = colores[0])
    ax.plot(date_days, deaths_1, alpha = 1, zorder = 1, color = colores[1])
    ax.plot(date_days, deaths_2, alpha = 1, zorder = 1, color = colores[2])
    ax.plot(date_days, deaths_3, alpha = 1, zorder = 1, color = colores[3])
    ax.plot(date_days, deaths_4, alpha = 1, zorder = 1, color = colores[4])
    ax.plot(date_days, deaths_5, alpha = 1, zorder = 1, color = colores[5])

P = []
legend_label = []

p1 = ax.plot(pd.np.NaN, pd.np.NaN, color = colores[0])
p2 = ax.plot(pd.np.NaN, pd.np.NaN, color = colores[1])
p3 = ax.plot(pd.np.NaN, pd.np.NaN, color = colores[2])
p4 = ax.plot(pd.np.NaN, pd.np.NaN, color = colores[3])
p5 = ax.plot(pd.np.NaN, pd.np.NaN, color = colores[4])
p6 = ax.plot(pd.np.NaN, pd.np.NaN, color = colores[5])
    
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


ax.scatter(death_days, death_data, s = 10, c = 'black', label = 'Data-deaths', zorder = 2)

plt.xticks(pd.date_range(star_day, periods=70, freq='15d'), rotation=90, fontsize=12)
plt.yticks(fontsize=14)
ax.set_xlim(dt.datetime(2020,6,1),dt.datetime(2022,3,31))
ax.grid()

ax.legend(P, legend_label, loc='upper center', bbox_to_anchor=(0.5, 1.145),
              fancybox=True, shadow=True, ncol=7, fontsize = 13.5)

plt.savefig('../figures/deaths_var_contributions'+TEST+'.png', bbox_inches='tight', pad_inches=0.03)
