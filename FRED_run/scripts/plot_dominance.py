import pandas as pd
import matplotlib.pyplot as plt
import warnings
import matplotlib.colors as mcolors
import matplotlib.transforms as mtransforms
import datetime as dt
import subprocess

TEST = 'omicron_eq_delta_vax_out'

colores = []
for name, color in mcolors.TABLEAU_COLORS.items():
    colores.append(name[4:])
warnings.filterwarnings('ignore')

subprocess.run("Rscript ~/FRED_Implementation/fit_scripts/fit_test.R test_"+TEST, shell=True, check=True)

file = '../../output/fred_output_model_fit_TEST.csv'
df_fit = pd.read_csv(file)


## PLOT
plt.rcParams.update({'font.size': 20})
fig1, ax1 = plt.subplots()
fig1.set_size_inches(w=17, h=7)

P = []
legend_label = []

file = '../../input_files/Bogota_Covid_Variants_Dominance_IC.csv'
df_dom_data = pd.read_csv(file)

mask = df_dom_data['var'] == 'Alpha'
date_alpha = pd.to_datetime(df_dom_data[mask]['date'].to_numpy(), format='%m/%d/%Y', errors='coerce')
data_alpha = df_dom_data[mask]['PointEst'].to_numpy()
ax1.scatter(date_alpha, data_alpha, color =  colores[0])

mask = df_dom_data['var'] == 'Gamma'
date_alpha = pd.to_datetime(df_dom_data[mask]['date'].to_numpy(), format='%m/%d/%Y', errors='coerce')
data_alpha = df_dom_data[mask]['PointEst'].to_numpy()
ax1.scatter(date_alpha, data_alpha, color =  colores[1])

mask = df_dom_data['var'] == 'Mu'
date_alpha = pd.to_datetime(df_dom_data[mask]['date'].to_numpy(), format='%m/%d/%Y', errors='coerce')
data_alpha = df_dom_data[mask]['PointEst'].to_numpy()
ax1.scatter(date_alpha, data_alpha, color =  colores[2])

mask = df_dom_data['var'] == 'Delta'
date_alpha = pd.to_datetime(df_dom_data[mask]['date'].to_numpy(), format='%m/%d/%Y', errors='coerce')
data_alpha = df_dom_data[mask]['PointEst'].to_numpy()
ax1.scatter(date_alpha, data_alpha, color =  colores[3])

mask = df_dom_data['var'] == 'Omicron'
date_alpha = pd.to_datetime(df_dom_data[mask]['date'].to_numpy(), format='%m/%d/%Y', errors='coerce')
data_alpha = df_dom_data[mask]['PointEst'].to_numpy()
ax1.scatter(date_alpha, data_alpha, color =  colores[4])

scenario = 'omicron-eq-delta-vax_1-mov_base'
counter = 0

mask = df_fit['intervention_id'] == scenario

alpha_fit_low = df_fit[mask]['DomAlpha_low'].to_numpy()
alpha_fit_high = df_fit[mask]['DomAlpha_high'].to_numpy() 
alpha_fit_median = df_fit[mask]['DomAlpha_median'].to_numpy()

date_alpha_fit = pd.to_datetime(df_fit[mask]['Date'].to_numpy(), format='%Y-%m-%d', errors='coerce')
p1 = ax1.plot(date_alpha_fit, alpha_fit_median,  alpha = 1, zorder = 2, color = colores[counter], linewidth = 3)
ax1.fill_between(date_alpha_fit, alpha_fit_low, alpha_fit_high, where=alpha_fit_high >= alpha_fit_low, facecolor=colores[counter], interpolate=True, alpha = 0.1)

gamma_fit_low = df_fit[mask]['DomGamma_low'].to_numpy()
gamma_fit_high = df_fit[mask]['DomGamma_high'].to_numpy() 
gamma_fit_median = df_fit[mask]['DomGamma_median'].to_numpy()

date_gamma_fit = pd.to_datetime(df_fit[mask]['Date'].to_numpy(), format='%Y-%m-%d', errors='coerce')
p2 = ax1.plot(date_gamma_fit, gamma_fit_median,  alpha = 1, zorder = 2, color = colores[counter+1], linewidth = 3)
ax1.fill_between(date_gamma_fit, gamma_fit_low, gamma_fit_high, where=gamma_fit_high >= gamma_fit_low, facecolor=colores[counter+1], interpolate=True, alpha = 0.1)

mu_fit_low = df_fit[mask]['DomKappa_low'].to_numpy()
mu_fit_high = df_fit[mask]['DomKappa_high'].to_numpy() 
mu_fit_median = df_fit[mask]['DomKappa_median'].to_numpy()

date_mu_fit = pd.to_datetime(df_fit[mask]['Date'].to_numpy(), format='%Y-%m-%d', errors='coerce')
p3 = ax1.plot(date_mu_fit, mu_fit_median,  alpha = 1, zorder = 2, color = colores[counter+2], linewidth = 3)
ax1.fill_between(date_mu_fit, mu_fit_low, mu_fit_high, where=mu_fit_high >= mu_fit_low, facecolor=colores[counter+2], interpolate=True, alpha = 0.1)

delta_fit_low = df_fit[mask]['DomDelta_low'].to_numpy()
delta_fit_high = df_fit[mask]['DomDelta_high'].to_numpy() 
delta_fit_median = df_fit[mask]['DomDelta_median'].to_numpy()

date_delta_fit = pd.to_datetime(df_fit[mask]['Date'].to_numpy(), format='%Y-%m-%d', errors='coerce')
p4 = ax1.plot(date_delta_fit, delta_fit_median,  alpha = 1, zorder = 2, color = colores[counter+3], linewidth = 3)
ax1.fill_between(date_delta_fit, delta_fit_low, delta_fit_high, where=delta_fit_high >= delta_fit_low, facecolor=colores[counter+3], interpolate=True, alpha = 0.1)

omicron_fit_low = df_fit[mask]['DomOmicron_low'].to_numpy()
omicron_fit_high = df_fit[mask]['DomOmicron_high'].to_numpy() 
omicron_fit_median = df_fit[mask]['DomOmicron_median'].to_numpy()

date_omicron_fit = pd.to_datetime(df_fit[mask]['Date'].to_numpy(), format='%Y-%m-%d', errors='coerce')
p5 = ax1.plot(date_omicron_fit, omicron_fit_median,  alpha = 1, zorder = 2, color = colores[counter+4], linewidth = 3)
ax1.fill_between(date_omicron_fit, omicron_fit_low, omicron_fit_high, where=omicron_fit_high >= omicron_fit_low, facecolor=colores[counter+4], interpolate=True, alpha = 0.1)

p11 = ax1.fill(pd.np.NaN, pd.np.NaN, color = colores[counter], alpha=0.2)
p12 = ax1.fill(pd.np.NaN, pd.np.NaN, color = colores[counter+1], alpha=0.2)
p13 = ax1.fill(pd.np.NaN, pd.np.NaN, color = colores[counter+2], alpha=0.2)
p14 = ax1.fill(pd.np.NaN, pd.np.NaN, color = colores[counter+3], alpha=0.2)
p15 = ax1.fill(pd.np.NaN, pd.np.NaN, color = colores[counter+4], alpha=0.2)

P.append((p1[0], p11[0]))
P.append((p2[0], p12[0]))
P.append((p3[0], p13[0]))
P.append((p4[0], p14[0]))
P.append((p5[0], p15[0]))

legend_label.append('Alpha')
legend_label.append('Gamma')
legend_label.append('Mu')
legend_label.append('Delta')
legend_label.append('Omicron')
    

### Plot details
time_labels = pd.date_range(dt.datetime(2020,1,3), periods=130, freq='7d')
plt.xticks(time_labels, rotation=90, fontsize=13)
ax1.set_xlim(dt.datetime(2021,1,1),dt.datetime(2022,3,15))
ax1.yaxis.set_ticks(pd.np.arange(0, 1.1, 0.1))
ax1.set_ylim(0, 1.05)
ax1.set_ylabel('Proporción variantes')
plt.grid()
ax1.legend(P, legend_label, loc='upper center', bbox_to_anchor=(0.5, 1.12),
          fancybox=True, shadow=True, ncol=5, fontsize = 15)

plt.savefig('../figures/dominance_'+TEST+'.png', bbox_inches='tight', pad_inches=0.03)

########################
## AR
########################
plt.rcParams.update({'font.size': 20})
fig1, ax1 = plt.subplots()
fig1.set_size_inches(w=17, h=7)

counter = 0
P = []
legend_label = []

AR_fit_median = df_fit[mask]['AR_median'].to_numpy()
AR_1_fit_median = df_fit[mask]['AR_1_median'].to_numpy()
AR_2_fit_median = df_fit[mask]['AR_2_median'].to_numpy()
AR_3_fit_median = df_fit[mask]['AR_3_median'].to_numpy()
AR_4_fit_median = df_fit[mask]['AR_4_median'].to_numpy()
AR_5_fit_median = df_fit[mask]['AR_5_median'].to_numpy()

date_AR_fit = pd.to_datetime(df_fit[mask]['Date'].to_numpy(), format='%Y-%m-%d', errors='coerce')

p2 = ax1.plot(date_AR_fit, AR_fit_median,  alpha = 1, zorder = 2, color = colores[0], linewidth = 2)
p3 = ax1.plot(date_AR_fit, AR_1_fit_median,  alpha = 1, zorder = 2, color = colores[1], linewidth = 2)
p4 = ax1.plot(date_AR_fit, AR_2_fit_median,  alpha = 1, zorder = 2, color = colores[2], linewidth = 2)
p5 = ax1.plot(date_AR_fit, AR_3_fit_median,  alpha = 1, zorder = 2, color = colores[3], linewidth = 2)
p6 = ax1.plot(date_AR_fit, AR_4_fit_median,  alpha = 1, zorder = 2, color = colores[4], linewidth = 2)
p7 = ax1.plot(date_AR_fit, AR_5_fit_median,  alpha = 1, zorder = 2, color = colores[5], linewidth = 2)

P.append((p2[0],))
P.append((p3[0],))
P.append((p4[0],))
P.append((p5[0],))
P.append((p6[0],))
P.append((p7[0],))

legend_label.append('Original')
legend_label.append('Alpha')
legend_label.append('Gamma')
legend_label.append('Mu')
legend_label.append('Delta')
legend_label.append('Omicron')
        
### Plot details
time_labels = pd.date_range(dt.datetime(2020,1,3), periods=135, freq='7d')
plt.xticks(time_labels, rotation=90, fontsize=13)
ax1.set_xlim(dt.datetime(2021,3,1),dt.datetime(2022,3,1))
ax1.set_ylabel('Tasa de ataque')
#ax1.yaxis.set_ticks(pd.np.arange(0, 1.1, 0.1))
#ax1.set_ylim(0.48, 1.01)

## Area de validación
# trans = mtransforms.blended_transform_factory(ax1.transData, ax1.transAxes)
# ax1.axvline(x =projection_date, color = 'black', alpha = 0.6, linewidth = 2, linestyle = (0, (3, 1, 1, 1)), zorder = 7)
# ax1.axvline(x = fit_date, color = 'black', alpha = 0.6, linewidth = 2, linestyle =  (0, (3, 1, 1, 1)), zorder = 7)

# ax1.fill_between(date_AR_fit, 0, 1, where = (date_AR_fit >= fit_date) & (date_AR_fit <= projection_date),
#                  facecolor=plt.cm.Pastel1(7), alpha=0.5, transform=trans, zorder = 1)
# ax1.fill_between(date_AR_fit, 0, 1, where = (date_AR_fit <= fit_date),
#                  facecolor=plt.cm.Pastel1(2), alpha=0.5, transform=trans, zorder = 1)

## Legend
ax1.legend(P, legend_label, loc='upper center', bbox_to_anchor=(0.5, 1.145),
          fancybox=True, shadow=True, ncol=7, fontsize = 13.5)
# ax2 = ax1.twinx()
# q1 = ax2.bar([0], [0], ls=(0, (3, 1, 1, 1)), color=plt.cm.Pastel1(2), ec='black')
# q2 = ax2.bar([0], [0], ls=(0, (3, 1, 1, 1)), color=plt.cm.Pastel1(7), ec='black')
# ax2.get_yaxis().set_visible(False)
# ax2.legend([(q1[0]),(q2[0])], ['Calibración', 'Validación'], loc='upper center', bbox_to_anchor=(0.5, -0.18),
#           fancybox=True, shadow=True, ncol=7, fontsize = 14)
ax1.grid() 
plt.savefig('../figures/AR_'+TEST+'.png', bbox_inches='tight', pad_inches=0.03)
