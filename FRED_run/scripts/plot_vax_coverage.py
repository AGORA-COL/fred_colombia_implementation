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

num_grupo = {1:'0-4',2:'5-9',3:'10-14',4:'15-19',5:'20-24',6:'25-29',7:'30-34',
             8:'35-39',9:'40-44',10:'45-49',11:'50-54',12:'55-59',13:'60-64',
             14:'65-69',15:'70-74',16:'75-79',17:'80-120'}

######################
### Simul settings ###
######################
forecast_day = dt.datetime(2022,4,1)
initial_day = dt.datetime(2020,6,1)

test_label = 'omicron_lineages'
folder_name = 'FRED_11001_projections_asymp_1.00_fm_0.73_ksus_10.00_var_1_vax_070_mov_'+test_label+'_out'
scenario = 'omicronBAX_X4_omicron-cross-070-transexc-1_2-severity-28-vax_1-mov_open_school'
title = scenario

######################
##### Load info ######
######################
## Simulations
file = '../../output/SHORT_FORECAST/'+folder_name+'/fred_output.csv'
df = pd.read_csv(file)
mask = df['job_id'] == df['job_id'].unique()[0]

## Date
fec_ini='2020/01/01'
days = df[mask].Day
simul_dates = [pd.to_datetime(fec_ini) + dt.timedelta(days=x) for x in days]

## Age groups
age_group = []
age_group_1 = []
for i in range(16):
    age_group.append(str(5*i)+'_'+str(5*(i+1) - 1))
    age_group_1.append(str(5*i)+'_'+str(5*(i+1)))
age_group.append('80_120')
age_group_1.append('80_120')

df_empt = pd.DataFrame(columns = ['Fecha', 'Vac', 'Edad'])
for variable in ['V1A','V2A','V3A']:
    counter = 1
    df_vax_map = df_empt.copy()
    df_aux = df_empt.copy()
    total_vax = pd.np.zeros(len(simul_dates))
    c = 0
    for age in age_group:
        df_aux_1 = df_aux 
        df_aux_1['Edad'] = counter*pd.np.ones(len(simul_dates))
        df_aux_1['Fecha'] = simul_dates
        df_aux_1['Vac'] = df[mask][variable+'_'+age+'_mean']
        df_vax_map = pd.concat([df_vax_map, df_aux_1])
        total_vax += df_aux_1['Vac'].to_numpy()
        counter += 1
        c += 1
    if variable == 'V1A':
        df_vax_map_dosis_1 = df_vax_map
        total_vax_1 = total_vax
    elif variable == 'V2A':
        df_vax_map_dosis_2 = df_vax_map
        total_vax_2 = total_vax
    elif variable == 'V3A':
        df_vax_map_dosis_refuerzo = df_vax_map
        total_vax_refuerzo = total_vax

## FRED pob info
pob_fred = []
for age in age_group_1:
    pob_fred.append(df[mask]['AgeNA'+age].unique()[0])

## FRED vax info
df_cobertura_FRED_dosis_1 = df_vax_map_dosis_1.groupby(by=['Edad']).sum().reset_index()
df_cobertura_FRED_dosis_1['pob'] = pob_fred
df_cobertura_FRED_dosis_1['cobertura'] = df_cobertura_FRED_dosis_1['Vac']/df_cobertura_FRED_dosis_1['pob']
df_cobertura_FRED_dosis_1.columns = ['cod_edad','Vac','pob','cobertura']

df_cobertura_FRED_dosis_2 = df_vax_map_dosis_2.groupby(by=['Edad']).sum().reset_index()
df_cobertura_FRED_dosis_2['pob'] = pob_fred
df_cobertura_FRED_dosis_2['cobertura'] = df_cobertura_FRED_dosis_2['Vac']/df_cobertura_FRED_dosis_2['pob']
df_cobertura_FRED_dosis_2.columns = ['cod_edad','Vac','pob','cobertura']

df_cobertura_FRED_dosis_refuerzo = df_vax_map_dosis_refuerzo.groupby(by=['Edad']).sum().reset_index()
df_cobertura_FRED_dosis_refuerzo['pob'] = pob_fred
df_cobertura_FRED_dosis_refuerzo['cobertura'] = df_cobertura_FRED_dosis_refuerzo['Vac']/df_cobertura_FRED_dosis_refuerzo['pob']
df_cobertura_FRED_dosis_refuerzo.columns = ['cod_edad','Vac','pob','cobertura']

## Load real vax info
df_primera_dosis = pd.read_csv('../../data/primera_dosis.csv')
df_primera_dosis['Fecha'] = pd.to_datetime(df_primera_dosis['Fecha'])
#del(df_primera_dosis['Unnamed: 0'])

df_segunda_dosis = pd.read_csv('../../data/segunda_dosis.csv')
df_segunda_dosis['Fecha'] = pd.to_datetime(df_segunda_dosis['Fecha'])
#del(df_segunda_dosis['Unnamed: 0'])

df_refuerzo_dosis = pd.read_csv('../../data/refuerzo_dosis.csv')
df_refuerzo_dosis['Fecha'] = pd.to_datetime(df_refuerzo_dosis['Fecha'])
#del(df_refuerzo_dosis['Unnamed: 0'])

#################################################################################
### Primera dosis
#################################################################################
for i in range(16):
    mask = (df_primera_dosis['Edad'] >= 5*i) & (df_primera_dosis['Edad'] < 5*(i+1))
    df_primera_dosis.loc[mask,'age_group'] = str(5*i)+'-'+str(5*(i+1) - 1)

mask = (df_primera_dosis['Edad'] >= 80)
df_primera_dosis.loc[mask,'age_group'] = '80-120'

df_primera_dosis = df_primera_dosis.groupby(by=['Fecha', 'age_group']).sum().reset_index()[['Fecha','age_group','Vacunas']]

counter = 1
for i in range(16):
    mask = df_primera_dosis['age_group'] == str(5*i)+'-'+str(5*(i+1) - 1)
    df_primera_dosis.loc[mask,'color'] = counter
    counter += 1

mask = df_primera_dosis['age_group'] == '80-120'
df_primera_dosis.loc[mask,'color'] = counter

total_primera_dosis = df_primera_dosis.groupby(by=['Fecha']).sum().reset_index()

#################################################################################
### Segunda dosis
#################################################################################
for i in range(16):
    mask = (df_segunda_dosis['Edad'] >= 5*i) & (df_segunda_dosis['Edad'] < 5*(i+1))
    df_segunda_dosis.loc[mask,'age_group'] = str(5*i)+'-'+str(5*(i+1) - 1)

mask = (df_segunda_dosis['Edad'] >= 80)
df_segunda_dosis.loc[mask,'age_group'] = '80-120'

df_segunda_dosis = df_segunda_dosis.groupby(by=['Fecha', 'age_group']).sum().reset_index()[['Fecha','age_group','Vacunas']]

counter = 1
for i in range(16):
    mask = df_segunda_dosis['age_group'] == str(5*i)+'-'+str(5*(i+1) - 1)
    df_segunda_dosis.loc[mask,'color'] = counter
    counter += 1

mask = df_segunda_dosis['age_group'] == '80-120'
df_segunda_dosis.loc[mask,'color'] = counter

total_segunda_dosis = df_segunda_dosis.groupby(by=['Fecha']).sum().reset_index()

#################################################################################
### Refuerzo dosis
#################################################################################
for i in range(16):
    mask = (df_refuerzo_dosis['Edad'] >= 5*i) & (df_refuerzo_dosis['Edad'] < 5*(i+1))
    df_refuerzo_dosis.loc[mask,'age_group'] = str(5*i)+'-'+str(5*(i+1) - 1)

mask = (df_refuerzo_dosis['Edad'] >= 80)
df_refuerzo_dosis.loc[mask,'age_group'] = '80-120'

df_refuerzo_dosis = df_refuerzo_dosis.groupby(by=['Fecha', 'age_group']).sum().reset_index()[['Fecha','age_group','Vacunas']]

counter = 1
for i in range(16):
    mask = df_refuerzo_dosis['age_group'] == str(5*i)+'-'+str(5*(i+1) - 1)
    df_refuerzo_dosis.loc[mask,'color'] = counter
    counter += 1

mask = df_refuerzo_dosis['age_group'] == '80-120'
df_refuerzo_dosis.loc[mask,'color'] = counter

total_refuerzo_dosis = df_refuerzo_dosis.groupby(by=['Fecha']).sum().reset_index()

### Plot total first sum
plt.rcParams.update({'font.size': 15})
fig1, ax1 = plt.subplots()
fig1.set_size_inches(w=18, h=7)

P = []
legend_label = []
##################################################################################################
p1 = ax1.scatter(total_primera_dosis['Fecha'], total_primera_dosis['Vacunas'], color = 'orange')
P.append((p1))
legend_label.append('Datos')

p2 = ax1.plot(simul_dates, total_vax_1, alpha = 1)
P.append((p2[0]))
legend_label.append('FRED')
##################################################################################################
ax1.legend(P, legend_label, loc='upper center', bbox_to_anchor=(0.5, 1.12),
              fancybox=True, shadow=True, ncol=5, fontsize = 13.5)

ax1.set_xlim(initial_day, forecast_day)
ax1.set_ylabel('Primera dosis')
ax1.grid()
plt.savefig('../../figures/total_vax_1_'+test_label+'.png', bbox_inches='tight', pad_inches=0.03)

### Plot total second sum
plt.rcParams.update({'font.size': 15})
fig1, ax1 = plt.subplots()
fig1.set_size_inches(w=18, h=7)

P = []
legend_label = []
##################################################################################################
p1 = ax1.scatter(total_segunda_dosis['Fecha'], total_segunda_dosis['Vacunas'], color = 'orange')
P.append((p1))
legend_label.append('Datos')

p2 = ax1.plot(simul_dates, total_vax_2, alpha = 1)
P.append((p2[0]))
legend_label.append('FRED')
##################################################################################################
ax1.legend(P, legend_label, loc='upper center', bbox_to_anchor=(0.5, 1.12),
              fancybox=True, shadow=True, ncol=5, fontsize = 13.5)

ax1.set_xlim(initial_day, forecast_day)
ax1.set_ylabel('Segunda dosis')
ax1.grid()
plt.savefig('../../figures/total_vax_2_'+test_label+'.png', bbox_inches='tight', pad_inches=0.03)

### Plot total third sum
plt.rcParams.update({'font.size': 15})
fig1, ax1 = plt.subplots()
fig1.set_size_inches(w=18, h=7)

P = []
legend_label = []
##################################################################################################
p1 = ax1.scatter(total_refuerzo_dosis['Fecha'], total_refuerzo_dosis['Vacunas'], color = 'orange')
P.append((p1))
legend_label.append('Datos')

p2 = ax1.plot(simul_dates, total_vax_refuerzo, alpha = 1)
P.append((p2[0]))
legend_label.append('FRED')
##################################################################################################
ax1.legend(P, legend_label, loc='upper center', bbox_to_anchor=(0.5, 1.12),
              fancybox=True, shadow=True, ncol=5, fontsize = 13.5)

ax1.set_xlim(initial_day, forecast_day)
ax1.set_ylabel('Refuerzo dosis')
ax1.grid()
plt.savefig('../../figures/total_vax_refuerzo_'+test_label+'.png', bbox_inches='tight', pad_inches=0.03)

################
### Pob info ###
################
df_pob = pd.read_csv('../../data/pob_edad_loc.csv', encoding = 'latin')
df_pob_aux = df_pob.loc[(df_pob['AREA'] == 'Total') & (df_pob['AÑO'] == 2020)]

pob_census = []
for age in age_group[:-1]:
    pob_census.append(df_pob_aux['Total_'+age.replace('_','-')].str.replace(',','').astype(int).sum())
    
pob_80_120 = 0
pob_80_84 = df_pob_aux['Total_80-84'].str.replace(',','').astype(int).sum()
pob_85_89 = df_pob_aux['Total_80-84'].str.replace(',','').astype(int).sum()
pob_90_94 = df_pob_aux['Total_80-84'].str.replace(',','').astype(int).sum()
pob_95_99 = df_pob_aux['Total_80-84'].str.replace(',','').astype(int).sum()
pob_100 = df_pob_aux['Total_80-84'].str.replace(',','').astype(int).sum()

for i in range(16,20):
    pob_80_120 += df_pob_aux['Total_'+str(5*i)+'-'+str(5*(i+1) - 1)].str.replace(',','').astype(int).sum()
pob_80_120 += df_pob_aux['Total_100+'].str.replace(',','').astype(int).sum()

pob_census.append(pob_80_120)

#################
### FRED info ###
#################
df_cobertura = df_primera_dosis.groupby(by=['color']).sum().reset_index()
df_cobertura['pob'] = pob_census
df_cobertura['Vacunas'] = df_cobertura['Vacunas'].fillna(0)
df_cobertura['cobertura'] = df_cobertura['Vacunas']/df_cobertura['pob']
df_cobertura.columns = ['cod_edad','Vac','pob','cobertura']

df_cobertura_comp = pd.DataFrame({})
df_cobertura_comp['Edad'] = age_group
df_cobertura_comp['cober_datos'] = df_cobertura['cobertura']
df_cobertura_comp['cober_FRED'] = df_cobertura_FRED_dosis_1['cobertura']
df_cobertura_comp['diff_porc'] = 100*pd.np.abs(df_cobertura_FRED_dosis_1['cobertura'] - df_cobertura['cobertura'])/df_cobertura['cobertura']

df_cobertura_comp_primera = df_cobertura_comp
df_primera_dosis_acum_FRED = df_cobertura_FRED_dosis_1
df_primera_dosis_acum = df_cobertura

df_cobertura = pd.DataFrame({})
df_cobertura['color'] = pd.np.arange(1,18,1)
df_cobertura = pd.merge(df_cobertura, df_segunda_dosis.groupby(by=['color']).sum().reset_index(), how='left', on='color')
df_cobertura['pob'] = pob_census
df_cobertura['Vacunas'] = df_cobertura['Vacunas'].fillna(0)
df_cobertura['cobertura'] = df_cobertura['Vacunas']/df_cobertura['pob']
df_cobertura.columns = ['cod_edad','Vac','pob','cobertura']

df_cobertura_comp = pd.DataFrame({})
df_cobertura_comp['Edad'] = age_group
df_cobertura_comp['cober_datos'] = df_cobertura['cobertura']
df_cobertura_comp['cober_FRED'] = df_cobertura_FRED_dosis_2['cobertura']
df_cobertura_comp['diff_porc'] = 100*pd.np.abs(df_cobertura_FRED_dosis_2['cobertura'] - df_cobertura['cobertura'])/df_cobertura['cobertura']

df_cobertura_comp_segunda = df_cobertura_comp
df_segunda_dosis_acum_FRED = df_cobertura_FRED_dosis_2
df_segunda_dosis_acum = df_cobertura

df_cobertura = pd.DataFrame({})
df_cobertura['color'] = pd.np.arange(1,18,1)
df_cobertura = pd.merge(df_cobertura, df_refuerzo_dosis.groupby(by=['color']).sum().reset_index(), how='left', on='color')
df_cobertura['pob'] = pob_census
df_cobertura['Vacunas'] = df_cobertura['Vacunas'].fillna(0)
df_cobertura['cobertura'] = df_cobertura['Vacunas']/df_cobertura['pob']
df_cobertura.columns = ['cod_edad','Vac','pob','cobertura']

df_cobertura_comp = pd.DataFrame({})
df_cobertura_comp['Edad'] = age_group
df_cobertura_comp['cober_datos'] = df_cobertura['cobertura']
df_cobertura_comp['cober_FRED'] = df_cobertura_FRED_dosis_2['cobertura']
df_cobertura_comp['diff_porc'] = 100*pd.np.abs(df_cobertura_FRED_dosis_2['cobertura'] - df_cobertura['cobertura'])/df_cobertura['cobertura']

df_cobertura_comp_refuerzo = df_cobertura_comp
df_refuerzo_dosis_acum_FRED = df_cobertura_FRED_dosis_2
df_refuerzo_dosis_acum = df_cobertura

## Ticks
ticks = []
for age in age_group_1:
    ticks.append(age.replace('_','-'))

###############################
### Cobertura primera dosis ###
###############################

plt.rcParams.update({'font.size': 15})
fig, ax = plt.subplots(5,4, sharex=True, tight_layout=True)
fig.suptitle('Coberturas primera dosis', y=-0.04, fontsize= 25)
fig.set_size_inches(w=17, h=17)

color_1 = 'blue'
color_2 = 'purple'

counter = 3
for grupo in range(17,0,-1):
    df_vax_map_dosis_1['Vac'] = df_vax_map_dosis_1['Vac'].fillna(0)
    df_vax_map_aux = df_vax_map_dosis_1.loc[df_vax_map_dosis_1['Edad'] == grupo]
    df_primera_dosis_aux = df_primera_dosis.loc[df_primera_dosis['color'] == grupo]
    q1 = ax[int(pd.np.floor(counter/4)),counter%4].scatter(df_vax_map_aux['Fecha'], pd.np.add.accumulate(df_vax_map_aux['Vac'])/df_cobertura_FRED_dosis_1['pob'].to_numpy()[grupo-1], s = 20, color = color_1)
    q2 = ax[int(pd.np.floor(counter/4)),counter%4].scatter(df_primera_dosis_aux['Fecha'], pd.np.add.accumulate(df_primera_dosis_aux['Vacunas'].to_numpy())/df_cobertura['pob'].to_numpy()[grupo-1], s = 20, color = color_2)
    ax[int(pd.np.floor(counter/4)),counter%4].grid()
    ax[int(pd.np.floor(counter/4)),counter%4].set_title(num_grupo[grupo].replace('_','-'))
    plt.setp(ax[int(pd.np.floor(counter/4)),counter%4].get_xticklabels(), rotation=90)
    ax[int(pd.np.floor(counter/4)),counter%4].set_xlim(dt.datetime(2021,1,1),dt.datetime(2022,3,1))
    counter += 1
    
fig.legend([(q1),(q2)], ['FRED', 'Datos'], loc='lower center', bbox_to_anchor=(0.5,-0.03),
          fancybox=True, shadow=True, ncol=5, fontsize = 14)
    
fig.delaxes(ax[0,0])
fig.delaxes(ax[0,1])
fig.delaxes(ax[0,2])

plt.savefig('../../figures/cobertura_primera_dosis.png', bbox_inches='tight', pad_inches=0.03)

###############################
### Cobertura segunda dosis ###
###############################

plt.rcParams.update({'font.size': 15})
fig, ax = plt.subplots(5,4, sharex=True, tight_layout=True)
fig.suptitle('Coberturas segunda dosis', y=-0.04, fontsize= 25)
fig.set_size_inches(w=17, h=17)

color_1 = 'red'
color_2 = 'orange'

counter = 3
for grupo in range(17,0,-1):
    df_vax_map_dosis_2['Vac'] = df_vax_map_dosis_2['Vac'].fillna(0)
    df_vax_map_aux = df_vax_map_dosis_2.loc[df_vax_map_dosis_2['Edad'] == grupo]
    df_segunda_dosis_aux = df_segunda_dosis.loc[df_segunda_dosis['color'] == grupo]
    q1 = ax[int(pd.np.floor(counter/4)),counter%4].scatter(df_vax_map_aux['Fecha'], pd.np.add.accumulate(df_vax_map_aux['Vac'])/df_cobertura_FRED_dosis_2['pob'].to_numpy()[grupo-1], s = 20, color = color_1)
    q2 = ax[int(pd.np.floor(counter/4)),counter%4].scatter(df_segunda_dosis_aux['Fecha'], pd.np.add.accumulate(df_segunda_dosis_aux['Vacunas'].to_numpy())/df_cobertura['pob'].to_numpy()[grupo-1], s = 20, color = color_2)
    ax[int(pd.np.floor(counter/4)),counter%4].grid()
    ax[int(pd.np.floor(counter/4)),counter%4].set_title(num_grupo[grupo].replace('_','-'))
    plt.setp(ax[int(pd.np.floor(counter/4)),counter%4].get_xticklabels(), rotation=90)
    ax[int(pd.np.floor(counter/4)),counter%4].set_xlim(dt.datetime(2021,1,1),dt.datetime(2022,3,1))
    counter += 1
    
fig.legend([(q1),(q2)], ['FRED', 'Datos'], loc='lower center', bbox_to_anchor=(0.5,-0.03),
          fancybox=True, shadow=True, ncol=5, fontsize = 14)
    
fig.delaxes(ax[0,0])
fig.delaxes(ax[0,1])
fig.delaxes(ax[0,2])

plt.savefig('../../figures/cobertura_segunda_dosis.png', bbox_inches='tight', pad_inches=0.03)

###########################
### Cobertura refuerzo  ###
###########################

plt.rcParams.update({'font.size': 15})
fig, ax = plt.subplots(5,4, sharex=True, tight_layout=True)
fig.suptitle('Coberturas refuerzo dosis', y=-0.04, fontsize= 25)
fig.set_size_inches(w=17, h=17)

color_1 = 'green'
color_2 = 'yellow'

counter = 3
for grupo in range(17,0,-1):
    df_vax_map_dosis_refuerzo['Vac'] = df_vax_map_dosis_refuerzo['Vac'].fillna(0)
    df_vax_map_aux = df_vax_map_dosis_refuerzo.loc[df_vax_map_dosis_refuerzo['Edad'] == grupo]
    df_refuerzo_dosis_aux = df_refuerzo_dosis.loc[df_refuerzo_dosis['color'] == grupo]
    q1 = ax[int(pd.np.floor(counter/4)),counter%4].scatter(df_vax_map_aux['Fecha'], pd.np.add.accumulate(df_vax_map_aux['Vac'])/df_cobertura_FRED_dosis_refuerzo['pob'].to_numpy()[grupo-1], s = 20, color = color_1)
    q2 = ax[int(pd.np.floor(counter/4)),counter%4].scatter(df_refuerzo_dosis_aux['Fecha'], pd.np.add.accumulate(df_refuerzo_dosis_aux['Vacunas'].to_numpy())/df_cobertura['pob'].to_numpy()[grupo-1], s = 20, color = color_2)
    ax[int(pd.np.floor(counter/4)),counter%4].grid()
    ax[int(pd.np.floor(counter/4)),counter%4].set_title(num_grupo[grupo].replace('_','-'))
    plt.setp(ax[int(pd.np.floor(counter/4)),counter%4].get_xticklabels(), rotation=90)
    ax[int(pd.np.floor(counter/4)),counter%4].set_xlim(dt.datetime(2021,1,1),dt.datetime(2022,3,1))
    counter += 1
    
fig.legend([(q1),(q2)], ['FRED', 'Datos'], loc='lower center', bbox_to_anchor=(0.5,-0.03),
          fancybox=True, shadow=True, ncol=5, fontsize = 14)
    
fig.delaxes(ax[0,0])
fig.delaxes(ax[0,1])
fig.delaxes(ax[0,2])

plt.savefig('../../figures/cobertura_refuerzo_dosis.png', bbox_inches='tight', pad_inches=0.03)
