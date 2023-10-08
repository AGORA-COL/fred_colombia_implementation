import pandas as pd
import numpy as np
import statsmodels.formula.api as sm
import seaborn as sns
import datetime as dtm
from datetime import datetime as dt
from importlib import reload
from google.cloud import bigquery
client = bigquery.Client()

today=pd.to_datetime(dt.today().date())

dicc_decenios={'01- 0 a 9 años':[0,9], 
               '02- 10 a 19 años':[10,19],
               '03- 20 a 29 años':[20,29], 
               '04- 30 a 39 años':[30,39], 
               '05- 40 a 49 años':[40,49], 
               '06- 50 a 59 años':[50,59],
               '07- 60 a 69 años':[60,69],
               '08- 70 a 79 años':[70,79], 
               '09- 80 años y más':[80,9999]}
dicc_etapas  ={'Etapa 5':[0,39], 
               'Etapa 4':[40,49], 
               'Etapa 3':[50,59],
               'Etapa 2':[60,79],
               'Etapa 1':[80,9999]}

def asosciar_grupo_etario(df,var_ref,dicc=dicc_decenios,nombre='GRUPO_ETARIO'):
    conditions=[(df[var_ref]>= dicc[ge][0]) & (df[var_ref]<= dicc[ge][1]) for ge in list(dicc.keys())]
    choices=[ge for ge in list(dicc.keys())]
    df[nombre]=np.select(conditions,choices,default=np.nan)
    return df

query="""
SELECT T_ID, N_IDENTIFICACION, EDAD, SEXO, FECHA_DE_VACUNACION_DDMMAAAA, DOSIS_APLICADA, LABORATORIO_BIOLOGICO
FROM `proyecto-covid-sds.vacunacion_cons_sis150_20_08_2021.SIS150*`
"""
df_vac_act=client.query(query).to_dataframe(progress_bar_type='tqdm')

df_vac_act['FECHA_DE_VACUNACION_DDMMAAAA']=pd.to_datetime(df_vac_act['FECHA_DE_VACUNACION_DDMMAAAA'])
df_vac_act['FECHA_DE_VACUNACION_DDMMAAAA']=df_vac_act['FECHA_DE_VACUNACION_DDMMAAAA'].dt.date
df_vac_act['EDAD']=pd.to_numeric(df_vac_act.EDAD,errors='coerce')
df_vac_act=asosciar_grupo_etario(df_vac_act,'EDAD')
df_vac_act=asosciar_grupo_etario(df_vac_act,'EDAD',dicc=dicc_etapas,nombre='ETAPA')

conditions=[df_vac_act['DOSIS_APLICADA'].str.contains('1'),df_vac_act['DOSIS_APLICADA'].str.contains('2'),df_vac_act['DOSIS_APLICADA'].str.contains('3'),df_vac_act['DOSIS_APLICADA'].str.contains('4')]
choices=['1 = Primera Dosis', '2 = Segunda Dosis', '3 = Única', '4 = Refuerzo']
df_vac_act['DOSIS_APLICADA']=np.select(conditions, choices, default='Sin información')

df_vac_act['T_ID']=df_vac_act.T_ID.str.upper().str.strip().str.replace('.','')
df_vac_act['LLAVE']=df_vac_act.T_ID+df_vac_act.N_IDENTIFICACION

df_vac_act=df_vac_act.drop_duplicates(subset=['T_ID','N_IDENTIFICACION','DOSIS_APLICADA'])

df_vac_t=df_vac_act[df_vac_act.DOSIS_APLICADA!='Sin información']\
                    .pivot(values=['FECHA_DE_VACUNACION_DDMMAAAA'],
                           index=['T_ID','N_IDENTIFICACION'],
                           columns='DOSIS_APLICADA').reset_index()
df_vac_t.columns = [' '.join(col).strip() for col in df_vac_t.columns]
df_vac_t=\
df_vac_t.merge(df_vac_act.drop_duplicates(subset=['LLAVE'])[['T_ID','N_IDENTIFICACION','EDAD','GRUPO_ETARIO','ETAPA']],
               how='left',on=['T_ID','N_IDENTIFICACION'])


df_vac_t['FECHA_COMPLETO_INMUNIZADO']=np.where(df_vac_t['FECHA_DE_VACUNACION_DDMMAAAA 3 = Única'].notnull(),
                                               df_vac_t['FECHA_DE_VACUNACION_DDMMAAAA 3 = Única']+dtm.timedelta(days=15),
                                               df_vac_t['FECHA_DE_VACUNACION_DDMMAAAA 2 = Segunda Dosis']+dtm.timedelta(days=15)
                                              )

df_vac_t['Vacunas'] = 1
print('Procesando informacion primera dosis')
df_aux_1 = df_vac_t.groupby(by =['FECHA_DE_VACUNACION_DDMMAAAA 1 = Primera Dosis', 'EDAD']).sum().reset_index()

print('Procesando informacion primera dosis')
df_aux_2 = df_vac_t.groupby(by =['FECHA_DE_VACUNACION_DDMMAAAA 2 = Segunda Dosis', 'EDAD']).sum().reset_index()

print('Procesando informacion primera dosis')
df_aux_3 = df_vac_t.groupby(by =['FECHA_DE_VACUNACION_DDMMAAAA 3 = Única', 'EDAD']).sum().reset_index()

print('Procesando informacion primera dosis')
df_aux_4 = df_vac_t.groupby(by =['FECHA_DE_VACUNACION_DDMMAAAA 4 = Refuerzo', 'EDAD']).sum().reset_index()

df_aux_1.columns = ['Fecha', 'Edad', 'Vacunas']
df_aux_2.columns = ['Fecha', 'Edad', 'Vacunas']
df_aux_3.columns = ['Fecha', 'Edad', 'Vacunas']
df_aux_4.columns = ['Fecha', 'Edad', 'Vacunas']

df_aux_1.to_csv('../data/primera_dosis.csv', index=False)
df_aux_2.to_csv('../data/segunda_dosis.csv', index=False)
df_aux_3.to_csv('../data/unica_dosis.csv', index=False)
df_aux_4.to_csv('../data/refuerzo_dosis.csv', index=False)
