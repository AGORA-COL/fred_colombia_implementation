from google.cloud import bigquery
import pandas as pd
import google.auth

# Create credentials with Drive & BigQuery API scopes.
# Both APIs must be enabled for your project before running this code.
#
# If you are using credentials from gcloud, you must authorize the
# application first with the following command:
#
# gcloud auth application-default login --scopes=https://www.googleapis.com/auth/drive,https://www.googleapis.com/auth/cloud-platform
credentials, project = google.auth.default(
    scopes=[
        "https://www.googleapis.com/auth/drive",
        "https://www.googleapis.com/auth/cloud-platform",
    ]
)

# Construct a BigQuery client object.
client = bigquery.Client(credentials=credentials, project=project)

def FRED_dominance():
    query = """
    SELECT * FROM `asis-sds.FRED.FRED_dominance`   
    """
    df = client.query(query).to_dataframe(progress_bar_type ='tqdm')
    return df

def FRED_cases():
    query = """
    SELECT * FROM `asis-sds.FRED.FRED_cases`   
    """
    df = client.query(query).to_dataframe(progress_bar_type ='tqdm')
    return df

def FRED_deaths():
    query = """
    SELECT * FROM `asis-sds.FRED.FRED_deaths`   
    """
    df = client.query(query).to_dataframe(progress_bar_type ='tqdm')
    return df

def FRED_solicitudes_UCI():
    query = """
    SELECT SAFE_CAST(FECHA_SOLICITUD AS DATE) AS FECHA_SOLICITUD, ESTADO, COUNT(*) AS SOLICITUDES
    FROM `asis-sds.SECRETARIA_SALUD.crue_uci_solicitudes_asis` 
    WHERE TIPO_UCI = 'COVID-19'
    GROUP BY FECHA_SOLICITUD, ESTADO 
    ORDER BY FECHA_SOLICITUD, ESTADO   
    """
    df = client.query(query).to_dataframe(progress_bar_type ='tqdm')
    return df

def FRED_solicitudes_UCI_edad():
    query = """
    SELECT SAFE_CAST(FECHA_SOLICITUD AS DATE) AS FECHA_SOLICITUD, ESTADO, COUNT(*) AS SOLICITUDES, EDAD
    FROM `asis-sds.SECRETARIA_SALUD.crue_uci_solicitudes_asis` 
    WHERE TIPO_UCI = 'COVID-19'
    GROUP BY FECHA_SOLICITUD, ESTADO, EDAD
    ORDER BY FECHA_SOLICITUD, ESTADO, EDAD
    """
    df = client.query(query).to_dataframe(progress_bar_type ='tqdm')
    return df

def FRED_vacunaciones():
    """
    OJO: Revisar este query
    """
    query = """
    SELECT 
      FechaAplicacion AS Date, 
      'Bogotá' AS Department,
      COUNT(*) AS VaccinesApplied 
    FROM `asis-sds.SECRETARIA_SALUD.vacunacion_PAI_asis` 
    WHERE FechaAplicacion IS NOT NULL AND FechaAplicacion BETWEEN '2021-01-01' AND CURRENT_DATE() - 14
    GROUP BY FechaAplicacion ORDER BY FechaAplicacion
    """
    df = client.query(query).to_dataframe(progress_bar_type ='tqdm')
    return df

def FRED_muertes():
    query = """
    WITH 
      CASOS AS (
        SELECT FECHA_DE_INGRESO AS Date, COUNT(*) AS Cases FROM `asis-sds.SECRETARIA_SALUD.positivos_asis` GROUP BY Date
      ),
      FALLECIDOS AS(
        SELECT FECHA_DE_MUERTE AS Date, COUNT(*) AS Deaths FROM `asis-sds.SECRETARIA_SALUD.positivos_asis` WHERE RECUPERADO = 'Fallecido' GROUP BY Date
      )

    SELECT 
      11001 AS MunCode,
      * EXCEPT(Deaths),
      (CASE WHEN Deaths IS NULL THEN 0
       ELSE Deaths END) AS Deaths
    FROM(
      SELECT A.*, B.Deaths FROM CASOS AS A
      LEFT JOIN
      (SELECT * FROM FALLECIDOS) AS B
      ON A.Date = B.Date
    ) ORDER BY Date
    """
    df = client.query(query).to_dataframe(progress_bar_type ='tqdm')
    return df

def FRED_muertes_edad():
    query = """
    SELECT 
      11001 AS MunCode,
      FECHA_DE_MUERTE AS Date, 
      (CASE WHEN EDAD_DECENIOS_ASIS = '01- 0 a 9 años' THEN 'ACF0_10' 
            WHEN EDAD_DECENIOS_ASIS = '02- 10 a 19 años' THEN 'ACF10_20' 
            WHEN EDAD_DECENIOS_ASIS = '03- 20 a 29 años' THEN 'ACF20_30' 
            WHEN EDAD_DECENIOS_ASIS = '04- 30 a 39 años' THEN 'ACF30_40' 
            WHEN EDAD_DECENIOS_ASIS = '05- 40 a 49 años' THEN 'ACF40_50' 
            WHEN EDAD_DECENIOS_ASIS = '06- 50 a 59 años' THEN 'ACF50_60' 
            WHEN EDAD_DECENIOS_ASIS = '07- 60 a 69 años' THEN 'ACF60_70' 
            WHEN EDAD_DECENIOS_ASIS = '08- 70 a 79 años' THEN 'ACF70_80' 
       ELSE 'ACF80_120' END) AS AgeGroup, 
      COUNT(*) AS Deaths 
    FROM `asis-sds.SECRETARIA_SALUD.positivos_asis` 
    WHERE RECUPERADO = 'Fallecido'
    GROUP BY FECHA_DE_MUERTE, EDAD_DECENIOS_ASIS
    ORDER BY Date, AgeGroup
    """
    df = client.query(query).to_dataframe(progress_bar_type ='tqdm')
    return df