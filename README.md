#  _FRED Agent-Based Model for COVID-19 – Implementación para Colombia _ 

Este repositorio contiene una implementación adaptada del modelo FRED (A Framework for Reconstructing Epidemiological Dynamics) para simular la dinámica del brote de COVID-19 en diferentes departamentos de Colombia. Este modelo permite representar de manera realista las interacciones y comportamientos individuales en la población, y evaluar el efecto de diversas intervenciones farmacéuticas y no farmacéuticas.

## Descripción
El modelo FRED es un sistema de simulación estocástico basado en agentes, desarrollado por el Public Health Dynamics Laboratory de la Universidad de Pittsburgh. Esta adaptación busca apoyar el análisis de la situación epidemiológica en Colombia durante la pandemia de COVID-19 y evaluar el impacto de distintas medidas de control a nivel local y nacional.

La presente implementación hace parte del producto:
“Modelo de aproximación estadística: soporte en Excel o código de programación que justifica los cálculos planteados en los diferentes modelos utilizados con el fin de simular y evaluar situación epidemiológica, y/o medidas farmacéuticas y/o no farmacéuticas en relación con la pandemia COVID-19 en Colombia.”

## Primeros pasos
### Requisitos
- Compilador de C++
- R (versión 4.0 o superior)
- (Opcional) Jupyter Lab para visualización adicional

### Instalación

1. Clona este repositorio:

```bash
git clone https://github.com/AGORA-COL/fred_colombia_implementation.git
```

2. Navegue hasta el repositorio clonado:

```bash
cd fred_colombia_implementation
```


3. Ejecuta el script de instalación de FRED:
```bash
make setup
```
Esto descargará y configurará la versión base de FRED desde
https://github.com/confunguido/FRED

4. Descarga las poblaciones sintéticas necesarias:
```bash
make download files="colombia_11001 colombia_1100101 ..."
```
Puedes modificar los códigos para seleccionar los municipios o departamentos que desees simular.

## Archivos de entrada (Input files)
Este repositorio incluye los insumos necesarios para correr las simulaciones en FRED:

- Dominancia de variantes por departamento: [Dominance plot per department](https://dvelozad.github.io/fred_widgets/dominance_widget.html)
- Ocupación de camas hospitalarias: [Beds usage plot per department](https://dvelozad.github.io/fred_widgets/bed_utilization_plot.html)
- Parámetros de movilidad y confinamiento: [Shelter and mobility input plot per department](https://dvelozad.github.io/fred_widgets/mobility_shelter_trends.html)
- Contactos por tipo de espacio: [Contacts input plot per department](https://dvelozad.github.io/fred_widgets/mobility_contacts_trends.html)


## Salidas del modelo (Outputs files)
El programa generará una serie de tablas de datos que muestran la progresión de la enfermedad, como el número de infectados, el número de recuperados, el número de fallecidos, las hospitalizaciones y el uso de la UCI a lo largo del tiempo.

⚠️ Nota: Este repositorio no incluye los archivos de salida debido a su gran tamaño. Sin embargo, se proveen todos los insumos necesarios para ejecutar el modelo y reproducir los análisis.

## Contribuciones
Las contribuciones son bienvenidas. Puedes:

- Hacer un fork del repositorio
- Crear una rama con tus cambios
- Enviar un Pull Request

## Autores
- Diego Veloza Díaz.
- Guido Camargo España.
- Daniel Santiago Bonilla Betancourth.
- Jennifer Murillo-Alvarado.
- Andrés Moreno.
- Jaime Pavlich-Mariscal.
- Zulma M. Cucunubá.


## Financiación
Esta investigación fue financiada por el Ministerio de Ciencia, Tecnología e Innovación de Colombia, proyecto ÁGORA: “Alianza para la Generación de Evidencia sobre COVID-19, su Respuesta y Lecciones Aprendidas para la Postpandemia y Futuras Epidemias” (Contrato N° 637-2022).

## Cómo citar este recurso
Si utilizas esta rutina o sus sistemas de agrupación en tus análisis o publicaciones, por favor citarlo de la siguiente manera:

Veloza Díaz, D., Camargo España, G., Bonilla Betancourth, D. S., Murillo-Alvarado, J., Moreno, A., Pavlich-Mariscal, J., & Cucunubá, Z. M. (2025). ÁGORA: Ciencia y decisiones en salud pública. Proyecto ÁGORA Colombia.
Disponible en: https://github.com/AGORA-COL/fred_colombia_implementation

También puedes exportar esta cita en formatos como BibTeX, RIS, APA y más desde el botón “Cite this repository” en la parte superior derecha de este repositorio (disponible si has agregado el archivo CITATION.cff).

## Contacto
Si tienes preguntas, sugerencias o comentarios, por favor crea un Issue en este repositorio.

## Agradecimientos
* FRED es un sistema de modelado basado en agentes desarrollado por el Laboratorio de Dinámica de Salud Pública de la Universidad de Pittsburgh.

  
