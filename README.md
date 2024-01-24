# FRED Agent Model for COVID-19 Outbreak

This repository contains the implementation of the FRED (A Framework for Reconstructing Epidemiological Dynamics) agent-based model for simulating the spread of the COVID-19 outbreak in various departments of Colombia. The model allows for a realistic representation of individual behaviors and interactions, allowing users to study the spread of the disease and assess the impact of various intervention strategies.

## Project Description

The FRED agent model is a stochastic, agent-based simulation package designed to study the spread of diseases through a population. The simulation follows individuals as they go about their daily routines and interact with each other, potentially spreading the disease. This implementation of the FRED model aims to analyze the COVID-19 pandemic situation across different departments in Colombia.

## Getting Started

### Dependencies

* C++ compiler
* R 4.0 or newer
* (Optional) Jupyter lab

### Installation

1. Clone this repository:

```bash
git clone https://github.com/AGORA-COL/fred_colombia_implementation.git
```

2. Navigate into the cloned repository:

```bash
cd fred_colombia_implementation
```

3. Run the FRED setup:
```bash
make setup
```
This will setup FRED from https://github.com/confunguido/FRED

4. Download the synthetic populations:
```bash
make download files="colombia_11001 colombia_1100101 ..."
```
You can select which populations to download and use.

## Input files

Here are listed the input files used in each simulation. Click at the link to access a visualization of the data.

[Dominance plot per department](https://dvelozad.github.io/fred_widgets/timeseries_plot.html)

[Beds usage plot per department](https://dvelozad.github.io/fred_widgets/bed_utilization_plot.html)

[Shelter in place input plot per department](https://dvelozad.github.io/fred_widgets/shelter_in_place_plot.html)


## Outputs

The program will output a series of data tables showing the progression of the disease, such as number of infected, number of recovered, number of deceased individuals, Hospitalizations, and ICU usage over time. 

This repo doesn't cointain any of FRED's outputs due to their large sizes. However here we include the input files neccesary to run the model.

## Contributing

We welcome contributions to this project. Please fork this repository and create a Pull Request with your changes.

## Contact

Please create an issue in this repository for any questions, suggestions, or comments. 

## Acknowledgments

* FRED is an agent-based modeling system developed by the Public Health Dynamics Laboratory at the University of Pittsburgh.