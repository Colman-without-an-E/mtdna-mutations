# Data Analysis
This directory contains Python files for analysing simulation results. 

## Directory Structure
These files can be categorised into different types:
1. **Utility functions**
    - `utils.py` contains all the utility functions for data analysis, including compiling and ensembling simulations data, visualising data, helper functions on Pandas dataframes, and performing statistical tests. 
2. **Data processing**
    - `compile_data.py` compiles either populations data or site frequency spectrum data across simulations with the same parameter set. 
    ```bash
    python compile_data.py <parameters_set> <sim_type> <data_type>
    ```
    parameters_set is the index of the parameters set set in `include/parameters_.h`. sim_type is either `ra`, `ssd`, or `wildtype`. data_type is either `pop` for populations data, or `sfs` for site frequency spectrum data. 
    - `combine_pop_sfs.py` combines the already compiled populations data and site frequency spectrum data into one Pandas dataframe. Usage:
    ```bash
    python combine_pop_sfs.py <parameters_set> <sim_type>
    ```
    - `ensemble_data.py` retains only the simulations where a chosen threshold of mutant fraction is reached and discards the others, and then ensembles data across simulations. For the ensembled mutant fractions, the mean and standard error of the mean are computed at each recorded time and each cell. For the ensembled site frequency spectrum, they are simply the *unweighted ensemble site frequency spectrum* (see report). The number of retained simulations are also included for both data types, mutant fraction and site frequency spectrum. Usage: 
    ```bash
    python ensemble_data.py <parameters_set> [h_threshold]
    ```
    h_threshold is the threshold of mutant fraction for retaining simulations. 
3. **Examples**
    - `example_plot_h.py` plots the ensemble mutant fraction against time under the RA model
    - `example_plot_sfs.py` plots the ensemble site frequency spectrum recorded when mutant fraction threshold is first reached under the RA model
    - `example_mann-whitney.py` prints the Mann-Whitney U test results of the difference of site frequency spectra under the RA and SSD models of the same parameters set 1. 