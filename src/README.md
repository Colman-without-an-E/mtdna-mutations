# SRC

This directory contains the C source files for simulating mtDNA mutations under different models. 

## Files

- `_debug.c` contains debugging and printing functions. Can be compiled together with other source files to create an executable, which can then be run using `gdb` for debugging purposes. 
- `lib_sim.c` provides functions for the general simulation framework. 
    - It implements a structure called Dict, which is somewhat equivalent to a dict type in Python. This is used to store site frequency spectrum information. 
    - It contains functions that update the state of the system based on which of the 3 events occur: replication, degradation, or diffusion. Occurence of these events are determined probabilistically according to the propensity of the state. Note that only the wildtype propensity update rule `wildtype_propensity_update` is given because the mutant propensity update rule is determined in the `*single_sim.c` source file under the specified model. 
    - It also contains functions that write data to text files. 
- `ra_single_sim.c` is the source file for simulating a system under the RA model. 
- `ssd_single_sim.c` is the source file for simulating a system under the SSD model. 
- `sss_single_sim.c` is the source file for simulating a system under the SSS model. 
- `wildtype_single_sim.c` is the source file for simulating a system without any mutants. 

## Usage
GCC is required. Without debugging, the `*single_sim.c` source file should be compiled together with `lib_sim.c`:
```bash
gcc [ra, ssd, sss, or wildtype]_single_sim.c lib_sim.c -o single_sim.exe -lgsl -lgslcblas -lm
./single_sim.exe 1
```
where an argument for the seed is required for the executable. 

## Important variable names and meaning
Most important variables are explained in the C file. However, we go through them here as well:
- `wildtype_state`: A 3-dimensional array storing the point mutational information of each wildtype individual of the system at current time. `wildtype_state[k][i][0]` is the number of point mutations which individual `i` possesses in cell `k`. `wildtype_state[k][i][1:]` are the identities of the point mutations which individual `i` possesses in cell `k`. `wildtype_state` is dynamically allocated such that the size of the array is exactly the number of individuals in a cell, and adjusts to the latest point mutation identity. 
- `ra_state`, `ssd_state`, `sss_state`: Same as above, but with deletion mutants under the specified model instead. Point mutational information of deletion mutants and wildtypes (individuals with no deletion mutations) are stored separately. 
- `mutant_counts`: A 2-dimensional array storing the number of individuals of each point mutation. Each newly arose point mutation is assigned a new identity (positive integer incrementing by one), and the `mutant_counts[0][0]` is the latest mutation identity across the whole system. `mutant_counts[1:][0]` = 0 are dummy entries. `mutant_counts[k][i]` is the number of individuals in cell `k` with mutation identity `i`. `mutant_counts` is dynamically allocated to have number of columns `mutant_counts[0][0]` + 1. 
- `wildtype_propensity`: A 2-dimensional array of wildtype propensity vector. `wildtype_propensity[k][0]` is the propensity of degradation of a wildtype individual in cell `k`. `wildtype_propensity[k][1]` is the propensity of replication of a wildtype individual in cell `k`. `wildtype_propensity[k][2]` is the propensity of diffusion of a  wildtype in cell `k`.
- `propensity`: Same as above, but the propensity of deletion mutants instead. 

## Interpreting simulation results
Running the executable will generate 3 text files in the specified directory (will make a new one if it does not exist) in `include/parameters.h`. The first text file generated is `parameters.txt`, which is just a list of parameters and values for that simulation. The actual simulation results come in two files:
1. `[sim_type]_sim_populations[seed].txt`
2. `[sim_type]_sim_site_frequency_spectrum[seed].txt`
where `[sim_type]`, either `ra`, `ssd`, or `sss` specifies the model, and `[seed]` is the seed input as the argument of the executable. 

### Populations data
`[sim_type]_sim_populations[seed].txt` contains 4 columns of headers `t`, `cell`, `w`, and `m`. `t` is time, the interval of which is determined by the `RECORDING_SPACE` macro in `include/parameters.h`. At each `t`, the population size of wildtypes and deletion mutants, `w` and `m` respectively, are recorded at each cell, indexed from 0 by the header `cell`. Notice that the `cell`=-1 encodes the population sizes across all cells. When `[sim_type]` is `wildtype`, column `m` is always 0. 

### Site frequency spectrum data
`[sim_type]_sim_site_frequency_spectrum[seed].txt` contains 3 columns of headers `t`, `cell`, and `sfs`. `t` is time and `cell` is the cell index as above. `sfs` is not actually the site frequency spectrum. It is a counter of point mutations based on their prevalence in the entire system. This can be explained more intuitively via an example. Consider `{"1":2,"4":1}`. This means that there are 2 unique point mutations (recall each point mutation is given its own unique identity when it first arises) that are present in exactly 1 individual in the system. Note that it does not matter if they are wildtypes or deletion mutants, and that the 2 point mutations need not be present in the same individual. In addition to that, there is 1 unique point mutation that is present in exactly 4 individuals. One can think of `{"1":2,"4":1}` as a shorthand for `{1, 1, 4}`, the distribution of the number of individuals containing each unique point mutation. To obtain the site frequency spectrum, it suffices to divide the key, that is the integer in double quotations, by the total number of inviduals. This will give us the distribution of heteroplasmy of each unique point mutation, precisely the definition of site frequency spectrum. This step is performed in `data_analysis/ensemble_data.py`. 