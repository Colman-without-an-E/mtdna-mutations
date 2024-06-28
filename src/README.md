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
- `wildtype_state`: A 3-dimensional array storing the point mutational information of each wildtype individual of the system at current time. `wildtype_state[k][i][0]` is the number of standard mutations which individual `i` possesses in cell `k`. `wildtype_state[k][i][1:]` are the identities of the standard mutations which individual `i` possesses in cell `k`. `wildtype_state` is dynamically allocated such that the size of the array is exactly the number of individuals in a cell, and also the 

## Interpreting simulation results