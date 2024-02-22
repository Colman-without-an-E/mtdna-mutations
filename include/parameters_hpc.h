#ifndef PARAMETERS_H
#define PARAMETERS_H

// Save file locations
#define DIR_LOC "../data/parameters_set1"

// Model parameters
#define CELLS 100
#define LOG_SITE_STD_MUTATION_RATE 4.0
#define DEGRADATION_RATE 7e-2
#define DIFFUSION_RATE 7e-2
#define NUCLEUS_CONTROL_FACTOR 2.5e-3
#define TARGET_POP 20
#define REPLICATIVE_ADVANTAGE 0.2
#define DENSITY 0.2
#define RATE_DIFFERENCE 7e-2
#define LEN_GENOME 16569

// Simulation parameters
#define SIM_LENGTH 10000.0
#define INTRODUCE_AFTER 1000.0
#define RECORDING_SPACE 100.0
#define N_BINS 50

// Break condition parameters
#define MAX_N_EVENTS 100000000
#define MAX_MUTANTS 100000

#endif