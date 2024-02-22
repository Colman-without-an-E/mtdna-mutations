#ifndef PARAMETERS_H
#define PARAMETERS_H

// Save file locations
#define DIR_LOC "../data"

// Model parameters
#define CELLS 5
#define LOG_SITE_STD_MUTATION_RATE 4.0
#define DEGRADATION_RATE 7e-2
// #define DIFFUSION_RATE 0.1
#define DIFFUSION_RATE 1e-3
// #define DIFFUSION_RATE 0.0
// #define NUCLEUS_CONTROL_FACTOR 2e-4
#define NUCLEUS_CONTROL_FACTOR 2.5e-3
// #define TARGET_POP 3000
#define TARGET_POP 100
// #define REPLICATIVE_ADVANTAGE 4.8
#define REPLICATIVE_ADVANTAGE 1.0
// #define REPLICATIVE_ADVANTAGE 0.1
#define DENSITY 0.2
#define LEN_GENOME 16569

// Simulation parameters
#define SEED 10
#define N_SIMS 1000
#define BATCH_SIZE 100
#define SIM_LENGTH 10000.0
#define INTRODUCE_AT 3000.0
#define RECORDING_SPACE 100.0
#define N_BINS 20

// Break condition parameters
#define MAX_N_EVENTS 100000000
#define MAX_MUTANTS 100000

#endif