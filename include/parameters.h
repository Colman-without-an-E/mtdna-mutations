#ifndef PARAMETERS_H
#define PARAMETERS_H

// Model parameters
#define CELLS 2
#define LOG_SITE_STD_MUTATION_RATE 4.0
#define DEGRADATION_RATE 1.0
#define DIFFUSION_RATE 0.1
#define NUCLEUS_CONTROL_FACTOR 2.5e-3
#define TARGET_POP 20
#define REPLICATIVE_ADVANTAGE 0.1
#define TARGET_RA_H 0.25
#define LEN_GENOME 16569

// Simulation parameters
#define SEED 1
#define N_SIMS 2
#define BATCH_SIZE 1
#define SIM_LENGTH 100.00
#define RECORDING_SPACE 10.0
#define N_BINS 20

// Break condition parameters
#define MAX_N_EVENTS 100000
#define MAX_MUTANTS 100000

// Save file locations
#define RA_POP "../data/ra_population_sims_%d.txt"
#define RA_SFS "../data/ra_site_frequency_spectrum_sims_%d.txt"
#define SSD_POP "../data/ssd_population_sims_%d.txt"
#define SSD_SFS "../data/ssd_site_frequency_spectrum_sims_%d.txt"

#endif