#ifndef PARAMETERS_H
#define PARAMETERS_H

// Save file locations
#define DIR_NAME "parameters_set1"

// Model parameters
#define CELLS 2 // number of cells
#define LOG_SITE_STD_MUTATION_RATE 5.64345 // 10^(-LOG_SITE_STD_MUTATION_RATE) is the true mutation rate
#define DEGRADATION_RATE 7e-2 // degradation rate per capita. mu in literature.
#define DIFFUSION_RATE 7e-2 // diffusion rate per capita. gamma in literature.
#define NUCLEUS_CONTROL_FACTOR 2.5e-3 // nucleus control factor. c in literature.
#define TARGET_POP 20 // target population. Nss in literature.
#define REPLICATIVE_ADVANTAGE 0.2 // replicate advantage coefficient in the RA model. kra in literature.
#define DENSITY 0.2 // inverse density in the SSD model. delta in literature
#define RATE_DIFFERENCE 7e-2 // rate difference in the SSS model. epsilon in literature. 
#define LEN_GENOME 16569 // length of mtDNA genome

// Simulation parameters
#define SIM_LENGTH 5000.0 // total simulation length in time
#define INTRODUCE_AFTER 2500.0 // introduce a single deletion mutant at this time
#define RECORDING_SPACE 10.0 // time interval between each data record

// Break condition parameters
#define MAX_N_EVENTS 10000000000 // maximum number of events allowed before breaking
#define MAX_MUTANTS 1000000 // maximum number of mutants allowed before breaking

#endif
