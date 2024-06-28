#ifndef PARAMETERS_H
#define PARAMETERS_H

// Save file locations
#define DIR_NAME "parameters_set100"

// Model parameters
#define CELLS 100 // set 1
// #define CELLS 2 // set 2,3,4,5
#define LOG_SITE_STD_MUTATION_RATE 5.64345
#define DEGRADATION_RATE 7e-2
#define DIFFUSION_RATE 7e-2 // set 1,2,5,7,11,12,13,14,15,21
// #define DIFFUSION_RATE 7e-4 // set 3,6,22
// #define DIFFUSION_RATE 7.0 // set 4,23
#define NUCLEUS_CONTROL_FACTOR 2.5e-3
// #define TARGET_POP 5 // set 1
// #define TARGET_POP 20
// #define TARGET_POP 100 // set 21,22,23
// #define TARGET_POP 96.6 // set set 11
// #define TARGET_POP 96.8 // set set 12
// #define TARGET_POP 97.0 // set set 13
// #define TARGET_POP 97.2 // set set 14
// #define TARGET_POP 97.4 // set set 15
#define TARGET_POP 13.50896 // set set 100
// #define TARGET_POP 15.424136 // set set 101
// #define TARGET_POP 17.339308 // set set 102
// #define TARGET_POP 19.25448 // set set 103
// #define TARGET_POP 21.169654 // set set 104
// #define REPLICATIVE ADVANTAGE 0.05 // set 1
// #define REPLICATIVE_ADVANTAGE 0.2 // set 2,3,4,6
// #define REPLICATIVE_ADVANTAGE 0.005374272690266158// set 5,7
// #define REPLICATIVE_ADVANTAGE 8.5e-3 // set 11
// #define REPLICATIVE_ADVANTAGE 8.0e-3 // set 12
// #define REPLICATIVE_ADVANTAGE 7.5e-3 // set 13
// #define REPLICATIVE_ADVANTAGE 7.0e-3 // set 14
// #define REPLICATIVE_ADVANTAGE 6.5e-3 // set 15
#define REPLICATIVE_ADVANTAGE 2.872760e-2 // set 100
// #define REPLICATIVE_ADVANTAGE 2.393966e-2 // set 101
// #define REPLICATIVE_ADVANTAGE 1.915173e-2 // set 102 (wave speed formula)
// #define REPLICATIVE_ADVANTAGE 1.436380e-2 // set 103
// #define REPLICATIVE_ADVANTAGE 9.575865e-3 // set 104
// #define REPLICATIVE_ADVANTAGE 0.0 // set 20,21,22,23
#define DENSITY 0.2
// #define DENSITY 1.0 // set 20
#define RATE_DIFFERENCE 7e-2
#define LEN_GENOME 16569

// Simulation parameters
#define SIM_LENGTH 20000.0 // set 1,6,7
// #define SIM_LENGTH 5000.0 // set 2,3,4,5
#define INTRODUCE_AFTER 10000.0 // set 1,6,7
// #define INTRODUCE_AFTER 2500.0 // set 2,3,4,5
#define RECORDING_SPACE 100.0 // set 1,6,7
// #define RECORDING_SPACE 10.0 // set 2,3,4,5

// Break condition parameters
#define MAX_N_EVENTS 10000000000
#define MAX_MUTANTS 1000000

#endif
