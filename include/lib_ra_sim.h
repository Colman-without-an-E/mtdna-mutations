#ifndef RA_SIM_H
#define RA_SIM_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <limits.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

void propensity_update(double** propensity, double* propensity_sums, int cell_idx, int* wildtype_populations, int* ra_or_ssd_populations, double degradation_rate, double diffusion_rate, double nucleus_control_factor, int target_population, double replicative_advantage);

void gillespie_event(const gsl_rng* rng, double** propensity, double* propensity_sums, int*** wildtype_state, int*** ra_state, int** mutant_counts, int* wildtype_populations, int* ra_populations, long double site_std_mutation_rate, double degradation_rate, double diffusion_rate, double nucleus_control_factor, int target_population, double replicative_advantage);

void compact_relabel_wildtype_mutations(int** mutant_counts, int*** wildtype_state, int* wildtype_populations);

void compact_relabel_mutations(int** mutant_counts, int*** wildtype_state, int*** ra_state, int* wildtype_populations, int* ra_populations);

void copy_state(int*** source, int*** dest, int* nrows);

void copy_mutant_counts(int** source, int** dest);

void copy_population(int* source, int* dest);

#endif