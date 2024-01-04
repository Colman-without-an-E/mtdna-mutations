#ifndef DEBUG_H
#define DEBUG_H

#include <stdio.h>
#include <math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "../include/parameters.h"

void print_state(int*** state, int* populations, int cell_idx);

void print_mutant_counts(int** mutant_counts);

void generate_state(const gsl_rng* rng, int*** state, int* populations, int pop_mu, double mut_p, int mut_id_mu);

void get_mutant_counts_wildtype(int** mutant_counts, int*** wildtype_state, int* wildtype_populations);

void inspect_true_mutant_counts_wildtype(int*** wildtype_state, int* wildtype_populations);

void get_mutant_counts(int** mutant_counts, int*** wildtype_state, int*** ra_state, int* wildtype_populations, int* ra_populations);

void inspect_true_mutant_counts(int*** wildtype_state, int*** ra_state, int* wildtype_populations, int* ra_populations);

#endif