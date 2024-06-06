#ifndef SIM_H
#define SIM_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


/* Functions used in Gillespie algorithm */

void state_update_degrade(int*** state, int cell_idx, int degrade_row, int* populations, int** mutant_counts);

void state_update_replicate(int*** state, int cell_idx, int replicate_row, int* populations, int** mutant_counts, int n_new_std_mutations);

void state_update_diffuse(int*** state, int cell_idx_from, int cell_idx_to, int diffuse_row, int* populations, int** mutant_counts);

int choose_neighbouring_cell(const gsl_rng* rng, int cell_idx);

void wildtype_propensity_update(double** propensity, double* propensity_sums, int cell_idx, int* wildtype_populations, double degradation_rate, double diffusion_rate, double nucleus_control_factor, double target_population);

void propensity_update(double** propensity, double* propensity_sums, int cell_idx, int* wildtype_populations, int* ra_or_ssd_populations, double degradation_rate, double diffusion_rate, double nucleus_control_factor, double target_population, double replicative_advantage);

int weighted_sample(const gsl_rng* rng, int n, double* weights);

void wildtype_gillespie_event(const gsl_rng* rng, double** propensity, double* propensity_sums, int*** wildtype_state, int* wildtype_populations, int** mutant_counts, long double site_std_mutation_rate, double degradation_rate, double diffusion_rate, double nucleus_control_factor, double target_population);

void gillespie_event(const gsl_rng* rng, double** propensity, double* propensity_sums, int*** wildtype_state, int*** ra_or_ssd_state, int** mutant_counts, int* wildtype_populations, int* ra_or_ssd_populations, long double site_std_mutation_rate, double degradation_rate, double diffusion_rate, double nucleus_control_factor, double target_population, double replicative_advantage);

void get_max_std_heteroplasmies(double* max_std_heteroplasmies, int** mutant_counts, int* wildtype_populations);

double get_max_ra_or_ssd_heteroplasmy(int* wildtype_populations, int* ra_or_ssd_populations);

/* End of Gillespie algorithm functions */


/* Data storage functions */

void introduce_ra_or_ssd(const gsl_rng* rng, int*** wildtype_state, int*** dest_state, int* wildtype_populations, int* dest_populations, int cell_idx);

void compact_relabel_wildtype_mutations(int** mutant_counts, int*** wildtype_state, int* wildtype_populations);

void compact_relabel_mutations(int** mutant_counts, int*** wildtype_state, int*** ra_or_ssd_state, int* wildtype_populations, int* ra_or_ssd_populations);

void copy_state(int*** source, int*** dest, int* nrows);

void copy_mutant_counts(int** source, int** dest);

void copy_population(int* source, int* dest);

void write_data_to_file_pre_introduce(int* wildtype_populations, int** mutant_counts, double recording_time, char* population_filename, char* sfs_filename);

void write_data_to_file(int* wildtype_populations, int* ra_or_ssd_populations, int** mutant_counts, double recording_time, char* population_filename, char* sfs_filename);

/* End of data storage */


/* Debug functions */

int argmax1D(double* array, int length);

void print_state(int*** state, int* populations, int cell_idx);

void print_mutant_counts(int** mutant_counts);

void generate_state(const gsl_rng* rng, int*** state, int* populations, int pop_mu, double mut_p, int mut_id_mu);

void get_mutant_counts_wildtype(int** mutant_counts, int*** wildtype_state, int* wildtype_population);

void get_mutant_counts(int** mutant_counts, int*** wildtype_state, int*** ra_state, int* wildtype_population, int* ra_population);

int identical_mutant_counts(int** mutant_counts1, int** mutant_counts2);

int argmax(double* array, int length);

int argmax_mutant_counts(int* args, int** mutant_counts);

/* End of debug functions */


#endif