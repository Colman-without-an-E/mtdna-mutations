#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "../include/parameters_hpc.h"
#include "../include/lib_sim.h"

#define RND gsl_rng_uniform_pos(rng) // generate number from Unif(0,1)

// Define propensity updating rule
void propensity_update(double** propensity, double* propensity_sums, int cell_idx, int* wildtype_populations, int* ssd_populations, double degradation_rate, double diffusion_rate, double nucleus_control_factor, double target_population, double density) {
	return;
}


int main(int argc, char *argv[]) {
    int seed;
    const double log_site_std_mutation_rate = LOG_SITE_STD_MUTATION_RATE;
    const double degradation_rate = DEGRADATION_RATE;
	const double diffusion_rate = DIFFUSION_RATE;
    const double nucleus_control_factor = NUCLEUS_CONTROL_FACTOR;
    const double target_population = TARGET_POP;
	const double density = DENSITY;

	if (argc==2) {
		seed = atoi(argv[1]);
	} else {
		printf("argc = %d\n", argc);
		printf("Usage: %s  seed\n", argv[0]);
		return 0;
	}

	/* set up GSL RNG */
	gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(rng, seed);
	/* end of GSL setup */

	// Generate replication rate and mutation rate for this run
    const long double site_std_mutation_rate = pow(10, - log_site_std_mutation_rate);

    // Change to directory which stores simulation results. If directory doesn't exist, create one. 
    if (chdir(DIR_NAME)) {
		mkdir(DIR_NAME, 0700); // Linux
		// mkdir(DIR_NAME); // Windows
		chdir(DIR_NAME);
    }

	// Set up file to save parameter values
	FILE *fp_parameters = fopen("parameters.txt", "w");
	fprintf(fp_parameters, "parameter,value\n");
	fprintf(fp_parameters, "cells,%d\n", CELLS);
	fprintf(fp_parameters, "log_site_std_mutation_rate,%e\n", LOG_SITE_STD_MUTATION_RATE);
	fprintf(fp_parameters, "degradation_rate,%e\n", DEGRADATION_RATE);
	fprintf(fp_parameters, "diffusion_rate,%e\n", DIFFUSION_RATE);
	fprintf(fp_parameters, "nucleus_control_factor,%e\n", NUCLEUS_CONTROL_FACTOR);
	fprintf(fp_parameters, "target_pop,%e\n", TARGET_POP);
	fprintf(fp_parameters, "density,%e\n", DENSITY);
	fprintf(fp_parameters, "len_genome,%d\n", LEN_GENOME);

	fprintf(fp_parameters, "sim_length,%e\n", SIM_LENGTH);
	fprintf(fp_parameters, "introduce_after,%e\n", INTRODUCE_AFTER);
	fprintf(fp_parameters, "recording_space,%e\n", RECORDING_SPACE);

	fprintf(fp_parameters, "max_n_events,%d\n", MAX_N_EVENTS);
	fprintf(fp_parameters, "max_mutants,%d\n", MAX_MUTANTS);
	fclose(fp_parameters);

	// Set up files to write population data in
	char wildtype_population_filename[50];
	sprintf(wildtype_population_filename, "wildtype_sim_populations%d.txt", seed);
	FILE *fp_wildtype_population = fopen(wildtype_population_filename, "w");
	fprintf(fp_wildtype_population, "t,cell,w,m\n");
	fclose(fp_wildtype_population);

	// Set up files to write site frequency spectrum data in
	char wildtype_sfs_filename[50];
	sprintf(wildtype_sfs_filename, "wildtype_sim_site_frequency_spectrum%d.txt", seed);
	FILE *fp_wildtype_sfs = fopen(wildtype_sfs_filename, "w");
	fprintf(fp_wildtype_sfs, "t,cell,sfs\n");
	fclose(fp_wildtype_sfs);

	// Declare variables to store standard mutational information of wildtype individuals
	/* wildtype_state contains mutational information about wildtype individuals
	wildtype_state[k][i][0] is the number of standard mutations which individual i possesses in cell k
	wildtype_state[k][i][1:] are the identities of the standard mutations which individual i possesses in cell k */
	int*** wildtype_state = malloc(CELLS * sizeof(int**));
	/* wildtype_populations[k] is the number of wildtype individuals in cell k */
	int* wildtype_populations = malloc(CELLS * sizeof(int));

	/* mutant_counts[0][0] is the latest mutation identity across the whole system
	mutant_counts[1:][0] = 0 are dummy entries
	mutant_counts[k][i] is the number of individuals in cell k with mutation identity i
	mutant_counts has number of columns mutant_counts[0][0] + 1 */
	int** mutant_counts = malloc(CELLS * sizeof(int*));

	/* wildtype_propensity[k][0] is the propensity of degradation of a wildtype individual in cell k
	wildtype_propensity[k][1] is the propensity of replication of a wildtype individual in cell k
	wildtype_propensity[k][2] is the propensity of diffusion of a  wildtype in cell k */
	double** wildtype_propensity = malloc(CELLS * sizeof(double*));
	/* wildtype_propensity_sums[k] is the sum of the wildtype propensity vector of cell k */
	double* wildtype_propensity_sums = malloc(CELLS * sizeof(double));

	for (int k=0; k<CELLS; ++k) {
		// Initialise wildtype population steady state
		wildtype_state[k] = malloc((int) target_population * sizeof(int*));
		for (int i=0; i<(int) target_population; ++i) {wildtype_state[k][i] = calloc(1, sizeof(int));}
		wildtype_populations[k] = (int) target_population;
	
		// Initialise mutant counts
		mutant_counts[k] = calloc(1, sizeof(int));

		// Compute initial wildtype propensity
		wildtype_propensity[k] = malloc(3 * sizeof(double));
		wildtype_propensity_update(wildtype_propensity, wildtype_propensity_sums, k, wildtype_populations, degradation_rate, diffusion_rate, nucleus_control_factor, target_population);
	}

	int total_wildtype_population;
	double wildtype_propensity_sum_across_cells;
	double current_time = 0.0;
	double recording_time = 0.0;
	int n_event = 0;
	
	// Record initial data
	write_data_to_file_pre_introduce(wildtype_populations, mutant_counts, recording_time, wildtype_population_filename, wildtype_sfs_filename);
	recording_time += RECORDING_SPACE;

	// Gillespie algorithm until time threshold reached
	while (current_time<SIM_LENGTH && n_event < MAX_N_EVENTS) {

		// Realise event according to propensity
		wildtype_propensity_sum_across_cells = 0;
		for (int k=0; k<CELLS; ++k) {wildtype_propensity_sum_across_cells += wildtype_propensity_sums[k];}
		current_time += -log(RND) / wildtype_propensity_sum_across_cells;
		wildtype_gillespie_event(rng, wildtype_propensity, wildtype_propensity_sums, wildtype_state, wildtype_populations, mutant_counts, site_std_mutation_rate, degradation_rate, diffusion_rate, nucleus_control_factor, target_population);
		n_event++;

		// Record data and end simulation if extinction
		total_wildtype_population = 0;
		for (int k=0; k<CELLS; ++k) {total_wildtype_population += wildtype_populations[k];}
		if (total_wildtype_population==0) {
			compact_relabel_wildtype_mutations(mutant_counts, wildtype_state, wildtype_populations);
			write_data_to_file_pre_introduce(wildtype_populations, mutant_counts, recording_time, wildtype_population_filename, wildtype_sfs_filename);
			break;
		}

		// Record data at recording time
		if (current_time>=recording_time) {
			compact_relabel_wildtype_mutations(mutant_counts, wildtype_state, wildtype_populations);
			write_data_to_file_pre_introduce(wildtype_populations, mutant_counts, recording_time, wildtype_population_filename, wildtype_sfs_filename);

			recording_time += RECORDING_SPACE;
		}
	}

	// Free memory
	for (int k=0; k<CELLS; ++k) {
		for (int i=0; i<wildtype_populations[k]; ++i){free(wildtype_state[k][i]);}
		free(wildtype_state[k]);
		free(wildtype_propensity[k]);
		free(mutant_counts[k]);
		
	}
	free(wildtype_state);
	free(wildtype_populations);

	free(wildtype_propensity);
	free(wildtype_propensity_sums);

	free(mutant_counts);

    gsl_rng_free(rng);
	return 0;
}