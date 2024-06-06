#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "../include/parameters_hpc.h"
#include "../include/lib_sss_sim.h"

#define RND gsl_rng_uniform_pos(rng) // generate number from Unif(0,1)

// Define propensity updating rule for SSS simulation
void propensity_update(double** propensity, double* propensity_sums, int cell_idx, int* wildtype_populations, int* sss_populations, double degradation_rate, double diffusion_rate, double nucleus_control_factor, int target_population, double rate_difference) {
	/* Updates the propensity of reactions in a cell
    
    Inputs
    ------
	propensity, 2d array of propensity vectors for each cell
	propensity_sums, 1d array of propensity sums for each cell
	cell_idx, index of cell to update propensity array of
	wildtype_populations, 1d array of wildtype population size of each cell
	sss_populations, 1d array of SSS population size of each cell
    degradation_rate, degradation rate
	diffusion_rate, diffusion rate
	nucleus_control_factor, nucleus conrol factor
    target_population, target population
	rate_difference, rate difference */

	int cell_wildtype_population = wildtype_populations[cell_idx];
	int cell_sss_population = sss_populations[cell_idx];
	double replication_rate = fmax(0.0, degradation_rate + nucleus_control_factor * (target_population - cell_wildtype_population - rate_difference*cell_sss_population));
	
	propensity[cell_idx][0] = (degradation_rate + rate_difference) * cell_wildtype_population;
	propensity[cell_idx][1] = (replication_rate + rate_difference) * cell_wildtype_population;
	propensity[cell_idx][2] = diffusion_rate * cell_wildtype_population;
	propensity[cell_idx][3] = degradation_rate * cell_sss_population;
	propensity[cell_idx][4] = replication_rate * cell_sss_population;
	propensity[cell_idx][5] = diffusion_rate * cell_sss_population;

	propensity_sums[cell_idx] = 0.00;
	for (int i=0; i<6; ++i) {propensity_sums[cell_idx] += propensity[cell_idx][i];}
	return;
}

int main(int argc, char *argv[]) {
    int seed;
    double log_site_std_mutation_rate = LOG_SITE_STD_MUTATION_RATE;
    double degradation_rate = DEGRADATION_RATE;
	double diffusion_rate = DIFFUSION_RATE;
    double nucleus_control_factor = NUCLEUS_CONTROL_FACTOR;
    int target_population = TARGET_POP;
	double rate_difference = RATE_DIFFERENCE;

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
    long double site_std_mutation_rate = pow(10, - log_site_std_mutation_rate);

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
	fprintf(fp_parameters, "target_pop,%d\n", TARGET_POP);
	fprintf(fp_parameters, "density,%e\n", DENSITY);
	fprintf(fp_parameters, "len_genome,%d\n", LEN_GENOME);

	fprintf(fp_parameters, "sim_length,%e\n", SIM_LENGTH);
	fprintf(fp_parameters, "introduce_after,%e\n", INTRODUCE_AFTER);
	fprintf(fp_parameters, "recording_space,%e\n", RECORDING_SPACE);

	fprintf(fp_parameters, "max_n_events,%d\n", MAX_N_EVENTS);
	fprintf(fp_parameters, "max_mutants,%d\n", MAX_MUTANTS);
	fclose(fp_parameters);

	// Set up files to write population data in
	char sss_population_filename[50];
	sprintf(sss_population_filename, "sss_sim_populations%d.txt", seed);
	FILE *fp_sss_population = fopen(sss_population_filename, "w");
	fprintf(fp_sss_population, "sim,cell,t,wildtype_population,sss_population\n");
	fclose(fp_sss_population);

	// Set up files to write site frequency spectrum data in
	char sss_sfs_filename[40];
	sprintf(sss_sfs_filename, "sss_sim_site_frequency_spectrum%d.txt", seed);
	FILE *fp_sss_sfs = fopen(sss_sfs_filename, "w");
	fprintf(fp_sss_sfs, "sim,t,mut_id1,mut_id2,...\n");
	fclose(fp_sss_sfs);

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
		wildtype_state[k] = malloc(target_population * sizeof(int*));
		for (int i=0; i<target_population; ++i) {wildtype_state[k][i] = calloc(1, sizeof(int));}
		wildtype_populations[k] = target_population;
	
		// Initialise mutant counts
		mutant_counts[k] = calloc(1, sizeof(int));

		// Compute initial wildtype propensity
		wildtype_propensity[k] = malloc(3 * sizeof(double));
		wildtype_propensity_update(wildtype_propensity, wildtype_propensity_sums, k, wildtype_populations, degradation_rate, diffusion_rate, nucleus_control_factor, target_population, rate_difference);
	}

	double wildtype_propensity_sum_across_cells;
	double current_time = 0.0;
	double recording_time = 0.0;
	int n_event = 0;

	// Record initial data
	write_data_to_file_pre_introduce(wildtype_populations, mutant_counts, 0, recording_time, sss_population_filename, sss_sfs_filename);
	recording_time += RECORDING_SPACE;

	// Gillespie algorithm until time threshold reached
	while (current_time<INTRODUCE_AFTER && n_event < MAX_N_EVENTS) {
		
		// Realise event according to propensity
		wildtype_propensity_sum_across_cells = 0;
		for (int k=0; k<CELLS; ++k) {wildtype_propensity_sum_across_cells += wildtype_propensity_sums[k];}
		current_time += -log(RND) / wildtype_propensity_sum_across_cells;
		wildtype_gillespie_event(rng, wildtype_propensity, wildtype_propensity_sums, wildtype_state, wildtype_populations, mutant_counts, site_std_mutation_rate, degradation_rate, diffusion_rate, nucleus_control_factor, target_population, rate_difference);
		n_event++;

		// Record data at recording time
		if (current_time>=recording_time) {
			compact_relabel_wildtype_mutations(mutant_counts, wildtype_state, wildtype_populations);
			write_data_to_file_pre_introduce(wildtype_populations, mutant_counts, 0, recording_time, sss_population_filename, sss_sfs_filename);

			recording_time += RECORDING_SPACE;
		}
	}

	int introduce_cell_idx = CELLS / 2;
	while (wildtype_populations[introduce_cell_idx]==0) {
		// Realise event according to propensity
		wildtype_propensity_sum_across_cells = 0;
		for (int k=0; k<CELLS; ++k) {wildtype_propensity_sum_across_cells += wildtype_propensity_sums[k];}
		current_time += -log(RND) / wildtype_propensity_sum_across_cells;
		wildtype_gillespie_event(rng, wildtype_propensity, wildtype_propensity_sums, wildtype_state, wildtype_populations, mutant_counts, site_std_mutation_rate, degradation_rate, diffusion_rate, nucleus_control_factor, target_population, rate_difference);
		n_event++;
	}

	// Free memory
	for (int k=0; k<CELLS; ++k) {free(wildtype_propensity[k]);}
	free(wildtype_propensity);
	free(wildtype_propensity_sums);

	/* sss_state contains mutational information about RA individuals
	sss_state[k][i][0] is the number of standard mutations which individual i possesses in cell k
	sss_state[k][i][1:] are the identities of the standard mutations which individual i possesses in cell k */
	int*** sss_state = malloc(CELLS * sizeof(int**));
	/* sss_populations[k] is the population size of RA individuals in cell k */
	int* sss_populations = malloc(CELLS * sizeof(int));

	// Introduce 1 SSS individual at central cell
	introduce_ra_or_ssd(rng, wildtype_state, sss_state, wildtype_populations, sss_populations, introduce_cell_idx);			

	// Initialise propensity
	/* propensity[k][0] is the propensity of degradation of a wildtype individual in cell k
	propensity[k][1] is the propensity of replication of a wildtype individual in cell k
	propensity[k][2] is the propensity of diffusion of a wildtype individual in cell k
	propensity[k][3] is the propensity of degradation of an RA/SSS individual in cell k
	propensity[k][4] is the propensity of replication of an RA/SSS individual in cell k
	propensity[k][5] is the propensity of diffusion of an RA/SSS individual in cell k*/
	double** propensity = malloc(CELLS * sizeof(double*));
	/* propensity_sums[k] is the sum of the propensity vector of cell k */
	double* propensity_sums = malloc(CELLS * sizeof(double));
	for (int k=0; k<CELLS; ++k) {
		propensity[k] = malloc(6 * sizeof(double));
		propensity_update(propensity, propensity_sums, k, wildtype_populations, sss_populations, degradation_rate, diffusion_rate, nucleus_control_factor, target_population, rate_difference);
	}

	double propensity_sum_across_cells;
	
	// Simulate
	while (current_time<SIM_LENGTH && n_event<MAX_N_EVENTS) {
		if (n_event%100==0) {printf("n_event = %d\n", n_event);}
		propensity_sum_across_cells = 0;
		for (int k=0; k<CELLS; ++k) {propensity_sum_across_cells += propensity_sums[k];}
		current_time += -log(RND) / propensity_sum_across_cells;
		gillespie_event(rng, propensity, propensity_sums, wildtype_state, sss_state, mutant_counts, wildtype_populations, sss_populations, site_std_mutation_rate, degradation_rate, diffusion_rate, nucleus_control_factor, target_population, rate_difference);
		n_event++;

		// In case of extinction of SSS individuals, record data before beginning new simulation
		int total_sss_population = 0;
		for (int k=0; k<CELLS; ++k) {total_sss_population += sss_populations[k];}
		if (total_sss_population==0) {
			compact_relabel_mutations(mutant_counts, wildtype_state, sss_state, wildtype_populations, sss_populations);
			write_data_to_file(wildtype_populations, sss_populations, mutant_counts, 0, recording_time, sss_population_filename, sss_sfs_filename);
			break;
		}

		// Record data at recording time
		if (current_time>=recording_time){
			compact_relabel_mutations(mutant_counts, wildtype_state, sss_state, wildtype_populations, sss_populations);
			write_data_to_file(wildtype_populations, sss_populations, mutant_counts, 0, recording_time, sss_population_filename, sss_sfs_filename);

			recording_time += RECORDING_SPACE;
		}
	}

	// Free memory
	for (int k=0; k<CELLS; ++k){
		for (int i=0; i<wildtype_populations[k]; ++i){free(wildtype_state[k][i]);}
		free(wildtype_state[k]);
		for (int i=0; i<sss_populations[k]; ++i){free(sss_state[k][i]);}
		free(sss_state[k]);

		free(mutant_counts[k]);

		free(propensity[k]);
	}
	free(wildtype_state);
	free(wildtype_populations);
	free(sss_state);
	free(sss_populations);

	free(mutant_counts);

	free(propensity);
	free(propensity_sums);

    gsl_rng_free(rng);
	return 0;
}