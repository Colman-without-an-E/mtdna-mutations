#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/stat.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "../include/parameters.h"
#include "../include/lib_sim.h"
#include "../include/_debug.h"

#define RND gsl_rng_uniform_pos(rng) // generate number from Unif(0,1)

// Define propensity updating rule for RA simulation
void propensity_update(double** propensity, double* propensity_sums, int cell_idx, int* wildtype_populations, int* ra_populations, double degradation_rate, double diffusion_rate, double nucleus_control_factor, int target_population, double replicative_advantage) {
	/* Updates the propensity of reactions in a cell
    
    Inputs
    ------
	propensity, 2d array of propensity vectors for each cell
	propensity_sums, 1d array of propensity sums for each cell
	cell_idx, index of cell to update propensity array of
	wildtype_populations, 1d array of wildtype population size of each cell
	ra_populations, 1d array of RA population size of each cell
    degradation_rate, degradation rate
	diffusion_rate, diffusion rate
    density, density
    target_population, target population
	replicative_advantage, replicative advantage */

	int cell_wildtype_population = wildtype_populations[cell_idx];
	int cell_ra_population = ra_populations[cell_idx];
	double wildtype_replication_rate = degradation_rate + nucleus_control_factor * (target_population - cell_wildtype_population - cell_ra_population);
	double ra_replication_rate = wildtype_replication_rate + replicative_advantage;
	
	propensity[cell_idx][0] = degradation_rate * cell_wildtype_population;
	propensity[cell_idx][1] = fmax(0.0, wildtype_replication_rate * cell_wildtype_population);
	propensity[cell_idx][2] = diffusion_rate * cell_wildtype_population;
	propensity[cell_idx][3] = degradation_rate * cell_ra_population;
	propensity[cell_idx][4] = fmax(0.0, ra_replication_rate * cell_ra_population);
	propensity[cell_idx][5] = diffusion_rate * cell_ra_population;

	propensity_sums[cell_idx] = 0.00;
	for (int i=0; i<6; ++i) {propensity_sums[cell_idx] += propensity[cell_idx][i];}
	return;
}

int main(int argc, char *argv[]) {
	/* set up GSL RNG */
	gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);	
	/* end of GSL setup */

    int seed = SEED;
    double log_site_std_mutation_rate = LOG_SITE_STD_MUTATION_RATE;
    double degradation_rate = DEGRADATION_RATE;
	double diffusion_rate = DIFFUSION_RATE;
    double nucleus_control_factor = NUCLEUS_CONTROL_FACTOR;
    int target_population = TARGET_POP;
	double replicative_advantage = REPLICATIVE_ADVANTAGE;

	// Generate replication rate and mutation rate for this run
    long double site_std_mutation_rate = pow(10, - log_site_std_mutation_rate);

	// Create directory with name current time, to store simulation results
	time_t t = time(NULL);
    struct tm *tm_info = localtime(&t);
    char dir_name[20];
    strftime(dir_name, 20, "%Y-%m-%d_%H-%M-%S", tm_info);
    char dir_loc[256];
    sprintf(dir_loc, "%s/%s", DIR_LOC, dir_name);
    puts(dir_loc);
    if(mkdir(dir_loc)) {
        perror("Error");
        return 1;
    }

    // Change to directory which stores simulation results
    if (chdir(dir_loc)) {
        perror("Error");
        return 1;
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
	fprintf(fp_parameters, "replicative_advantage,%e\n", REPLICATIVE_ADVANTAGE);
	fprintf(fp_parameters, "len_genome,%d\n", LEN_GENOME);

	fprintf(fp_parameters, "seed,%d\n", SEED);
	fprintf(fp_parameters, "n_sims,%d\n", N_SIMS);
	fprintf(fp_parameters, "sim_length,%e\n", SIM_LENGTH);
	fprintf(fp_parameters, "recording_space,%e\n", RECORDING_SPACE);
	fprintf(fp_parameters, "n_bins,%d\n", N_BINS);

	fprintf(fp_parameters, "max_n_events,%d\n", MAX_N_EVENTS);
	fprintf(fp_parameters, "max_mutants,%d\n", MAX_MUTANTS);
	fclose(fp_parameters);

	// Set up files to write population data in
	char ra_population_filename[30] = "ra_sim_populations.txt";
	FILE *fp_ra_population = fopen(ra_population_filename, "w");
	fprintf(fp_ra_population, "sim,cell,t,wildtype_population,ra_population\n");
	fclose(fp_ra_population);

	// Set up files to write site frequency spectrum data in
	char ra_sfs_filename[40] = "ra_sim_site_frequency_spectrum.txt";
	FILE *fp_ra_sfs = fopen(ra_sfs_filename, "w");
	fprintf(fp_ra_sfs, "sim,cell,t");
	double bin_lb = 0.0; // bin lower bound
	double bin_width = 1.0 / N_BINS;
	for (int bin=0; bin<N_BINS; ++bin) {
		fprintf(fp_ra_sfs, ",%.2f", bin_lb);
		bin_lb += bin_width;
	}
	fprintf(fp_ra_sfs, "\n");
	fclose(fp_ra_sfs);

	// Benchmark each batch of simulations
	clock_t tick, tock;
	double time_elapsed;

	for(int sim = 0; sim < N_SIMS; ++sim) {
		if (sim % BATCH_SIZE == 0){tick = clock();}
	
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
			wildtype_propensity_update(wildtype_propensity, wildtype_propensity_sums, k, wildtype_populations, degradation_rate, diffusion_rate, nucleus_control_factor, target_population);
		}

		int max_std_mutant_count;
		int cell_with_highest_heteroplasmy;
		double wildtype_propensity_sum_across_cells;
		double current_time = 0.0;
		double recording_time = 0.0;
		int n_event = 0;

		// Gillespie algorithm until time threshold reached
		while (current_time<INTRODUCE_AT && n_event < MAX_N_EVENTS) {
			// Realise event according to propensity
			// printf("t = %.3f n_event = %d\n");
			wildtype_propensity_sum_across_cells = 0;
			for (int k=0; k<CELLS; ++k) {wildtype_propensity_sum_across_cells += wildtype_propensity_sums[k];}
			current_time += -log(RND) / wildtype_propensity_sum_across_cells;
			wildtype_gillespie_event(rng, wildtype_propensity, wildtype_propensity_sums, wildtype_state, wildtype_populations, mutant_counts, site_std_mutation_rate, degradation_rate, diffusion_rate, nucleus_control_factor, target_population);
			n_event++;

			// Record data at recording time
			if (current_time>(recording_time - 1E-12)) {
				// printf("t = %.1f\n", current_time);
				compact_relabel_wildtype_mutations(mutant_counts, wildtype_state, wildtype_populations);
				write_data_to_file_pre_introduce(wildtype_populations, mutant_counts, sim, recording_time, ra_population_filename, ra_sfs_filename);

				recording_time += RECORDING_SPACE;
			}
		}

		// Free memory
		for (int k=0; k<CELLS; ++k) {free(wildtype_propensity[k]);}
		free(wildtype_propensity);
		free(wildtype_propensity_sums);

		/* ra_state contains mutational information about RA individuals
		ra_state[k][i][0] is the number of standard mutations which individual i possesses in cell k
		ra_state[k][i][1:] are the identities of the standard mutations which individual i possesses in cell k */
		int*** ra_state = malloc(CELLS * sizeof(int**));
		/* ra_populations[k] is the population size of RA individuals in cell k */
		int* ra_populations = malloc(CELLS * sizeof(int));

		// Introduce 1 RA individual at central cell
		int introduce_cell_idx = CELLS / 2;
		introduce_ra_or_ssd(rng, wildtype_state, ra_state, wildtype_populations, ra_populations, cell_with_highest_heteroplasmy);			

		// Initialise propensity
		/* propensity[k][0] is the propensity of degradation of a wildtype individual in cell k
		propensity[k][1] is the propensity of replication of a wildtype individual in cell k
		propensity[k][2] is the propensity of diffusion of a wildtype individual in cell k
		propensity[k][3] is the propensity of degradation of an RA/SSD individual in cell k
		propensity[k][4] is the propensity of replication of an RA/SSD individual in cell k
		propensity[k][5] is the propensity of diffusion of an RA/SSD individual in cell k*/
		double** propensity = malloc(CELLS * sizeof(double*));
		/* propensity_sums[k] is the sum of the propensity vector of cell k */
		double* propensity_sums = malloc(CELLS * sizeof(double));
		for (int k=0; k<CELLS; ++k) {
			propensity[k] = malloc(6 * sizeof(double));
			propensity_update(propensity, propensity_sums, k, wildtype_populations, ra_populations, degradation_rate, diffusion_rate, nucleus_control_factor, target_population, replicative_advantage);
		}

		double propensity_sum_across_cells;

		// Simulate
		while (current_time<SIM_LENGTH && n_event<MAX_N_EVENTS) {
			propensity_sum_across_cells = 0;
			for (int k=0; k<CELLS; ++k) {propensity_sum_across_cells += propensity_sums[k];}
			current_time += -log(RND) / propensity_sum_across_cells;
			gillespie_event(rng, propensity, propensity_sums, wildtype_state, ra_state, mutant_counts, wildtype_populations, ra_populations, site_std_mutation_rate, degradation_rate, diffusion_rate, nucleus_control_factor, target_population, replicative_advantage);
			n_event++;

			// In case of extinction of RA individuals, record data before beginning new simulation
			int total_ra_population = 0;
			for (int k=0; k<CELLS; ++k) {total_ra_population += ra_populations[k];}
			if (total_ra_population==0) {
				compact_relabel_mutations(mutant_counts, wildtype_state, ra_state, wildtype_populations, ra_populations);
				write_data_to_file(wildtype_populations, ra_populations, mutant_counts, sim, recording_time, ra_population_filename, ra_sfs_filename);
				break;
			}

			// Record data at recording time
			if (current_time>(recording_time - 1E-12)){
				// printf("t = %.1f\n", current_time);
				compact_relabel_mutations(mutant_counts, wildtype_state, ra_state, wildtype_populations, ra_populations);
				write_data_to_file(wildtype_populations, ra_populations, mutant_counts, sim, recording_time, ra_population_filename, ra_sfs_filename);

				recording_time += RECORDING_SPACE;
			}
		}

		// Free memory
		for (int k=0; k<CELLS; ++k){
			for (int i=0; i<wildtype_populations[k]; ++i){free(wildtype_state[k][i]);}
			free(wildtype_state[k]);
			for (int i=0; i<ra_populations[k]; ++i){free(ra_state[k][i]);}
			free(ra_state[k]);

			free(mutant_counts[k]);

			free(propensity[k]);
		}
		free(wildtype_state);
		free(wildtype_populations);
		free(ra_state);
		free(ra_populations);

		free(mutant_counts);

		free(propensity);
		free(propensity_sums);
		
		// Print time usage
		if ((sim+1) % BATCH_SIZE == 0) {
			tock = clock();
			time_elapsed = ((double) (tock - tick)) / CLOCKS_PER_SEC;
			printf("Simulation %d complete. Batch %d - Time taken: %f seconds\n", sim+1, (sim+1)/BATCH_SIZE, time_elapsed);
		}
	}
    gsl_rng_free(rng);
	return 0;
}