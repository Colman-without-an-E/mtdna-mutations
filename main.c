#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <limits.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics.h>

#define RND gsl_rng_uniform_pos(rng) // generate number from Unif(0,1)
#define CELLS 2
#define N_SIMS 2
#define BATCH_SIZE 1
#define SIM_LENGTH 100.00
#define RECORDING_SPACE 10.0
#define N_BINS 20
#define YEARS 365.00
#define WEEKS 7.00
#define MAX_MUTANTS 100000
#define LEN_GENOME 16569

void print_state(int** state, int nrow) {
	/* Print the state of the system for debugging
	
	Inputs
	------
	state, 2d array of the state
	nrow, number of rows */

	int n_mutations_to_print;
	for (int i=0; i<nrow; ++i) {
		n_mutations_to_print = state[i][0];
		for (int j=0; j<=n_mutations_to_print; ++j) {
			printf("%d ", state[i][j]);
		}
		printf("\n");
	}
	return;
}

void print_mutant_counts(int** mutant_counts) {
	/* Print mutant counts for debugging
	
	Inputs
	------
	mutant_counts, array of mutant counts */

	int ncols = mutant_counts[0][0] + 1;
	for (int k=0; k<CELLS; ++k) {
		for (int j=0; j<ncols; ++j) {printf("%d ", mutant_counts[k][j]);}
		printf("\n");
	}
	return;
}

int negative_int_in_mutant_counts(int** mutant_counts, int mutation_counter) {
	/* Returns column index of first negative integer in mutant_counts, otherwise returns 0. For debugging
	
	Inputs
	------
	mutant_counts, mutant counts
	mutation_counter, latest mutation identity */
	for (int k=0; k<CELLS; ++k) {
		for (int j=0; j<=mutation_counter; ++j) {
			if (mutant_counts[k][j]<0) {return j;}
		}
	}
	return 0;
}

void wildtype_propensity_update(double** propensity, double* propensity_sums, int cell_idx, int* wildtype_populations, double degradation_rate, double diffusion_rate, double nucleus_control_factor, int target_population) {
	/* Updates the wildtype propensity of reactions in a cell
    
    Inputs
    ------
    propensity, 2d array of propensity vectors for each cell
	propensity_sums, 1d array of propensity sums for each cell
	cell_idx, index of cell to update propensity array of
	wildtype_populations, 1d array of wildtype population size of each cell
    degredation_rate, degradation rate
	diffusion_rate, diffusion rate
    nucleus_control_factor, nucleus control factor
    target_population, target population */

	int cell_wildtype_population = wildtype_populations[cell_idx];
	double wildtype_replication_rate = degradation_rate + nucleus_control_factor * (target_population - cell_wildtype_population);

	propensity[cell_idx][0] = degradation_rate * cell_wildtype_population;
	propensity[cell_idx][1] = wildtype_replication_rate * cell_wildtype_population;
	propensity[cell_idx][2] = diffusion_rate * cell_wildtype_population;

	propensity_sums[cell_idx] = 0.00;
	for (int i=0; i<3; ++i){propensity_sums[cell_idx] += propensity[cell_idx][i];}
	return;
}

void ra_propensity_update(double** propensity, double* propensity_sums, int cell_idx, int* wildtype_populations, int* ra_populations, double degradation_rate, double diffusion_rate, double nucleus_control_factor, int target_population, double replicative_advantage) {
	/* Updates the propensity of reactions in a cell
    
    Inputs
    ------
	propensity, 2d array of propensity vectors for each cell
	propensity_sums, 1d array of propensity sums for each cell
	cell_idx, index of cell to update propensity array of
	wildtype_populations, 1d array of wildtype population size of each cell
	ra_or_ssd_populations, 1d array of RA population size of each cell
    ssd_populations, 1d array of SSD population size of each cell
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
	propensity[cell_idx][1] = wildtype_replication_rate * cell_wildtype_population;
	propensity[cell_idx][2] = diffusion_rate * cell_wildtype_population;
	propensity[cell_idx][3] = degradation_rate * cell_ra_population;
	propensity[cell_idx][4] = ra_replication_rate * cell_ra_population;
	propensity[cell_idx][5] = diffusion_rate * cell_ra_population;

	propensity_sums[cell_idx] = 0.00;
	for (int i=0; i<6; ++i) {propensity_sums[cell_idx] += propensity[cell_idx][i];}
	return;
}

void ssd_propensity_update(double** propensity, double* propensity_sums, int cell_idx, int* wildtype_populations, int* ssd_populations, double degradation_rate, double diffusion_rate, double nucleus_control_factor, double density, int target_population) {
	/* Updates the propensity of reactions in a cell
    
    Inputs
    ------
	propensity, 2d array of propensity vectors for each cell
	propensity_sums, 1d array of propensity sums for each cell
	cell_idx, index of cell to update propensity array of
	wildtype_populations, 1d array of wildtype population size of each cell
    ssd_populations, 1d array of SSD population size of each cell
    degradation_rate, degradation rate
	diffusion_rate, diffusion rate
    nucleus_control_factor, nucleus control factor
    density, density
    target_population, target population */

	int cell_wildtype_population = wildtype_populations[cell_idx];
	int cell_ssd_population = ssd_populations[cell_idx];
	double wildtype_replication_rate = degradation_rate + nucleus_control_factor * (target_population - cell_wildtype_population - density * cell_ssd_population);
	double ssd_replication_rate = wildtype_replication_rate;

	propensity[cell_idx][0] = degradation_rate * cell_wildtype_population;
	propensity[cell_idx][1] = wildtype_replication_rate * cell_wildtype_population;
	propensity[cell_idx][2] = diffusion_rate * cell_wildtype_population;
	propensity[cell_idx][3] = degradation_rate * cell_ssd_population;
	propensity[cell_idx][4] = ssd_replication_rate * cell_ssd_population;
	propensity[cell_idx][5] = diffusion_rate * cell_ssd_population;

	propensity_sums[cell_idx] = 0.00;
	for (int i=0; i<6; ++i) {propensity_sums[cell_idx] += propensity[cell_idx][i];}
	return;
}

int weighted_sample(const gsl_rng* rng, int n, double* weights) {
	/* Performs weighted sampling
	
	Inputs
	------
	rng, random number generator
    n, number of samples
    weights, array of weights
	*/
    double weight_sum = 0.00;
    for (int i=0; i<n; ++i){weight_sum += weights[i];}
    
    double random_point = RND * weight_sum;
    double cum_sum=0;
    
    for (int i=0; i<n-1; ++i) {
        cum_sum += weights[i];
        if (cum_sum > random_point) {
            return i;
        }
    }
    return n-1;
}

void state_update_degrade(int*** state, int cell_idx, int degrade_row, int* populations, int** mutant_counts) {
    /* Updates state when degradation of one individual occurs
    
    Inputs
    ------
    state, SSD or non-SSD state of system
	cell_idx, index of cell where degradation occurs
    degrade_row, row index of individual to be degraded
	populations, population count of each cell
    mutant_counts, 2d array of mutant counts for each cell */

	int cell_population = populations[cell_idx];

   	// Decrease mutant counts
	int n_mutations_to_die = state[cell_idx][degrade_row][0];
	int mutation_id;
	for (int j=1; j<=n_mutations_to_die; j++) {
		mutation_id = state[cell_idx][degrade_row][j];
		mutant_counts[cell_idx][mutation_id]--;
	}
	
	// Copy last row of state to degrade_row
	int n_mutations_to_copy = state[cell_idx][cell_population-1][0];
	state[cell_idx][degrade_row] = realloc(state[cell_idx][degrade_row], (n_mutations_to_copy+1) * sizeof(int));
	for (int j=0; j<=n_mutations_to_copy; j++) {state[cell_idx][degrade_row][j] = state[cell_idx][cell_population-1][j];}

	// Empty last row of state
	free(state[cell_idx][cell_population-1]);
	state[cell_idx] = realloc(state[cell_idx], (cell_population-1) * sizeof(int*));
	populations[cell_idx]--;
	return;
}

void state_update_replicate(int*** state, int cell_idx, int replicate_row, int* populations, int** mutant_counts, int n_new_std_mutations) {
    /* Updates state when replication of one individual occurs

    Inputs
    ------
	state, SSD or non-SSD state of system
	cell_idx, index of cell where degradation occurs
    replicate_row, row index of individual to be replicated
	populations, population count of each cell
    mutant_counts, 2d array of mutant counts of each cell
    n_new_std_mutations, number of new standard mutations that occur during replication */

	int cell_population = populations[cell_idx];

    int n_existing_std_mutations = state[cell_idx][replicate_row][0];
    int n_total_std_mutations = n_existing_std_mutations + n_new_std_mutations;
    if (n_total_std_mutations > MAX_MUTANTS) {
        printf("Maximum number of standard mutations allowed reached. Increase MAX_MUTANTS. Exiting.\n");
        exit(99);
    }

	// Reallocate state for new row
    state[cell_idx] = realloc(state[cell_idx], (cell_population + 1) * sizeof(int*));
    state[cell_idx][cell_population] = malloc((n_total_std_mutations+1) * sizeof(int));

	// Reallocate and initialise to 0 n_new_std_mutations new columns for mutant_counts
	if (n_new_std_mutations) {
		int ncol_new = mutant_counts[0][0] + 1 + n_new_std_mutations;
		for (int k=0; k<CELLS; ++k) {
			mutant_counts[k] = realloc(mutant_counts[k], ncol_new * sizeof(int));
			for (int j=mutant_counts[0][0]+1; j<ncol_new; ++j) {mutant_counts[k][j] = 0;}
		}
	}

    // Copy content to the newly allocated row
	state[cell_idx][cell_population][0] = n_total_std_mutations;
    for (int j=1; j<=n_existing_std_mutations; ++j) { // copy existing mutations id
		int mutation_id = state[cell_idx][replicate_row][j];
		state[cell_idx][cell_population][j] = mutation_id;

		mutant_counts[cell_idx][mutation_id]++;
	}
    for (int j=n_existing_std_mutations+1; j<=n_total_std_mutations; ++j) { // include new mutations
		mutant_counts[0][0]++; // mutation id
		state[cell_idx][cell_population][j] = mutant_counts[0][0];
		mutant_counts[cell_idx][mutant_counts[0][0]] = 1;
	}
	populations[cell_idx]++;
	return;
}

void state_update_diffuse(int*** state, int cell_idx_from, int cell_idx_to, int diffuse_row, int* populations, int** mutant_counts) {
    /* Updates state when diffusion of one individual occurs

    Inputs
    ------
	state, SSD or non-SSD state of system
	cell_idx_from, index of cell where diffusion occurs
	cell_idx_to, index of cell where the individual diffuses to
    diffuse_row, row index of individual to diffuse
	populations, population count of each cell
    mutant_counts, 2d array of mutant counts of each cell */
	
	int cell_population_from = populations[cell_idx_from];
	int cell_population_to = populations[cell_idx_to];

	// Reallocate for new row of state to be diffused into
	int n_mutations_to_copy = state[cell_idx_from][diffuse_row][0];
    state[cell_idx_to] = realloc(state[cell_idx_to], (cell_population_to + 1) * sizeof(int*));
    state[cell_idx_to][cell_population_to] = malloc((n_mutations_to_copy+1) * sizeof(int));

    // Copy content to the newly allocated row
    for (int j=0; j<=n_mutations_to_copy; ++j) {
		int mutation_id = state[cell_idx_from][diffuse_row][j];
		state[cell_idx_to][cell_population_to][j] = mutation_id;
		if (mutation_id) {mutant_counts[cell_idx_to][mutation_id]++;}
	}
	populations[cell_idx_to]++;

	// Remove original row
	state_update_degrade(state, cell_idx_from, diffuse_row, populations, mutant_counts);
	return;
}

int choose_neighbouring_cell(const gsl_rng* rng, int cell_idx) {
	/* Returns index of neighbouring cell with equal probability
	
	Inputs
	------
	rng, rng
	cell_idx, index of cell whose neighbouring cell index is to be found */

	int neighbour_cell_idx = cell_idx;
	switch (cell_idx) {
		case 0: // left boundary case
			neighbour_cell_idx++;
			break;
		case CELLS-1: // right boundary case
			neighbour_cell_idx--;
			break;
		default: // generic case
			if(gsl_rng_uniform_int(rng, 2)) {neighbour_cell_idx++;} else {neighbour_cell_idx--;}
			break;
	}
	return neighbour_cell_idx;
}

void non_ssd_gillespie_event(const gsl_rng* rng, double** propensity, double* propensity_sums, int*** wildtype_state, int* wildtype_populations, int** mutant_counts, long double site_std_mutation_rate, double degradation_rate, double diffusion_rate, double nucleus_control_factor, int target_population) {
	/* Realises chosen event and updates system 
	
	Inputs
	------
	rng, random number generator
    propensity, 2d array of propensity vectors for each cell
    propensity_sums, 1d array of propensity sums for each cell
    wildtype_state, 3d array of standard mutant state of the system
	wildtype_populations, 1d array of wildtype population size of each cell
	mutant_counts, 2d array of mutant counts
	site_std_mutation_rate, standard mutation rate at each site
	degradation_rate, degradation rate
	diffusion_rate, diffusion rate
	nucleus_control_factor, nucleus control factor
	target_population, target population */

	// Choose which cell event takes place in
	int cell_idx = weighted_sample(rng, CELLS, propensity_sums);

	// Choose event
	int event_id = weighted_sample(rng, 3, propensity[cell_idx]);

	int index_die;
	int index_born;
	int index_diffuse;
	int n_new_std_mutations;
	int cell_idx_diffuse_to = -1;
	switch (event_id) {
		case 0: // wildtype degradation
			index_die = floor(RND * wildtype_populations[cell_idx]);
			state_update_degrade(wildtype_state, cell_idx, index_die, wildtype_populations, mutant_counts);
			break;
		case 1: // wildtype replication
			index_born = floor(RND * wildtype_populations[cell_idx]);
			n_new_std_mutations = gsl_ran_binomial(rng, site_std_mutation_rate, LEN_GENOME);
			state_update_replicate(wildtype_state, cell_idx, index_born, wildtype_populations, mutant_counts, n_new_std_mutations);
			break;
		case 2: // wildtype diffusion
			// Choose cell to diffuse to
			cell_idx_diffuse_to = choose_neighbouring_cell(rng, cell_idx);
			index_diffuse = floor(RND * wildtype_populations[cell_idx]);
			state_update_diffuse(wildtype_state, cell_idx, cell_idx_diffuse_to, index_diffuse, wildtype_populations, mutant_counts);
			break;
	}

	// Update propensity vectors and sum at cells where reaction occured, or with individuals diffused to
	wildtype_propensity_update(propensity, propensity_sums, cell_idx, wildtype_populations, degradation_rate, diffusion_rate, nucleus_control_factor, target_population);
	if (cell_idx_diffuse_to>=0) {wildtype_propensity_update(propensity, propensity_sums, cell_idx_diffuse_to, wildtype_populations, degradation_rate, diffusion_rate, nucleus_control_factor, target_population);}
	return;
}

void gillespie_event(const gsl_rng* rng, int ssd_sim, double** propensity, double* propensity_sums, int*** wildtype_state, int*** ra_or_ssd_state, int** mutant_counts, int* wildtype_populations, int* ra_or_ssd_populations, long double site_std_mutation_rate, double degradation_rate, double diffusion_rate, double nucleus_control_factor, double density, int target_population, double replicative_advantage) {
	/* Realise chosen event and updates system with RAs
	
	Inputs
	------
	rng, random number generator
	ssd_sim, 0 or 1; 0 for RA simulation and 1 for SSD simulation
    propensity, 2d array of propensity vectors for each cell
    propensity_sums, 1d array of propensity sums for each cell
    wildtype_state, 3d array of wildtype state of the system
    ra_or_ssd_state, 3d array of RA/SSD state of the system
	mutant_counts, 2d array of mutant counts of each cell
	wildtype_populations, 1d array of wildtype population size of each cell
	ra_or_ssd_populations, 1d array of RA/SSD population size of each cell
	site_std_mutation_rate, standard mutation rate at each site
	degradation_rate, degradation rate
	diffusion_rate, diffusion rate
	nucleus_control_factor, nucleus control factor
	density, density
	target_population, target population
	replicative_advantage, replicative advantage */

	// Choose which cell event takes place in
	int cell_idx = weighted_sample(rng, CELLS, propensity_sums);

	// Choose event
	int event_id = weighted_sample(rng, 6, propensity[cell_idx]);

	int index_die;
	int index_born;
	int index_diffuse;
	int n_new_std_mutations;
	int cell_idx_diffuse_to = -1;
	switch (event_id) {
		case 0: // wildtype degradation
			index_die = floor(RND * wildtype_populations[cell_idx]);
			state_update_degrade(wildtype_state, cell_idx, index_die, wildtype_populations, mutant_counts);
			break;
		case 1: // wildtype replication
			index_born = floor(RND * wildtype_populations[cell_idx]);
			n_new_std_mutations = gsl_ran_binomial(rng, site_std_mutation_rate, LEN_GENOME);
			state_update_replicate(wildtype_state, cell_idx, index_born, wildtype_populations, mutant_counts, n_new_std_mutations);
			break;
		case 2: // wildtype diffusion
			// Choose cell to diffuse to
			cell_idx_diffuse_to = choose_neighbouring_cell(rng, cell_idx);
			index_diffuse = floor(RND * wildtype_populations[cell_idx]);
			state_update_diffuse(wildtype_state, cell_idx, cell_idx_diffuse_to, index_diffuse, wildtype_populations, mutant_counts);
			break;
		case 3: // RA degradation
			index_die = floor(RND * ra_or_ssd_populations[cell_idx]);
			state_update_degrade(ra_or_ssd_state, cell_idx, index_die, ra_or_ssd_populations, mutant_counts);
			break;
		case 4: // RA replication
			index_born = floor(RND * ra_or_ssd_populations[cell_idx]);
			n_new_std_mutations = gsl_ran_binomial(rng, site_std_mutation_rate, LEN_GENOME);
			state_update_replicate(ra_or_ssd_state, cell_idx, index_born, ra_or_ssd_populations, mutant_counts, n_new_std_mutations);
			break;
		case 5: // RA diffusion
			// Choose cell to diffuse to
			cell_idx_diffuse_to = choose_neighbouring_cell(rng, cell_idx);
			index_diffuse = floor(RND * ra_or_ssd_populations[cell_idx]);
			state_update_diffuse(ra_or_ssd_state, cell_idx, cell_idx_diffuse_to, index_diffuse, ra_or_ssd_populations, mutant_counts);
			break;
	}

	// Update propensity vectors and sum at cells where reaction occured, and with individuals diffused to
	if (ssd_sim) {
		ssd_propensity_update(propensity, propensity_sums, cell_idx, wildtype_populations, ra_or_ssd_populations, degradation_rate, diffusion_rate, nucleus_control_factor, density, target_population);
		if (cell_idx_diffuse_to>=0) {ssd_propensity_update(propensity, propensity_sums, cell_idx_diffuse_to, wildtype_populations, ra_or_ssd_populations, degradation_rate, diffusion_rate, nucleus_control_factor, density, target_population);}
	} else {
		ra_propensity_update(propensity, propensity_sums, cell_idx, wildtype_populations, ra_or_ssd_populations, degradation_rate, diffusion_rate, nucleus_control_factor, target_population, replicative_advantage);
		if (cell_idx_diffuse_to>=0) {ra_propensity_update(propensity, propensity_sums, cell_idx_diffuse_to, wildtype_populations, ra_or_ssd_populations, degradation_rate, diffusion_rate, nucleus_control_factor, target_population, replicative_advantage);}
	}
	return;
}

void ssd_gillespie_event(const gsl_rng* rng, double** propensity, double* propensity_sums, int*** wildtype_state, int*** ssd_state, int** mutant_counts, int* wildtype_populations, int* ssd_populations, long double site_std_mutation_rate, double degradation_rate, double diffusion_rate, double nucleus_control_factor, double density, int target_population) {
	/* Realise chosen event and updates system with SSDs
	
	Inputs
	------
	rng, random number generator
    propensity, 2d array of propensity vectors for each cell
    propensity_sums, 1d array of propensity sums for each cell
    wildtype_state, 3d array of wildtype state of the system
    ssd_state, 3d array of SSD state of the system
	mutant_counts, 2d array of mutant counts of each cell
	wildtype_populations, 1d array of wildtype population size of each cell
	ssd_populations, 1d array of SSD population size of each cell
	site_std_mutation_rate, standard mutation rate at each site
	degradation_rate, degradation rate
	diffusion_rate, diffusion rate
	nucleus_control_factor, nucleus control factor
	density, density
	target_population, target population */

	// Choose which cell event takes place in
	int cell_idx = weighted_sample(rng, CELLS, propensity_sums);

	// Choose event
	int event_id = weighted_sample(rng, 6, propensity[cell_idx]);

	int index_die;
	int index_born;
	int index_diffuse;
	int n_new_std_mutations;
	int cell_idx_diffuse_to = -1;
	switch (event_id) {
		case 0: // wildtype degradation
			index_die = floor(RND * wildtype_populations[cell_idx]);
			state_update_degrade(wildtype_state, cell_idx, index_die, wildtype_populations, mutant_counts);
			break;
		case 1: // wildtype replication
			index_born = floor(RND * wildtype_populations[cell_idx]);
			n_new_std_mutations = gsl_ran_binomial(rng, site_std_mutation_rate, LEN_GENOME);
			state_update_replicate(wildtype_state, cell_idx, index_born, wildtype_populations, mutant_counts, n_new_std_mutations);
			break;
		case 2: // wildtype diffusion
			// Choose cell to diffuse to
			cell_idx_diffuse_to = choose_neighbouring_cell(rng, cell_idx);
			index_diffuse = floor(RND * wildtype_populations[cell_idx]);
			state_update_diffuse(wildtype_state, cell_idx, cell_idx_diffuse_to, index_diffuse, wildtype_populations, mutant_counts);
			break;
		case 3: // SSD degradation
			index_die = floor(RND * ssd_populations[cell_idx]);
			state_update_degrade(ssd_state, cell_idx, index_die, ssd_populations, mutant_counts);
			break;
		case 4: // SSD replication
			index_born = floor(RND * ssd_populations[cell_idx]);
			n_new_std_mutations = gsl_ran_binomial(rng, site_std_mutation_rate, LEN_GENOME);
			state_update_replicate(ssd_state, cell_idx, index_born, ssd_populations, mutant_counts, n_new_std_mutations);
			break;
		case 5: // SSD diffusion
			// Choose cell to diffuse to
			cell_idx_diffuse_to = choose_neighbouring_cell(rng, cell_idx);
			index_diffuse = floor(RND * ssd_populations[cell_idx]);
			state_update_diffuse(ssd_state, cell_idx, cell_idx_diffuse_to, index_diffuse, ssd_populations, mutant_counts);
			break;
	}

	// Update propensity vectors and sum at cells where reaction occured, and with individuals diffused to
	ssd_propensity_update(propensity, propensity_sums, cell_idx, wildtype_populations, ssd_populations, degradation_rate, diffusion_rate, nucleus_control_factor, density, target_population);
	if (cell_idx_diffuse_to>=0) {ssd_propensity_update(propensity, propensity_sums, cell_idx_diffuse_to, wildtype_populations, ssd_populations, degradation_rate, diffusion_rate, nucleus_control_factor, density, target_population);}
	return;
}

int argmax(double* array, int length) {
	/* Returns the argument of the maximum of an array
	
	Inputs
	------
	array, 1d array
	length, length of array */

	int idx = 0;
	for (int i=0; i<length; ++i) {
		if (array[i] > array[idx]) {idx = i;}
	}
	return idx;
}

void get_max_std_heteroplasmies(double* max_std_heteroplasmies, int** mutant_counts, int* wildtype_populations) {
	/* Get maximum standard heteroplasmy of each cell
	
	Inputs
	------
	max_std_heteroplasmies, 1d array of maximum standard heteroplasmy of each cell
	mutant_counts, 2d array of mutant counts of each cell
	wildtype_populations, 1d array of wildtype population size of each cell
	ra_populations, 1d array of RA population size of each cell */

	int mutation_counter = mutant_counts[0][0];
	int max_mutant_count = 0;
	for (int k=0; k<CELLS; ++k) {
		for (int i=1; i<=mutation_counter; ++i){max_mutant_count = fmax(max_mutant_count, mutant_counts[k][i]);}
		max_std_heteroplasmies[k] = (double) max_mutant_count / wildtype_populations[k];
	}
	return;
}

double get_max_ra_or_ssd_heteroplasmy(int* wildtype_populations, int* ra_or_ssd_populations) {
	/* Returns the maximum RA or SSD heteroplasmy level across cells
	
	Inputs
	------
	wildtype_populations, 1d array of wildtype population size of each cell
	ra_or_ssd_populations, 1d array of RA or SSD population size of each cell */
	double heteroplasmy;
	double max_value = 0.00;
	for (int k=0; k<CELLS; ++k) {
		heteroplasmy = (double) ra_or_ssd_populations[k] / (wildtype_populations[k] + ra_or_ssd_populations[k]);
		max_value = fmax(heteroplasmy, max_value);
	}
	return max_value;
}

void introduce_ra_or_ssd(const gsl_rng* rng, int*** wildtype_state, int*** dest_state, int* wildtype_populations, int* dest_populations, int cell_idx) {
    /* Introduces RA or SSD to a wildtype individual. Also updates population size
	
	Inputs
	------
	rng, random number generator
	wildtype_state, 3d array of the wildtype state of the system
	dest_state, dynamically allocated empty 3d array of state of the system
	wildtype_populations, 1d array of wildtype population size of each cell
	dest_populations, 1d array of population size of each cell,, ie array of zeroes
	cell_idx, index of cell to which RA or SSD mutation is introduced */

	// Initialise RA/SSD state and populations of zero individuals
	for (int k=0; k<CELLS; ++k) {
		dest_state[k] = malloc(0);
		dest_populations[k] = 0;
	}

	// Allocate one row in dest_state[cell_idx] for the RA/SSD individual
	dest_state[cell_idx] = realloc(dest_state[cell_idx], sizeof(int*));

	// Copy row to dest_state
	int cell_wildtype_population = wildtype_populations[cell_idx];
	int row_idx = floor(RND * cell_wildtype_population);
	int n_mutations_to_copy = wildtype_state[cell_idx][row_idx][0];
	dest_state[cell_idx][0] = malloc((n_mutations_to_copy+1) * sizeof(int));
	for (int j=0; j<=n_mutations_to_copy; ++j) {
		dest_state[cell_idx][0][j] = wildtype_state[cell_idx][row_idx][j];
	}

	// Replace the copied row of wildtype state with its last row
	n_mutations_to_copy = wildtype_state[cell_idx][cell_wildtype_population-1][0];
	wildtype_state[cell_idx][row_idx] = realloc(wildtype_state[cell_idx][row_idx], (n_mutations_to_copy+1) * sizeof(int));
	for (int j=0; j<=n_mutations_to_copy; j++) {wildtype_state[cell_idx][row_idx][j] = wildtype_state[cell_idx][cell_wildtype_population-1][j];}
	free(wildtype_state[cell_idx][cell_wildtype_population-1]);
	wildtype_state[cell_idx] = realloc(wildtype_state[cell_idx], (cell_wildtype_population-1) * sizeof(int*));

	// Update population sizes
	wildtype_populations[cell_idx]--;
	dest_populations[cell_idx]++;
    return;
}

void compact_relabel_mutations(int** mutant_counts, int*** wildtype_state, int*** ra_or_ssd_state, int* wildtype_populations, int* ra_or_ssd_populations) {
    /* Relabels mutation identities of mutant_counts, wildtype_state and ra_or_ssd_state
	
	Inputs
	------
	mutant_counts, 2d array of mutant counts of each cell
	wildtype_state, 3d array state of the wildtype individuals of each cell
	ra_or_ssd_state, 3d array state of the RA/SSD individuals of each cell
	wildtype_populations, 1d array of wildtype population size in each cell
	ra_or_ssd_populations, 1d array of RA/SSD population size in each cell */
	
	int old_mutation_counter = mutant_counts[0][0];
	// Create mapping from old mutation id to new mutation id
    int new_id = 0;
    int* id_map = calloc(old_mutation_counter + 1, sizeof(int)); // id_map[old mutation id] = new mutation id
    for (int i=1; i<=old_mutation_counter; ++i) {
        for (int k=0; k<CELLS; ++k) {
            if (mutant_counts[k][i]) {
                id_map[i] = ++new_id;
                break;
            }
        }
    }

    // Update wildtype state using the mapping
    int n_mutations_to_copy;
    for (int k=0; k<CELLS; ++k) {
        for (int i=0; i<wildtype_populations[k]; ++i) {
            n_mutations_to_copy = wildtype_state[k][i][0];
            for (int j=1; j<=n_mutations_to_copy; ++j) {
                wildtype_state[k][i][j] = id_map[wildtype_state[k][i][j]];
            }
        }
    }

    // Update RA/SSD state using the mapping
    for (int k=0; k<CELLS; ++k) {
        for (int i=0; i<ra_or_ssd_populations[k]; ++i) {
            n_mutations_to_copy = ra_or_ssd_state[k][i][0];
            for (int j=1; j<=n_mutations_to_copy; ++j) {
                ra_or_ssd_state[k][i][j] = id_map[ra_or_ssd_state[k][i][j]];
            }
        }
    }

    // Update mutant_counts according to new mutation identity
	int new_id_count = 0;
	for (int old_id=1; old_id<=old_mutation_counter; ++old_id) {
        if (id_map[old_id]) {
            for (int k=0; k<CELLS; ++k) {
                mutant_counts[k][id_map[old_id]] = mutant_counts[k][old_id]; // since id_map[old_id] <= old_id
            }
        }
	}
    mutant_counts[0][0] = new_id;
    for (int k=0; k<CELLS; ++k) {mutant_counts[k] = realloc(mutant_counts[k], (new_id+1) * sizeof(int));}

    // Free allocated memory
    free(id_map);
    return;
}

void copy_state(int*** source, int*** dest, int* nrows) {
	/* Copies the state of the system
	
	Inputs
	------
	source, 3d array of state of the system
	dest, dynamically allocated 3d array to be copied to
	nrows, 1d array of number of rows of source state */

	int n_mutation_to_copy;
	for (int k=0; k<CELLS; ++k) {
		dest[k] = malloc(nrows[k] * sizeof(int*));
		for (int i=0; i<nrows[k]; ++i) {
			n_mutation_to_copy = source[k][i][0];
			dest[k][i] = malloc((n_mutation_to_copy+1) * sizeof(int));
			for (int j=0; j<=n_mutation_to_copy; ++j) {
				dest[k][i][j] = source[k][i][j];
			}
		}
	}
	return;
}

void copy_mutant_counts(int** source, int** dest) {
	/* Copies the mutant counts of each cell of the system

	Inputs
	------
	source, 2d array of mutant counts
	dest, dynamically allocated 2d array to be copied to */

	int source_mutation_counter = source[0][0];
	for (int k=0; k<CELLS; ++k){
		dest[k] = malloc((source_mutation_counter+1) * sizeof(int));
		for (int j=0; j<=source_mutation_counter; ++j){dest[k][j] = source[k][j];}
	}
	return;
}

void copy_population(int* source, int* dest) {
	/* Copies the population size of each cell of the system

	Inputs
	------
	source, 1d array of mutant counts
	dest, dynamically allocated 1d array to be copied to */

	for (int k=0; k<CELLS; ++k){dest[k] = source[k];}
	return;
}

void write_data_to_file(int* wildtype_populations, int* ra_or_ssd_populations, int** mutant_counts, int sim, double recording_time, FILE* fp_population, FILE* fp_sfs) {
	/* Write data to text file

	Inputs
	-------
    wildtype_populations, 1d array of population size of wildtype individuals of each cell
    ra_or_ssd_populations, 1d array of population size of RA/SSD individuals of each cell
	mutant_counts, 2d array of mutant counts of each cell
	sim, index of current simulation
	recording_time, recording timestamp
	fp_population, location of writing file for population data
	fp_sfs, location of writing file for site frequency spectrum data */

    int cell_total_population;
	int cell_mutant_count_i; // mutant count of mutations of id i in a cell
	int mutant_sum;
	int n_existing_mutation;
	int existing_mutation;

	// Compute site frequency spectrum of each cell
	double* site_frequency_spectrum;
	int bin_idx;
    for (int k=0; k<CELLS; ++k) {
        cell_total_population = wildtype_populations[k] + ra_or_ssd_populations[k];
		mutant_sum = 0;
		n_existing_mutation = 0;
		existing_mutation = 0;
		site_frequency_spectrum = calloc(N_BINS, sizeof(double));
		for (int i=1; i<=mutant_counts[0][0]; ++i) {
			cell_mutant_count_i = mutant_counts[k][i];
			if (cell_mutant_count_i) { // only count existing mutations
				mutant_sum += cell_mutant_count_i;

				// Count mutations of RA/SSD individuals with heteroplasmy level in appropriate bin
				bin_idx = (int) ((double) N_BINS * cell_mutant_count_i / cell_total_population);
				if (bin_idx==N_BINS) {bin_idx--;}
				site_frequency_spectrum[bin_idx]++;

				n_existing_mutation++;
			}

			// Normalise frequency to between 0 and 1, essentially giving the density
			for (bin_idx = 0; bin_idx<N_BINS; ++bin_idx) {
				site_frequency_spectrum[bin_idx] = site_frequency_spectrum[bin_idx] / n_existing_mutation;
			}
		}

		// Write population data to file
		fprintf(fp_population, "%d,%d,%.0lf,%d,%d\n", sim, k, recording_time, wildtype_populations[k], ra_or_ssd_populations[k]);

		// Write site frequency (density) data to file
		fprintf(fp_sfs, "%d,%d,%.0lf", sim, k, recording_time);
		for (bin_idx=0; bin_idx<N_BINS; ++bin_idx) {fprintf(fp_sfs, ",%.3f", site_frequency_spectrum[bin_idx]);}
		fprintf(fp_sfs, "\n");

		free(site_frequency_spectrum);
	}
	return;
}

int main(int argc, char *argv[]) {
	/* set up GSL RNG */
	gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);	
	/* end of GSL setup */

    int seed;
    double log_site_std_mutation_rate;
    double degradation_rate;
	double diffusion_rate;
    double nucleus_control_factor;
    double density;
    int target_population;
	double replicative_advantage;
	double target_ra_heteroplasmy;
	double target_ssd_heteroplasmy;

    long double site_std_mutation_rate;

	if(argc == 11){
		seed = atoi(argv[1]); 
		log_site_std_mutation_rate = atof(argv[2]);
		degradation_rate = atof(argv[3]);
		diffusion_rate = atof(argv[4]);
        nucleus_control_factor = atof(argv[5]);
        density = atof(argv[6]);
		target_population = atoi(argv[7]);
		replicative_advantage = atof(argv[8]);
		target_ra_heteroplasmy = atof(argv[9]);
		target_ssd_heteroplasmy = atof(argv[10]);

		gsl_rng_set(rng, seed);
	} else {
		printf("argc = %d\n", argc);
		printf("Usage: %s  seed site_std_mutation_rate(power) degradation_rate diffusion_rate nucleus_control_factor density target_population replicative_advantage target_ra_heteroplasmy target_ssd_heteroplasmy\n", argv[0]);
		return 0;	
	}

	// Set up files to write population data in
	char ra_population_filename[100];
	sprintf(ra_population_filename, "data/ra_population_sims_%d.txt", seed);
	FILE *fp_ra_population;
	fp_ra_population = fopen(ra_population_filename, "w");
	fprintf(fp_ra_population, "sim,cell,t,wildtype_population,ra_population\n");
	fclose(fp_ra_population);
	char ssd_population_filename[100];
	sprintf(ssd_population_filename, "data/ssd_population_sims_%d.txt", seed);
	FILE *fp_ssd_population;
	fp_ssd_population = fopen(ssd_population_filename, "w");
	fprintf(fp_ssd_population, "sim,cell,t,wildtype_population,ssd_population\n");
	fclose(fp_ssd_population);

	// Set up files to write site frequency spectrum data in
	char ra_sfs_filename[100];
	sprintf(ra_sfs_filename, "data/ra_site_frequency_spectrum_sims_%d.txt", seed);
	FILE *fp_ra_sfs;
	fp_ra_sfs = fopen(ra_sfs_filename, "w");
	fprintf(fp_ra_sfs, "sim,cell,t");
	double bin_lb = 0.0; // bin lower bound
	double bin_width = 1.0 / N_BINS;
	for (int bin=0; bin<N_BINS; ++bin) {
		fprintf(fp_ra_sfs, ",%.2f", bin_lb);
		bin_lb += bin_width;
	}
	fprintf(fp_ra_sfs, "\n");
	fclose(fp_ra_sfs);
	char ssd_sfs_filename[100];
	sprintf(ssd_sfs_filename, "data/ssd_site_frequency_spectrum_sims_%d.txt", seed);
	FILE *fp_ssd_sfs;
	fp_ssd_sfs = fopen(ssd_sfs_filename, "w");
	fprintf(fp_ssd_sfs, "sim,cell,t");
	bin_lb = 0.0; // bin lower bound
	bin_width = 1.0 / N_BINS;
	for (int bin=0; bin<N_BINS; ++bin) {
		fprintf(fp_ssd_sfs, ",%.2f", bin_lb);
		bin_lb += bin_width;
	}
	fprintf(fp_ssd_sfs, "\n");
	fclose(fp_ssd_sfs);

	// Generate replication rate and mutation rate for this run
	site_std_mutation_rate = pow(10, - log_site_std_mutation_rate);

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

		double* max_std_heteroplasmies = malloc(CELLS * sizeof(double));
		int max_std_mutant_count;
		int cell_with_highest_heteroplasmy;
		int n_events = 0;

		// Gillespie algorithm until a cell reaches homoplasmy
		printf("Trying to reach standard mutation homoplasmy...\n");
		do {
			// Realise event according to propensity
			non_ssd_gillespie_event(rng, wildtype_propensity, wildtype_propensity_sums, wildtype_state, wildtype_populations, mutant_counts, site_std_mutation_rate, degradation_rate, diffusion_rate, nucleus_control_factor, target_population);

			get_max_std_heteroplasmies(max_std_heteroplasmies, mutant_counts, wildtype_populations); 

			cell_with_highest_heteroplasmy = argmax(max_std_heteroplasmies, CELLS);
			n_events++;
		} while (max_std_heteroplasmies[cell_with_highest_heteroplasmy] < 1);
		printf("Standard mutation homoplasmy reached\n");

		// Free wildtype propensity memory
		for (int k=0; k<CELLS; ++k) {free(wildtype_propensity[k]);}
		free(wildtype_propensity);
		free(wildtype_propensity_sums);

		/* ra_or_ssd_state contains mutational information about RA/SSD individuals
		ra_or_ssd_state[k][i][0] is the number of standard mutations which individual i possesses in cell k
		ra_or_ssd_state[k][i][1:] are the identities of the standard mutations which individual i possesses in cell k */
		int*** ra_or_ssd_state = malloc(CELLS * sizeof(int**));
		/* ra_or_ssd_populations[k] is the population size of RA/SSD individuals in cell k */
		int* ra_or_ssd_populations = malloc(CELLS * sizeof(int));

		// Simulate systems with RA and SSD separately
		for (int ssd_sim=0; ssd_sim<=1; ++ssd_sim) {

			// Store states, populations, and mutant counts to initialise for when RA/SSD population dies out
			int*** initial_wildtype_state = malloc(CELLS * sizeof(int**));
			copy_state(wildtype_state, initial_wildtype_state, wildtype_populations);
			int* initial_wildtype_populations = malloc(CELLS * sizeof(int));
			copy_population(wildtype_populations, initial_wildtype_populations);
			int** initial_mutant_counts = malloc(CELLS * sizeof(int*));
			copy_mutant_counts(mutant_counts, initial_mutant_counts);

			// Introduce 1 RA/SSD individual to system
			introduce_ra_or_ssd(rng, wildtype_state, ra_or_ssd_state, wildtype_populations, ra_or_ssd_populations, cell_with_highest_heteroplasmy);			
			compact_relabel_mutations(mutant_counts, wildtype_state, ra_or_ssd_state, wildtype_populations, ra_or_ssd_populations);
			
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
			if (ssd_sim) {
				for (int k=0; k<CELLS; ++k) {
					propensity[k] = malloc(6 * sizeof(double));
					ssd_propensity_update(propensity, propensity_sums, k, wildtype_populations, ra_or_ssd_populations, degradation_rate, diffusion_rate, nucleus_control_factor, density, target_population);
				}
			} else {
				for (int k=0; k<CELLS; ++k) {
					propensity[k] = malloc(6 * sizeof(double));
					ra_propensity_update(propensity, propensity_sums, k, wildtype_populations, ra_or_ssd_populations, degradation_rate, diffusion_rate, nucleus_control_factor, target_population, replicative_advantage);
				}
			}

			double target_ra_or_ssd_heteroplasmy;
			if (ssd_sim) {target_ra_or_ssd_heteroplasmy = target_ssd_heteroplasmy;} else {target_ra_or_ssd_heteroplasmy = target_ra_heteroplasmy;}
			double max_ra_or_ssd_heteroplasmy = get_max_ra_or_ssd_heteroplasmy(wildtype_populations, ra_or_ssd_populations);
			int n_event2 = 0;

			if (ssd_sim) {printf("Trying to reach target SSD heteroplasmy...\n");}
			else {printf("Trying to reach target RA heteroplasmy...\n");}

			while (max_ra_or_ssd_heteroplasmy<target_ra_or_ssd_heteroplasmy) {
				printf("n_event2 = %d\n", n_event2);
				// Realise event and time of occurence according to propensity
				gillespie_event(rng, ssd_sim, propensity, propensity_sums, wildtype_state, ra_or_ssd_state, mutant_counts, wildtype_populations, ra_or_ssd_populations, site_std_mutation_rate, degradation_rate, diffusion_rate, nucleus_control_factor, density, target_population, replicative_advantage);
				
				max_ra_or_ssd_heteroplasmy = get_max_ra_or_ssd_heteroplasmy(wildtype_populations, ra_or_ssd_populations);
				n_event2++;

				// Extinction of RA/SSD individuals
				if (max_ra_or_ssd_heteroplasmy < 1e-12) {
					printf("EXTINCTION!\n");

					// Re-setup initial states
					for (int k=0; k<CELLS; ++k) {
						for (int i=0; i<wildtype_populations[k]; ++i) {free(wildtype_state[k][i]);}
						free(wildtype_state[k]);
						for (int i=0; i<ra_or_ssd_populations[k]; ++i) {free(ra_or_ssd_state[k][i]);}
						free(ra_or_ssd_state[k]);
						free(mutant_counts[k]);
					}
					copy_state(initial_wildtype_state, wildtype_state, initial_wildtype_populations);
					copy_population(initial_wildtype_populations, wildtype_populations);
					copy_mutant_counts(initial_mutant_counts, mutant_counts);

					// Re-introduce RA/SSD individual
					introduce_ra_or_ssd(rng, wildtype_state, ra_or_ssd_state, wildtype_populations, ra_or_ssd_populations, cell_with_highest_heteroplasmy);			
					compact_relabel_mutations(mutant_counts, wildtype_state, ra_or_ssd_state, wildtype_populations, ra_or_ssd_populations);
					max_ra_or_ssd_heteroplasmy = get_max_ra_or_ssd_heteroplasmy(wildtype_populations, ra_or_ssd_populations);

					// Re-setup propensity
					if (ssd_sim) {
						for (int k=0; k<CELLS; ++k) {
							ssd_propensity_update(propensity, propensity_sums, k, wildtype_populations, ra_or_ssd_populations, degradation_rate, diffusion_rate, nucleus_control_factor, density, target_population);
						}
					} else {
						for (int k=0; k<CELLS; ++k) {
							ra_propensity_update(propensity, propensity_sums, k, wildtype_populations, ra_or_ssd_populations, degradation_rate, diffusion_rate, nucleus_control_factor, target_population, replicative_advantage);
						}
					}
					n_event2 = 0;
					// printf("w = (%d %d)\n", wildtype_populations[0], wildtype_populations[1]);
					// printf("m = (%d %d)\n", ra_or_ssd_populations[0], ra_or_ssd_populations[1]);
					// for (int k=0; k<CELLS; ++k) {
					// 	printf("Cell %d:\nwildtype_state=\n", k);
					// 	print_state(wildtype_state[k], wildtype_populations[k]);
					// 	printf("ra_state=\n");
					// 	print_state(ra_or_ssd_state[k], ra_or_ssd_populations[k]);
					// }
				}
			}
			printf("Target RA/SSD heteroplasmy reached\n");
			compact_relabel_mutations(mutant_counts, wildtype_state, ra_or_ssd_state, wildtype_populations, ra_or_ssd_populations);

			// Free initial states memory
			for (int k=0; k<CELLS; ++k) {
				for (int i=0; i<initial_wildtype_populations[k]; ++i) {free(initial_wildtype_state[k][i]);}
				free(initial_wildtype_state[k]);
				free(initial_mutant_counts[k]);
			}
			free(initial_mutant_counts);
			free(initial_wildtype_state);
			printf("Initial states memory freed\n");

			double propensity_sum_across_cells;
			double current_time = 0.0;
			double recording_time = 0.0;
			int n_event3 = 0;

			// Record data at time 0
			if (current_time>(recording_time - 1E-12)){
				compact_relabel_mutations(mutant_counts, wildtype_state, ra_or_ssd_state, wildtype_populations, ra_or_ssd_populations);
				if (ssd_sim) {
					fp_ssd_population = fopen(ssd_population_filename, "a");
					fp_ssd_sfs = fopen(ssd_sfs_filename, "a");
					write_data_to_file(wildtype_populations, ra_or_ssd_populations, mutant_counts, sim, recording_time, fp_ssd_population, fp_ssd_sfs);
					fclose(fp_ssd_population);
					fclose(fp_ssd_sfs);
				} else {
					fp_ra_population = fopen(ra_population_filename, "a");
					fp_ra_sfs = fopen(ra_sfs_filename, "a");
					write_data_to_file(wildtype_populations, ra_or_ssd_populations, mutant_counts, sim, recording_time, fp_ra_population, fp_ra_sfs);
					fclose(fp_ra_population);
					fclose(fp_ra_sfs);
				}
				recording_time += RECORDING_SPACE;
			}
			
			// Simulate
			printf("Simulating...\n");
			while (current_time<SIM_LENGTH && n_event3<10000) {
				propensity_sum_across_cells = 0;
				for (int k=0; k<CELLS; ++k) {propensity_sum_across_cells += propensity_sums[k];}
				current_time += -log(RND) / propensity_sum_across_cells;
				gillespie_event(rng, ssd_sim, propensity, propensity_sums, wildtype_state, ra_or_ssd_state, mutant_counts, wildtype_populations, ra_or_ssd_populations, site_std_mutation_rate, degradation_rate, diffusion_rate, nucleus_control_factor, density, target_population, replicative_advantage);
				n_event3++;

				// In case of extinction of RA/SSD individuals, record data before beginning new simulation
				int total_ra_or_ssd_population = 0;
				for (int k=0; k<CELLS; ++k) {total_ra_or_ssd_population += ra_or_ssd_populations[k];}
				if (total_ra_or_ssd_population==0) {
					compact_relabel_mutations(mutant_counts, wildtype_state, ra_or_ssd_state, wildtype_populations, ra_or_ssd_populations);
					if (ssd_sim) {
						fp_ssd_population = fopen(ssd_population_filename, "a");
						fp_ssd_sfs = fopen(ssd_sfs_filename, "a");
						write_data_to_file(wildtype_populations, ra_or_ssd_populations, mutant_counts, sim, recording_time, fp_ssd_population, fp_ssd_sfs);
						fclose(fp_ssd_population);
						fclose(fp_ssd_sfs);
						printf("EXTINCTION!\n");
					} else {
						fp_ra_population = fopen(ra_population_filename, "a");
						fp_ra_sfs = fopen(ra_sfs_filename, "a");
						write_data_to_file(wildtype_populations, ra_or_ssd_populations, mutant_counts, sim, recording_time, fp_ra_population, fp_ra_sfs);
						fclose(fp_ra_population);
						fclose(fp_ra_sfs);
						printf("EXTINCTION!\n");
					}
				}

				// Record data at recording time
				if (current_time>(recording_time - 1E-12)){
					compact_relabel_mutations(mutant_counts, wildtype_state, ra_or_ssd_state, wildtype_populations, ra_or_ssd_populations);
					if (ssd_sim) {
						fp_ssd_population = fopen(ssd_population_filename, "a");
						fp_ssd_sfs = fopen(ssd_sfs_filename, "a");
						write_data_to_file(wildtype_populations, ra_or_ssd_populations, mutant_counts, sim, recording_time, fp_ssd_population, fp_ssd_sfs);
						fclose(fp_ssd_population);
						fclose(fp_ssd_sfs);
					} else {
						fp_ra_population = fopen(ra_population_filename, "a");
						fp_ra_sfs = fopen(ra_sfs_filename, "a");
						write_data_to_file(wildtype_populations, ra_or_ssd_populations, mutant_counts, sim, recording_time, fp_ra_population, fp_ra_sfs);
						fclose(fp_ra_population);
						fclose(fp_ra_sfs);
					}
					recording_time += RECORDING_SPACE;
				}
			}
			// Free memory
			for (int k=0; k<CELLS; ++k){
				free(propensity[k]);

				for (int i=0; i<wildtype_populations[k]; ++i){free(wildtype_state[k][i]);}
				free(wildtype_state[k]);
				for (int i=0; i<ra_or_ssd_populations[k]; ++i){free(ra_or_ssd_state[k][i]);}
				free(ra_or_ssd_state[k]);

				free(mutant_counts[k]);
			}
		}
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