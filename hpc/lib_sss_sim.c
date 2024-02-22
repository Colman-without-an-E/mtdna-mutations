#include "../include/parameters_hpc.h"
#include "../include/lib_sss_sim.h"

#define RND gsl_rng_uniform_pos(rng) // generate number from Unif(0,1)

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

	// Absorbing boundary case
	if (cell_idx_to==-1 || cell_idx_to==CELLS) {
		state_update_degrade(state, cell_idx_from, diffuse_row, populations, mutant_counts);
		return;
	}

	// General non-boundary case
	else {
		// Reallocate for new row of state to be diffused into
		int n_mutations_to_copy = state[cell_idx_from][diffuse_row][0];
		state[cell_idx_to] = realloc(state[cell_idx_to], (cell_population_to + 1) * sizeof(int*));
		state[cell_idx_to][cell_population_to] = malloc((n_mutations_to_copy+1) * sizeof(int));

		// Copy content to the newly allocated row
		state[cell_idx_to][cell_population_to][0] = n_mutations_to_copy;
		for (int j=1; j<=n_mutations_to_copy; ++j) {
			int mutation_id = state[cell_idx_from][diffuse_row][j];
			state[cell_idx_to][cell_population_to][j] = mutation_id;
			if (mutation_id) {mutant_counts[cell_idx_to][mutation_id]++;}
		}
		populations[cell_idx_to]++;

		// Remove original row
		state_update_degrade(state, cell_idx_from, diffuse_row, populations, mutant_counts);
		return;
	}
}

int choose_neighbouring_cell(const gsl_rng* rng, int cell_idx) {
	/* Returns index of neighbouring cell with equal probability
	
	Inputs
	------
	rng, rng
	cell_idx, index of cell whose neighbouring cell index is to be found */

	int neighbour_cell_idx = cell_idx;
	(gsl_rng_uniform_int(rng, 2)) ? neighbour_cell_idx++ : neighbour_cell_idx--;
	return neighbour_cell_idx;
}

void wildtype_propensity_update(double** propensity, double* propensity_sums, int cell_idx, int* wildtype_populations, double degradation_rate, double diffusion_rate, double nucleus_control_factor, int target_population, double rate_difference) {
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
    target_population, target population
	rate_difference, rate difference */

	int cell_wildtype_population = wildtype_populations[cell_idx];
	double wildtype_replication_rate = degradation_rate + nucleus_control_factor * (target_population - cell_wildtype_population);

	propensity[cell_idx][0] = (degradation_rate + rate_difference) * cell_wildtype_population;
	propensity[cell_idx][1] = fmax(0.0, (wildtype_replication_rate + rate_difference) * cell_wildtype_population);
	propensity[cell_idx][2] = diffusion_rate * cell_wildtype_population;

	propensity_sums[cell_idx] = 0.00;
	for (int i=0; i<3; ++i){propensity_sums[cell_idx] += propensity[cell_idx][i];}
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
    
    double random_point = gsl_rng_uniform_pos(rng) * weight_sum;
    double cum_sum=0;
    
    for (int i=0; i<n-1; ++i) {
        cum_sum += weights[i];
        if (cum_sum > random_point) {
            return i;
        }
    }
    return n-1;
}

void wildtype_gillespie_event(const gsl_rng* rng, double** propensity, double* propensity_sums, int*** wildtype_state, int* wildtype_populations, int** mutant_counts, long double site_std_mutation_rate, double degradation_rate, double diffusion_rate, double nucleus_control_factor, int target_population, double rate_difference) {
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
	wildtype_propensity_update(propensity, propensity_sums, cell_idx, wildtype_populations, degradation_rate, diffusion_rate, nucleus_control_factor, target_population, rate_difference);
	if (0<=cell_idx_diffuse_to && cell_idx_diffuse_to<CELLS) {wildtype_propensity_update(propensity, propensity_sums, cell_idx_diffuse_to, wildtype_populations, degradation_rate, diffusion_rate, nucleus_control_factor, target_population, rate_difference);}
	return;
}

void gillespie_event(const gsl_rng* rng, double** propensity, double* propensity_sums, int*** wildtype_state, int*** ra_or_ssd_state, int** mutant_counts, int* wildtype_populations, int* ra_or_ssd_populations, long double site_std_mutation_rate, double degradation_rate, double diffusion_rate, double nucleus_control_factor, int target_population, double rate_difference) {
	/* Realise chosen event and updates system with RAs
	
	Inputs
	------
	rng, random number generator
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
	target_population, target population
	rate_difference, rate difference */

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
		case 3: // RA/SSD degradation
			index_die = floor(RND * ra_or_ssd_populations[cell_idx]);
			state_update_degrade(ra_or_ssd_state, cell_idx, index_die, ra_or_ssd_populations, mutant_counts);
			break;
		case 4: // RA/SSD replication
			index_born = floor(RND * ra_or_ssd_populations[cell_idx]);
			n_new_std_mutations = gsl_ran_binomial(rng, site_std_mutation_rate, LEN_GENOME);
			state_update_replicate(ra_or_ssd_state, cell_idx, index_born, ra_or_ssd_populations, mutant_counts, n_new_std_mutations);
			break;
		case 5: // RA/SSD diffusion
			// Choose cell to diffuse to
			cell_idx_diffuse_to = choose_neighbouring_cell(rng, cell_idx);
			index_diffuse = floor(RND * ra_or_ssd_populations[cell_idx]);
			state_update_diffuse(ra_or_ssd_state, cell_idx, cell_idx_diffuse_to, index_diffuse, ra_or_ssd_populations, mutant_counts);
			break;
	}

	// Update propensity vectors and sum at cells where reaction occured, and with individuals diffused to
	propensity_update(propensity, propensity_sums, cell_idx, wildtype_populations, ra_or_ssd_populations, degradation_rate, diffusion_rate, nucleus_control_factor, target_population, rate_difference);
	if (0<=cell_idx_diffuse_to && cell_idx_diffuse_to<CELLS) {propensity_update(propensity, propensity_sums, cell_idx_diffuse_to, wildtype_populations, ra_or_ssd_populations, degradation_rate, diffusion_rate, nucleus_control_factor, target_population, rate_difference);}
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
	ra_or_ssd_populations, 1d array of RA population size of each cell */

	int mutation_counter = mutant_counts[0][0];
	int max_mutant_count = 0;
	for (int k=0; k<CELLS; ++k) {
		for (int i=1; i<=mutation_counter; ++i){max_mutant_count = fmax(max_mutant_count, mutant_counts[k][i]);}
		max_std_heteroplasmies[k] = (double) max_mutant_count / wildtype_populations[k];
	}
	return;
}

double get_max_ra_or_ssd_heteroplasmy(int* wildtype_populations, int* ra_or_ssd_populations) {
	/* Returns the maximum RA/SSD heteroplasmy level across cells
	
	Inputs
	------
	wildtype_populations, 1d array of wildtype population size of each cell
	ra_or_ssd_populations, 1d array of RA/SSD population size of each cell */
	double heteroplasmy;
	double max_value = 0.00;
	for (int k=0; k<CELLS; ++k) {
		heteroplasmy = (double) ra_or_ssd_populations[k] / (wildtype_populations[k] + ra_or_ssd_populations[k]);
		max_value = fmax(heteroplasmy, max_value);
	}
	return max_value;
}

void introduce_ra_or_ssd(const gsl_rng* rng, int*** wildtype_state, int*** dest_state, int* wildtype_populations, int* dest_populations, int cell_idx) {
    /* Introduces RA to a wildtype individual. Also updates population size
	
	Inputs
	------
	rng, random number generator
	wildtype_state, 3d array of the wildtype state of the system
	dest_state, dynamically allocated empty 3d array of state of the system
	wildtype_populations, 1d array of wildtype population size of each cell
	dest_populations, 1d array of population size of each cell,, ie array of zeroes
	cell_idx, index of cell to which RA mutation is introduced */

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

void compact_relabel_wildtype_mutations(int** mutant_counts, int*** wildtype_state, int* wildtype_populations) {
    /* Relabels mutation identities of mutant_counts and wildtype_state
	
	Inputs
	------
	mutant_counts, 2d array of mutant counts of each cell
	wildtype_state, 3d array state of the wildtype individuals of each cell
	wildtype_populations, 1d array of wildtype population size in each cell */
	
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

void write_data_to_file_pre_introduce(int* wildtype_populations, int** mutant_counts, int sim, double recording_time, char* population_filename, char* sfs_filename) {
	/* Write data pre-introduction of RA/SSD to text file

	Inputs
	-------
    wildtype_populations, 1d array of population size of wildtype individuals of each cell
	mutant_counts, 2d array of mutant counts of each cell
	sim, index of current simulation
	recording_time, recording timestamp
	fp_population, location of writing file for population data
	fp_sfs, location of writing file for site frequency spectrum data */

	// Write population data to file
	FILE* fp_population = fopen(population_filename, "a");
    for (int k=0; k<CELLS; ++k) {
		fprintf(fp_population, "%d,%d,%.0lf,%d,%d\n", sim, k, recording_time, wildtype_populations[k], 0);
	}
	fclose(fp_population);

	// Write site frequency data to file
	FILE* fp_sfs = fopen(sfs_filename, "a");
	fprintf(fp_sfs, "%d,%.0lf", sim, recording_time);
	for (int mut_id=1; mut_id<=mutant_counts[0][0]; ++mut_id) {
		int mutant_count_across_cells = 0;
		for (int k=0; k<CELLS; ++k) {
			mutant_count_across_cells += mutant_counts[k][mut_id];
		}
		fprintf(fp_sfs, ",%d", mutant_count_across_cells);
	}
	fprintf(fp_sfs, "\n");
	fclose(fp_sfs);
	return;
}

void write_data_to_file(int* wildtype_populations, int* ra_or_ssd_populations, int** mutant_counts, int sim, double recording_time, char* population_filename, char* sfs_filename) {
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

	// Write population data to file
	FILE* fp_population = fopen(population_filename, "a");
    for (int k=0; k<CELLS; ++k) {
		fprintf(fp_population, "%d,%d,%.0lf,%d,%d\n", sim, k, recording_time, wildtype_populations[k], ra_or_ssd_populations[k]);
	}
	fclose(fp_population);

	// Write site frequency data to file
	FILE* fp_sfs = fopen(sfs_filename, "a");
	fprintf(fp_sfs, "%d,%.0lf", sim, recording_time);
	for (int mut_id=1; mut_id<mutant_counts[0][0]; ++mut_id) {
		int mutant_count_across_cells = 0;
		for (int k=0; k<CELLS; ++k) {
			mutant_count_across_cells += mutant_counts[k][mut_id];
		}
		fprintf(fp_sfs, ",%d", mutant_count_across_cells);
	}
	fprintf(fp_sfs, "\n");
	fclose(fp_sfs);
	return;
}

// void write_data_to_file_pre_introduce(int* wildtype_populations, int** mutant_counts, int sim, double recording_time, char* population_filename, char* sfs_filename) {
// 	/* Write data pre-introduction of RA/SSD to text file

// 	Inputs
// 	-------
//     wildtype_populations, 1d array of population size of wildtype individuals of each cell
// 	mutant_counts, 2d array of mutant counts of each cell
// 	sim, index of current simulation
// 	recording_time, recording timestamp
// 	fp_population, location of writing file for population data
// 	fp_sfs, location of writing file for site frequency spectrum data */

//     int cell_total_population;
// 	int cell_mutant_count_i; // mutant count of mutations of id i in a cell
// 	int mutant_sum;
// 	int n_existing_mutation;
// 	int existing_mutation;

// 	// Compute site frequency spectrum of each cell
// 	int bin_idx;
//     for (int k=0; k<CELLS; ++k) {
// 		mutant_sum = 0;
// 		n_existing_mutation = 0;
// 		existing_mutation = 0;
// 		int site_frequency_spectrum[N_BINS] = {0};
// 		for (int i=1; i<=mutant_counts[0][0]; ++i) {
// 			cell_mutant_count_i = mutant_counts[k][i];
// 			if (cell_mutant_count_i) { // only count existing mutations
// 				mutant_sum += cell_mutant_count_i;

// 				// Count mutations of RA/SSD individuals with heteroplasmy level in appropriate bin
// 				bin_idx = (int) ((double) N_BINS * cell_mutant_count_i / wildtype_populations[k]);
// 				if (bin_idx==N_BINS) {bin_idx--;}
// 				site_frequency_spectrum[bin_idx]++;

// 				n_existing_mutation++;
// 			}
// 		}

// 		// Write population data to file
// 		FILE* fp_population = fopen(population_filename, "a");
// 		fprintf(fp_population, "%d,%d,%.0lf,%d,%d\n", sim, k, recording_time, wildtype_populations[k], 0);
//         fclose(fp_population);

// 		// Write site frequency (density) data to file
// 		FILE* fp_sfs = fopen(sfs_filename, "a");
// 		fprintf(fp_sfs, "%d,%d,%.0lf", sim, k, recording_time);
// 		for (bin_idx=0; bin_idx<N_BINS; ++bin_idx) {fprintf(fp_sfs, ",%d", site_frequency_spectrum[bin_idx]);}
// 		fprintf(fp_sfs, "\n");
// 		fclose(fp_sfs);
// 	}
// 	return;
// }

// void write_data_to_file(int* wildtype_populations, int* ra_or_ssd_populations, int** mutant_counts, int sim, double recording_time, char* population_filename, char* sfs_filename) {
// 	/* Write data to text file

// 	Inputs
// 	-------
//     wildtype_populations, 1d array of population size of wildtype individuals of each cell
//     ra_or_ssd_populations, 1d array of population size of RA/SSD individuals of each cell
// 	mutant_counts, 2d array of mutant counts of each cell
// 	sim, index of current simulation
// 	recording_time, recording timestamp
// 	fp_population, location of writing file for population data
// 	fp_sfs, location of writing file for site frequency spectrum data */

//     int cell_total_population;
// 	int cell_mutant_count_i; // mutant count of mutations of id i in a cell
// 	int mutant_sum;
// 	int n_existing_mutation;
// 	int existing_mutation;

// 	// Compute site frequency spectrum of each cell
// 	int bin_idx;
//     for (int k=0; k<CELLS; ++k) {
//         cell_total_population = wildtype_populations[k] + ra_or_ssd_populations[k];
// 		mutant_sum = 0;
// 		n_existing_mutation = 0;
// 		existing_mutation = 0;
// 		int site_frequency_spectrum[N_BINS] = {0};
// 		for (int i=1; i<=mutant_counts[0][0]; ++i) {
// 			cell_mutant_count_i = mutant_counts[k][i];
// 			if (cell_mutant_count_i) { // only count existing mutations
// 				mutant_sum += cell_mutant_count_i;

// 				// Count mutations of RA/SSD individuals with heteroplasmy level in appropriate bin
// 				bin_idx = (int) ((double) N_BINS * cell_mutant_count_i / cell_total_population);
// 				if (bin_idx==N_BINS) {bin_idx--;}
// 				site_frequency_spectrum[bin_idx]++;

// 				n_existing_mutation++;
// 			}
// 		}

// 		// Write population data to file
// 		FILE* fp_population = fopen(population_filename, "a");
// 		fprintf(fp_population, "%d,%d,%.0lf,%d,%d\n", sim, k, recording_time, wildtype_populations[k], ra_or_ssd_populations[k]);
//         fclose(fp_population);

// 		// Write site frequency (density) data to file
// 		FILE* fp_sfs = fopen(sfs_filename, "a");
// 		fprintf(fp_sfs, "%d,%d,%.0lf", sim, k, recording_time);
// 		for (bin_idx=0; bin_idx<N_BINS; ++bin_idx) {fprintf(fp_sfs, ",%d", site_frequency_spectrum[bin_idx]);}
// 		fprintf(fp_sfs, "\n");
// 		fclose(fp_sfs);
// 	}
// 	return;
// }
