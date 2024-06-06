#include "../include/parameters_hpc.h"
#include "../include/lib_sim.h"

#define HASHSIZE 211
#define RND gsl_rng_uniform_pos(rng) // generate number from Unif(0,1)

/* Counter: Build a dictionary structure and counting function for succinctly recording distributional data */

/* KeyValuePair.key is an integer key
KeyValuePair.value is an integer value
KeyValuePair.next is a pointer to the next KeyValuePair whose key has the same hash value */
struct KeyValuePair {
    int key;
    int value;
    struct KeyValuePair *next;
};

/* Dict.hash_table[i] is a pointer to a KeyValuePair whose hashed key is i */
typedef struct {
    unsigned hash_size;
    struct KeyValuePair** hash_table;
} Dict;

unsigned hash(int k, unsigned hash_size) {
    /* Simple hash function. Modular arithmetic. */
    return k % hash_size;
}

struct KeyValuePair* lookup(Dict dict, int key) {
    /* Performs hash table lookup of a dictionary
    
    Inputs
    ------
    dict, pointer to dictionary
    key, lookup key
    
    Returns
    -------
    pointer to the required key-value pair */

    struct KeyValuePair* item;
    for (item=dict.hash_table[hash(key, dict.hash_size)]; item!=NULL; item=item->next) {
        if (key==item->key) {
            return item; // found
        }
    }
    return NULL; // not found
}

void print_dictionary(Dict dict) {
    /* Prints all key-value pairs of a dictionary in the format of key:value.
    
    Inputs
    ------
    dict, dictionary */
    for (int i=0; i<dict.hash_size; ++i) {
        struct KeyValuePair* item;
        for (item=dict.hash_table[i]; item!=NULL; item=item->next) {
            printf("%d:%d\n", item->key, item->value);
        }
    }
}

Dict counter(int *array, int length, int hash_size) {
    /* Returns a dictionary which counts the occurences of unique integers in array
    
    Inputs
    ------
    array, 1d array of ints
    length, length of array */

    // Initialize dictionary with hash_table pointers to NULL
    Dict dict;
    dict.hash_size = hash_size;
    dict.hash_table = malloc(dict.hash_size * sizeof(struct KeyValuePair*));
    for (int h=0; h<dict.hash_size; ++h) {dict.hash_table[h] = NULL;}

    for (int i=0; i<length; ++i) {
        struct KeyValuePair* item = lookup(dict, array[i]);

        // Install new key-value pair to dictionary if array[i] not found
        if (item == NULL) {
            unsigned hash_value = hash(array[i], dict.hash_size);
            item = malloc(sizeof(*item));
            item->key = array[i];
            item->value = 1;
            item->next = dict.hash_table[hash_value];
            dict.hash_table[hash_value] = item;
        } else {
            item->value++; // increment count by 1 when array[i] is found
        }
    }
    return dict;
}

/* End of Counter */


/* Functions used in Gillespie algorithm */

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

	// Closed boundary case
	if (cell_idx_to==-1 || cell_idx_to==CELLS) {
		// state_update_degrade(state, cell_idx_from, diffuse_row, populations, mutant_counts);
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
	propensity[cell_idx][1] = fmax(0.0, wildtype_replication_rate * cell_wildtype_population);
	propensity[cell_idx][2] = 2 * diffusion_rate * cell_wildtype_population;

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

void wildtype_gillespie_event(const gsl_rng* rng, double** propensity, double* propensity_sums, int*** wildtype_state, int* wildtype_populations, int** mutant_counts, long double site_std_mutation_rate, double degradation_rate, double diffusion_rate, double nucleus_control_factor, int target_population) {
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
	if (0<=cell_idx_diffuse_to && cell_idx_diffuse_to<CELLS) {wildtype_propensity_update(propensity, propensity_sums, cell_idx_diffuse_to, wildtype_populations, degradation_rate, diffusion_rate, nucleus_control_factor, target_population);}
	return;
}

void gillespie_event(const gsl_rng* rng, double** propensity, double* propensity_sums, int*** wildtype_state, int*** ra_or_ssd_state, int** mutant_counts, int* wildtype_populations, int* ra_or_ssd_populations, long double site_std_mutation_rate, double degradation_rate, double diffusion_rate, double nucleus_control_factor, int target_population, double replicative_advantage) {
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
	propensity_update(propensity, propensity_sums, cell_idx, wildtype_populations, ra_or_ssd_populations, degradation_rate, diffusion_rate, nucleus_control_factor, target_population, replicative_advantage);
	if (0<=cell_idx_diffuse_to && cell_idx_diffuse_to<CELLS) {propensity_update(propensity, propensity_sums, cell_idx_diffuse_to, wildtype_populations, ra_or_ssd_populations, degradation_rate, diffusion_rate, nucleus_control_factor, target_population, replicative_advantage);}
	return;
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

/* End of Gillespie algorithm functions */


/* Data storage functions */

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

void write_data_to_file_pre_introduce(int* wildtype_populations, int** mutant_counts, double recording_time, char* population_filename, char* sfs_filename) {
	/* Write data pre-introduction of RA/SSD to text file

	Inputs
	-------
    wildtype_populations, 1d array of population size of wildtype individuals of each cell
	mutant_counts, 2d array of mutant counts of each cell
	recording_time, recording timestamp
	fp_population, location of writing file for population data
	fp_sfs, location of writing file for site frequency spectrum data */

	// Write aggregated population data to file
	FILE* fp_population = fopen(population_filename, "a");
	int w_tot = 0;
    for (int k=0; k<CELLS; ++k) {
		w_tot += wildtype_populations[k];
	}
	fprintf(fp_population, "%.0lf,%d,%d,%d\n", recording_time, -1, w_tot, 0);

	// Write cell-wise population data to file
    for (int k=0; k<CELLS; ++k) {
		fprintf(fp_population, "%.0lf,%d,%d,%d\n", recording_time, k, wildtype_populations[k], 0);
	}
	fclose(fp_population);

	// Write aggregated site frequency spectrum data across cells to file, encoded as cell index = -1
	FILE* fp_sfs = fopen(sfs_filename, "a");
	int* tot_mutant_counts = calloc(mutant_counts[0][0], sizeof(int));
	for (int mut_id=1; mut_id<=mutant_counts[0][0]; ++mut_id) {
		for (int k=0; k<CELLS; ++k) {
			tot_mutant_counts[mut_id-1] += mutant_counts[k][mut_id];
		}
	}
	Dict tot_counts = counter(tot_mutant_counts, mutant_counts[0][0], HASHSIZE);
	fprintf(fp_sfs, "%.0lf,%d,{", recording_time, -1);
	int comma_flag = 0;
	for (int h=0; h<tot_counts.hash_size; ++h) {
		for (struct KeyValuePair* item=tot_counts.hash_table[h]; item!=NULL; item=item->next) {
			if (comma_flag) {
				fprintf(fp_sfs, ",\"%d\":%d", item->key, item->value);
			} else {
				fprintf(fp_sfs, "\"%d\":%d", item->key, item->value);
			}
			comma_flag = 1;
		}
	}
	fprintf(fp_sfs, "}\n");
	free(tot_mutant_counts);
	
	// Write cell-wise site frequency spectrum data to file
	int* cell_mutant_counts = malloc(mutant_counts[0][0]*sizeof(int));
	for (int k=0; k<CELLS; ++k) {
		for (int mut_id=1; mut_id<=mutant_counts[0][0]; ++mut_id) {
			cell_mutant_counts[mut_id-1] = mutant_counts[k][mut_id];
		}
		Dict counts = counter(cell_mutant_counts, mutant_counts[0][0], HASHSIZE);
		fprintf(fp_sfs, "%.0lf,%d,{", recording_time, k);
		int comma_flag = 0;
		for (int h=0; h<counts.hash_size; ++h) {
			for (struct KeyValuePair* item=counts.hash_table[h]; item!=NULL; item=item->next) {
				if (comma_flag) {
					fprintf(fp_sfs, ",\"%d\":%d", item->key, item->value);
				} else {
					fprintf(fp_sfs, "\"%d\":%d", item->key, item->value);
				}
				comma_flag = 1;
			}
		}
		fprintf(fp_sfs, "}\n");
	}
	fclose(fp_sfs);
	free(cell_mutant_counts);
	return;
}

void write_data_to_file(int* wildtype_populations, int* ra_or_ssd_populations, int** mutant_counts, double recording_time, char* population_filename, char* sfs_filename) {
	/* Write data to text file

	Inputs
	-------
    wildtype_populations, 1d array of population size of wildtype individuals of each cell
    ra_or_ssd_populations, 1d array of population size of RA/SSD individuals of each cell
	mutant_counts, 2d array of mutant counts of each cell
	recording_time, recording timestamp
	fp_population, location of writing file for population data
	fp_sfs, location of writing file for site frequency spectrum data */

	// Write aggregated population data to file
	FILE* fp_population = fopen(population_filename, "a");
	int w_tot = 0;
	int m_tot = 0;
    for (int k=0; k<CELLS; ++k) {
		w_tot += wildtype_populations[k];
		m_tot += ra_or_ssd_populations[k];
	}
	fprintf(fp_population, "%.0lf,%d,%d,%d\n", recording_time, -1, w_tot, m_tot);

	// Write population data to file
    for (int k=0; k<CELLS; ++k) {
		fprintf(fp_population, "%.0lf,%d,%d,%d\n", recording_time, k, wildtype_populations[k], ra_or_ssd_populations[k]);
	}
	fclose(fp_population);

	// Write aggregated site frequency spectrum data across cells to file, encoded as cell index = -1
	FILE* fp_sfs = fopen(sfs_filename, "a");
	int* tot_mutant_counts = calloc(mutant_counts[0][0], sizeof(int));
	for (int mut_id=1; mut_id<=mutant_counts[0][0]; ++mut_id) {
		for (int k=0; k<CELLS; ++k) {
			tot_mutant_counts[mut_id-1] += mutant_counts[k][mut_id];
		}
	}
	Dict tot_counts = counter(tot_mutant_counts, mutant_counts[0][0], HASHSIZE);
	fprintf(fp_sfs, "%.0lf,%d,{", recording_time, -1);
	int comma_flag = 0;
	for (int h=0; h<tot_counts.hash_size; ++h) {
		for (struct KeyValuePair* item=tot_counts.hash_table[h]; item!=NULL; item=item->next) {
			if (comma_flag) {
				fprintf(fp_sfs, ",\"%d\":%d", item->key, item->value);
			} else {
				fprintf(fp_sfs, "\"%d\":%d", item->key, item->value);
			}
			comma_flag = 1;
		}
	}
	fprintf(fp_sfs, "}\n");
	free(tot_mutant_counts);
	
	// Write cell-wise site frequency spectrum data to file
	int* cell_mutant_counts = malloc(mutant_counts[0][0]*sizeof(int));
	for (int k=0; k<CELLS; ++k) {
		for (int mut_id=1; mut_id<=mutant_counts[0][0]; ++mut_id) {
			cell_mutant_counts[mut_id-1] = mutant_counts[k][mut_id];
		}
		Dict counts = counter(cell_mutant_counts, mutant_counts[0][0], HASHSIZE);
		fprintf(fp_sfs, "%.0lf,%d,{", recording_time, k);
		int comma_flag = 0;
		for (int h=0; h<counts.hash_size; ++h) {
			for (struct KeyValuePair* item=counts.hash_table[h]; item!=NULL; item=item->next) {
				if (comma_flag) {
					fprintf(fp_sfs, ",\"%d\":%d", item->key, item->value);
				} else {
					fprintf(fp_sfs, "\"%d\":%d", item->key, item->value);
				}
				comma_flag = 1;
			}
		}
		fprintf(fp_sfs, "}\n");
	}
	fclose(fp_sfs);
	free(cell_mutant_counts);
	return;
}

/* End of data storage */


/* Debug functions */

void print_state(int*** state, int* populations, int cell_idx) {
	/* Print the state of the system for debugging
	
	Inputs
	------
	state, 3d array of the state
	populations, population size of each cell
    cell_idx, index of cell of which to print the state */

    int nrow_to_print = populations[cell_idx];
	int n_mutations_to_print;
	for (int i=0; i<nrow_to_print; ++i) {
		n_mutations_to_print = state[cell_idx][i][0];
		for (int j=0; j<=n_mutations_to_print; ++j) {
			printf("%d ", state[cell_idx][i][j]);
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

void generate_state(const gsl_rng* rng, int*** state, int* populations, int pop_mu, double mut_p, int mut_id_mu) {
    /* Generates state randomly
    
    Inputs
    ------
    rng, rng
    state, empty 3d array to store state
    populations, empty 1d array to store population size of each cell
    pop_mu, mean population for generation
    mut_p, probability of a mutation of each individual for generation
    mut_id_mu, mean mutation id (mutation counter) for generation */
    
    int cell_population;
    int mutation_counter;
    int n_mutations;
    // int mut_id;

    // Generate non-zero overall-maximum mutation id, mutation_counter
    do {
        mutation_counter = gsl_ran_poisson(rng, mut_id_mu);
    } while (mutation_counter==0);
    int* mut_ids_asc = malloc(mutation_counter * sizeof(int));
    for (int i=0; i<mutation_counter; ++i) {mut_ids_asc[i]=i+1;}

    for (int k=0; k<CELLS; ++k) {

        // Generate non-zero population size 
        do {
            cell_population = gsl_ran_poisson(rng, pop_mu);
        } while (cell_population==0);
        populations[k] = cell_population;

        // Allocate memory for individuals
        state[k] = malloc(cell_population * sizeof(int*));

        for (int mut_id=1; mut_id<=mutation_counter; ++mut_id) {
            
            // Generate number of mutations of id mut_id for each individual
            for (int i=0; i<cell_population; ++i) {
                n_mutations = gsl_ran_binomial(rng, mut_p, mutation_counter);
                state[k][i] = malloc((n_mutations+1) * sizeof(int));
                state[k][i][0] = n_mutations;
                int* mut_state = malloc(n_mutations * sizeof(int));
                gsl_ran_choose(rng, mut_state, n_mutations, mut_ids_asc, mutation_counter, sizeof(int));
                for (int j=0; j<n_mutations; ++j) {state[k][i][j+1] = mut_state[j];}
                free(mut_state);
            }
        }
    }
    return;
}

void get_mutant_counts_wildtype(int** mutant_counts, int*** wildtype_state, int* wildtype_population) {
    /* Stores mutant counts of wildtype system in the 2D array mutant_counts */
    
	int mut_id;
    int n_mut;
    int mutation_counter = 0;
    for (int k=0; k<CELLS; ++k) {
        for (int i=0; i<wildtype_population[k]; ++i) {
            n_mut = wildtype_state[k][i][0];
            for (int j=1; j<=n_mut; ++j) {
                mut_id = wildtype_state[k][i][j];
                mutant_counts[k][mut_id]++;
                mutation_counter = fmax(mutation_counter, mut_id);
            }
        }
    }

    mutant_counts[0][0] = mutation_counter;
    for (int k=0; k<CELLS; ++k) {
        mutant_counts[k] = realloc(mutant_counts[k], (mutation_counter+1) * sizeof(int));
    }
    return;
}

void get_mutant_counts(int** mutant_counts, int*** wildtype_state, int*** ra_state, int* wildtype_population, int* ra_population) {
    /* Stores mutant counts of system in the 2D array mutant_counts */
	
	int mut_id;
    int n_mut;
    int mutation_counter = 0;
    for (int k=0; k<CELLS; ++k) {
        for (int i=0; i<wildtype_population[k]; ++i) {
            n_mut = wildtype_state[k][i][0];
            for (int j=1; j<=n_mut; ++j) {
                mut_id = wildtype_state[k][i][j];
                mutant_counts[k][mut_id]++;
                mutation_counter = fmax(mutation_counter, mut_id);
            }
        }
        for (int i=0; i<ra_population[k]; ++i) {
            n_mut = ra_state[k][i][0];
            for (int j=1; j<=n_mut; ++j) {
                mut_id = ra_state[k][i][j];
                mutant_counts[k][mut_id]++;
                mutation_counter = fmax(mutation_counter, mut_id);
            }
        }
    }

    mutant_counts[0][0] = mutation_counter;
    for (int k=0; k<CELLS; ++k) {
        mutant_counts[k] = realloc(mutant_counts[k], (mutation_counter+1) * sizeof(int));
    }
    return;
}

int identical_mutant_counts(int** mutant_counts1, int** mutant_counts2) {
	/* Compares mutant counts. Returns 1 if they are identical, 0 otherwise. */

	int mutation_counter = mutant_counts1[0][0];
	if (mutant_counts2[0][0]!=mutation_counter) {
		return 0;
	}

	for (int k=0; k<CELLS; ++k) {
		for (int j=1; j<=mutation_counter; ++j) {
			if (mutant_counts1[k][j]!=mutant_counts2[k][j]) {
				return 0;
			}
		}
	}
	return 1;
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

int argmax_mutant_counts(int* args, int** mutant_counts) {
	/* Returns the maximum of mutant counts. Also stores the corresponding arguments. 
	
	Inputs
	------
	args, 1d array of 2 entries to store indices of argmax
	mutant_counts, 2d array of mutant counts */

	// Initialise max and argmax
	int max = 0;
	args[0] = 0;
	args[1] = 0;

	// Iterate over every element in mutant_counts and update max and argmax
	for (int k=0; k<CELLS; ++k) {
		for (int j=1; j<=mutant_counts[0][0]; ++j) {
			if (mutant_counts[k][j]>max) {
				max = mutant_counts[k][j];
				args[0] = k;
				args[1] = j;
			}
		}
	}
	return max;
}

/* End of debug function */