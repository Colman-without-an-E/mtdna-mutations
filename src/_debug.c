#include "../include/_debug.h"

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

void get_mutant_counts_wildtype(int** mutant_counts, int*** wildtype_state, int* wildtype_populations) {
    int mut_id;
    int n_mut;
    int mutation_counter = 0;
    for (int k=0; k<CELLS; ++k) {
        mutant_counts[k] = calloc(1000, sizeof(int));
        
        for (int i=0; i<wildtype_populations[k]; ++i) {
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

void inspect_true_mutant_counts_wildtype(int*** wildtype_state, int* wildtype_populations) {
    int** true_mutant_counts = malloc(CELLS * sizeof(int*));
    get_mutant_counts_wildtype(true_mutant_counts, wildtype_state, wildtype_populations);
    print_mutant_counts(true_mutant_counts);
    for (int k=0; k<CELLS; ++k) {free(true_mutant_counts[k]);}
    free(true_mutant_counts);
    return;
}

void get_mutant_counts(int** mutant_counts, int*** wildtype_state, int*** ra_state, int* wildtype_populations, int* ra_populations) {
    int mut_id;
    int n_mut;
    int mutation_counter = 0;
    for (int k=0; k<CELLS; ++k) {
        mutant_counts[k] = calloc(1000, sizeof(int));
        
        for (int i=0; i<wildtype_populations[k]; ++i) {
            n_mut = wildtype_state[k][i][0];
            for (int j=1; j<=n_mut; ++j) {
                mut_id = wildtype_state[k][i][j];
                mutant_counts[k][mut_id]++;
                mutation_counter = fmax(mutation_counter, mut_id);
            }
        }

        for (int i=0; i<ra_populations[k]; ++i) {
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

void inspect_true_mutant_counts(int*** wildtype_state, int*** ra_state, int* wildtype_populations, int* ra_populations) {
    int** true_mutant_counts = malloc(CELLS * sizeof(int*));
    get_mutant_counts(true_mutant_counts, wildtype_state, ra_state, wildtype_populations, ra_populations);
    print_mutant_counts(true_mutant_counts);
    for (int k=0; k<CELLS; ++k) {free(true_mutant_counts[k]);}
    free(true_mutant_counts);
    return;
}
