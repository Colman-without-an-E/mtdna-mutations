#include <stdio.h>
#include <stdlib.h>
// #include <unistd.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#define RND gsl_rng_uniform_pos(rng) // generate number from Unif(0,1)

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

double propensity_update(double propensity[12], int state[4], double degradation_rate, double nucleus_control_factor, double target_population, double density, double diffusion_rate) {
	int w1 = state[0];
	int m1 = state[1];
	int w2 = state[2];
	int m2 = state[3];
	
	double replication_rate1 = fmax(0.0, degradation_rate + nucleus_control_factor * (target_population - w1 - density*m1));
	double replication_rate2 = fmax(0.0, degradation_rate + nucleus_control_factor * (target_population - w2 - density*m2));
	
	propensity[0] = degradation_rate * w1;
	propensity[1] = replication_rate1 * w1;
	propensity[2] = diffusion_rate * w1;
	propensity[3] = degradation_rate * m1;
	propensity[4] = replication_rate1 * m1;
	propensity[5] = diffusion_rate * m1;
	propensity[6] = degradation_rate * w2;
	propensity[7] = replication_rate2 * w2;
	propensity[8] = diffusion_rate * w2;
	propensity[9] = degradation_rate * m2;
	propensity[10] = replication_rate2 * m2;
	propensity[11] = diffusion_rate * m2;

	double propensity_sum = 0.0;
	for (int i=0; i<12; ++i) {propensity_sum += propensity[i];}
	return propensity_sum;
}

void gillespie_event(const gsl_rng* rng, double propensity[12], int state[4], double degradation_rate, double nucleus_control_factor, double target_population, double density, double diffusion_rate) {
	int event_id = weighted_sample(rng, 12, propensity);

	switch (event_id) {
		case 0: // w1 degradation
			state[0]--;
			break;
		case 1: // w1 replication
			state[0]++;
			break;
		case 2: // w1 diffusion
			state[0]--;
			state[2]++;
			break;
		case 3: // m1 degradation
			state[1]--;
			break;
		case 4: // m1 replication
			state[1]++;
			break;
		case 5: // m1 diffusion
			state[1]--;
			state[3]++;
			break;
		case 6: // w2 degradation
			state[2]--;
			break;
		case 7: // w2 replication
			state[2]++;
			break;
		case 8: // w2 diffusion
			state[2]--;
			state[0]++;
			break;
		case 9: // m2 degradation
			state[3]--;
			break;
		case 10: // m2 replication
			state[3]++;
			break;
		case 11: // m2 diffusion
			state[3]--;
			state[1]++;
			break;
	}
	return;
}

void write_data_to_file(int state[4], double recording_time, int seed, char* filename) {
	FILE* fp = fopen(filename, "a");
	fprintf(fp, "%d,%.01f", seed, recording_time);
	for (int i=0; i<4; ++i) {fprintf(fp, ",%d", state[i]);}
	fprintf(fp, "\n");
	fclose(fp);
	return;
}

int main(int argc, char *argv[]) {
    double degradation_rate = 7e-2;
    double nucleus_control_factor = 2.5e-3;
    double target_population = 30;
	double density = 0.2;

    unsigned int seed;
	double diffusion_rate;

	if (argc==3) {
		seed = atoi(argv[1]);
		diffusion_rate = atof(argv[2]);
	} else {
		printf("argc = %d\n", argc);
		printf("Usage : %s <seed> <diffusion_rate>\n", argv[0]);
		return 0;
	}

	/* set up GSL RNG */
	gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(rng, seed);
	/* end of GSL setup */

	// Set up files to write population data in
	char filename[50];
	sprintf(filename, "ssd_sim_no_sfs_%.2e_%d.txt", diffusion_rate, seed);
	FILE *fp_ssd_population = fopen(filename, "w");
	fprintf(fp_ssd_population, "sim,t,w1,m1,w2,m2\n");
	fclose(fp_ssd_population);

	for (unsigned int sim=0; sim<1000; ++sim) {
		if (sim%100==0) {printf("sim = %d\n", sim);}

		unsigned int state[4] = {25, 25, 25, 25};
		double propensity[12];
		double propensity_sum;

		double current_time = 0.0;
		double recording_time = 0.0;
		double recording_space = 0.5;
		double final_time = 50.0;
		write_data_to_file(state, recording_time, sim, filename);
		recording_time += recording_space;

		// Gillespie algorithm until time threshold reached
		while (current_time<final_time) {
			propensity_sum = propensity_update(propensity, state, degradation_rate, nucleus_control_factor, target_population, density, diffusion_rate);
			if (propensity_sum==0.0) {
				write_data_to_file(state, recording_time, sim, filename);
				recording_time += recording_space;
				break;
			}
			current_time += -log(RND) / propensity_sum;
			gillespie_event(rng, propensity, state, degradation_rate, nucleus_control_factor, target_population, density, diffusion_rate);

			// Record data at recording time
			if (current_time>=recording_time) {
				write_data_to_file(state, recording_time, sim, filename);
				recording_time += recording_space;
			}
		}
	}
    gsl_rng_free(rng);
	return 0;
}