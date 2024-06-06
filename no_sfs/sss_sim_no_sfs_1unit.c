#include <stdio.h>
#include <stdlib.h>
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

double propensity_update(double propensity[4], unsigned int state[2], double degradation_rate, double nucleus_control_factor, double target_population, double rate_diff) {
	unsigned int w = state[0];
	unsigned int m = state[1];
	
	double replication_rate = fmax(0.0, degradation_rate + nucleus_control_factor * (target_population - w - m));
	
	propensity[0] = (degradation_rate + rate_diff) * w;
	propensity[1] = (replication_rate + rate_diff) * w;
	propensity[2] = degradation_rate * m;
	propensity[3] = replication_rate * m;

	double propensity_sum = 0.0;
	for (int i=0; i<4; ++i) {propensity_sum += propensity[i];}
	return propensity_sum;
}

void gillespie_event(const gsl_rng* rng, double propensity[4], unsigned int state[2], double degradation_rate, double nucleus_control_factor, double target_population, double density) {
	int event_id = weighted_sample(rng, 4, propensity);

	switch (event_id) {
		case 0: // w degradation
			state[0]--;
			break;
		case 1: // w replication
			state[0]++;
			break;
		case 2: // m degradation
			state[1]--;
			break;
		case 3: // m replication
			state[1]++;
			break;
	}
	return;
}

void write_data_to_file(unsigned int state[2], double recording_time, int seed, char* filename) {
	FILE* fp = fopen(filename, "a");
	fprintf(fp, "%d,%.01f", seed, recording_time);
	for (int i=0; i<2; ++i) {fprintf(fp, ",%d", state[i]);}
	fprintf(fp, "\n");
	fclose(fp);
	return;
}

int main(int argc, char *argv[]) {
    double degradation_rate = 7e-2;
    double nucleus_control_factor = 2.5e-3;
    double target_population = 60;
	double rate_diff = 2.5e-2;

    unsigned int seed;

	if (argc==2) {
		seed = atoi(argv[1]);
	} else {
		printf("argc = %d\n", argc);
		printf("Usage : %s <seed>\n", argv[0]);
		return 0;
	}

	/* set up GSL RNG */
	gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(rng, seed);
	/* end of GSL setup */

	// Set up files to write population data in
	char filename[50];
	sprintf(filename, "sss_sim_no_sfs_1unit_%d.txt", seed);
	FILE *fp_ssd_population = fopen(filename, "w");
	fprintf(fp_ssd_population, "sim,t,w,m\n");
	fclose(fp_ssd_population);

	for (unsigned int sim=0; sim<5000; ++sim) {
		if (sim%100==0) {printf("sim = %d\n", sim);}

		unsigned int state[2] = {50, 50};
		double propensity[4];
		double propensity_sum;

		double current_time = 0.0;
		double recording_time = 0.0;
		double recording_space = 0.5;
		double final_time = 50.0;
		write_data_to_file(state, recording_time, sim, filename);
		recording_time += recording_space;

		// Gillespie algorithm until time threshold reached
		while (current_time<final_time) {
			propensity_sum = propensity_update(propensity, state, degradation_rate, nucleus_control_factor, target_population, rate_diff);
			if (propensity_sum==0.0) {
				write_data_to_file(state, recording_time, sim, filename);
				recording_time += recording_space;
				break;
			}
			current_time += -log(RND) / propensity_sum;
			gillespie_event(rng, propensity, state, degradation_rate, nucleus_control_factor, target_population, rate_diff);

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