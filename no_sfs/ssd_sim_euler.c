#include <stdio.h>
#include <stdlib.h>
// #include <unistd.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#define RND gsl_rng_uniform_pos(rng) // generate number from Unif(0,1)
#define RNORM gsl_ran_gaussian(rng, 1) // generate a normal random variate

// void drift_update(double drift[4], double state[4], double degradation_rate, double nucleus_control_factor, double target_population, double density, double diffusion_rate) {
// 	double w1 = state[0];
// 	double m1 = state[1];
// 	double w2 = state[2];
// 	double m2 = state[3];

// 	double replication_rate1 = fmax(0.0, degradation_rate + nucleus_control_factor * (target_population - w1 - density*m1));
// 	double replication_rate2 = fmax(0.0, degradation_rate + nucleus_control_factor * (target_population - w2 - density*m2));
	
// 	drift[0] = w1 * (replication_rate1 - degradation_rate - diffusion_rate) + w2 * diffusion_rate;
// 	drift[1] = m1 * (replication_rate1 - degradation_rate - diffusion_rate) + m2 * diffusion_rate;
// 	drift[2] = w2 * (replication_rate2 - degradation_rate - diffusion_rate) + w1 * diffusion_rate;
// 	drift[3] = m2 * (replication_rate2 - degradation_rate - diffusion_rate) + m1 * diffusion_rate;
// 	return;
// }

void euler_update(const gsl_rng* rng, double state[4], double time_step, double degradation_rate, double nucleus_control_factor, double target_population, double density, double diffusion_rate) {
	
	double replication_rate[4];
	replication_rate[0] = fmax(0.0, degradation_rate + nucleus_control_factor * (target_population - state[0] - density*state[1]));
	replication_rate[1] = replication_rate[0];
	replication_rate[2] = fmax(0.0, degradation_rate + nucleus_control_factor * (target_population - state[2] - density*state[3]));
	replication_rate[3] = replication_rate[2];
	
	double norms[4];
	for (int i=0; i<4; ++i) {norms[i] = RNORM;}
	for (int i=0; i<4; ++i) {
		int j = (i+2)%4;
		double drift = state[i]*(replication_rate[i]-degradation_rate-diffusion_rate)+state[j]*diffusion_rate;
		double noise = -sqrt(degradation_rate*state[i])*RNORM + sqrt(replication_rate[i]*state[i])*RNORM - sqrt(diffusion_rate*state[i])*norms[i] - sqrt(diffusion_rate*state[j])*norms[j];
		state[i] = fmax(0.0, state[i] + time_step*drift + sqrt(time_step)*noise);
	}
	return;
}

// int weighted_sample(const gsl_rng* rng, int n, double* weights) {
// 	/* Performs weighted sampling
	
// 	Inputs
// 	------
// 	rng, random number generator
//     n, number of samples
//     weights, array of weights
// 	*/
//     double weight_sum = 0.00;
//     for (int i=0; i<n; ++i){weight_sum += weights[i];}
    
//     double random_point = gsl_rng_uniform_pos(rng) * weight_sum;
//     double cum_sum=0;
    
//     for (int i=0; i<n-1; ++i) {
//         cum_sum += weights[i];
//         if (cum_sum > random_point) {
//             return i;
//         }
//     }
//     return n-1;
// }

void write_data_to_file(double state[4], double recording_time, int seed, char* filename) {
	FILE* fp = fopen(filename, "a");
	fprintf(fp, "%d,%.01f", seed, recording_time);
	for (int i=0; i<4; ++i) {fprintf(fp, ",%.2f", state[i]);}
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
	sprintf(filename, "ssd_sim_euler_%.2e_%d.txt", diffusion_rate, seed);
	FILE *fp_ssd_population = fopen(filename, "w");
	fprintf(fp_ssd_population, "sim,t,w1,m1,w2,m2\n");
	fclose(fp_ssd_population);

	for (unsigned int sim=0; sim<1000; ++sim) {
		if (sim%100==0) {printf("sim = %d\n", sim);}

		double state[4] = {25.0, 25.0, 25.0, 25.0};

		double current_time = 0.0;
		double time_step = 0.001;
		double recording_time = 0.0;
		double recording_space = 0.5;
		double final_time = 50.0;
		write_data_to_file(state, recording_time, sim, filename);
		recording_time += recording_space;

		// Gillespie algorithm until time threshold reached
		while (current_time<final_time) {
			current_time += time_step;
			euler_update(rng, state, time_step, degradation_rate, nucleus_control_factor, target_population, density, diffusion_rate);

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