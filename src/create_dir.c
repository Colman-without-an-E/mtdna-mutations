#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>

#define DIR_LOC "../data"

int main() {

    // Create directory, named with current time, to store simulation results
	// time_t t = time(NULL);
    // struct tm *tm_info = localtime(&t);
    // char dir_name[20];
    // strftime(dir_name, 20, "%Y-%m-%d_%H-%M-%S", tm_info);
    // char dir_loc[256];
    // sprintf(dir_loc, "%s/%s", DIR_LOC, dir_name);
    // puts(dir_loc);
    // if(mkdir(dir_loc)) {
    //     perror("Error");
    //     return 1;
    // }

    // // Change to directory which stores simulation results
    // if (chdir(dir_loc)) {
    //     perror("Error");
    //     return 1;
    // }

    // FILE *fp = fopen("textFile.txt" ,"w");
    // fprintf(fp, "hello\n");
    // fclose(fp);

    // char ra_population_filename[30] = "ra_sim_populations.txt";
	// char ra_sfs_filename[40] = "ra_sim_site_frequency_spectrum.txt";
    // puts(ra_population_filename);
    // puts(ra_sfs_filename);
	printf("n_sims,%f\n", 0.00000000253);
    return 0;
}
