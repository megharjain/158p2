#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

unsigned int matrix_checksum(int N, double *M);

void check_args(int argc, char **argv, int *resolution, double *fire_temp, double *wall_temp, double *epsilon) 
{
    if (argc != 5) {
        fprintf(stderr, "%s", "Usage: ./heat_distribution_omp N fire_temp wall_temp epsilon\n");
        exit(1);
    }

    *resolution = atoi(argv[1]);
    *fire_temp = atof(argv[2]);
    *wall_temp = atof(argv[3]);
    *epsilon = atof(argv[4]);

    if (*resolution < 3 || *resolution > 2000) {
        fprintf(stderr, "%s", "Error: wrong map order (3 <= N <= 2000)\n");
        exit(1);
    }
    if (*fire_temp < 0 || *fire_temp > 100) {
        fprintf(stderr, "%s", "Error: wrong fire temperature (0.000000 <= N <= 100.000000)\n");
        exit(1);
    }
    if (*wall_temp < 0 || *wall_temp > 100) {
        fprintf(stderr, "%s", "Error: wrong wall temperature (0.000000 <= N <= 100.000000)\n");
        exit(1);
    }
    if (*epsilon < 1e-6 || *epsilon > 100) {
        fprintf(stderr, "%s", "Error: wrong epsilon (0.000001 <= N <= 100.000000)\n");
        exit(1);
    }
    return;
}

void initialize_edges(double *room, int resolution, double fire_temp, double wall_temp)
{
    #pragma omp parallel for  
    for (int i = 0; i < resolution; i++) {
        room[i*resolution] = wall_temp;
        room[i*resolution+resolution-1] = wall_temp;
        room[i] = fire_temp;
        room[(resolution-1)*resolution+i] = wall_temp;
    }
    return;
}

double initialize_interior(double *room, int resolution)
{
   
    double mean_value = 0;

    #pragma omp parallel for reduction(+:mean_value) schedule (static,4)
    for (int i = 0; i < resolution; i++) {
        mean_value = mean_value + room[i] + room[(resolution-1)*resolution+i];
        
        if (i > 0 && i < resolution-1)
        {
            mean_value = mean_value + room[i*resolution] + room[i*resolution+(resolution-1)];
        }
    } 
        
    
    mean_value = mean_value/(4*resolution-4);
    
    #pragma omp parallel for schedule (static)
    for (int i = 1; i < resolution - 1; i++) {
        
        for (int j = 1; j < resolution -1; j++) {
            room[i*resolution+j] = mean_value;
        }
    }
    return mean_value;
}

void start_time(struct timespec *start)
{
    clock_gettime(CLOCK_MONOTONIC, start); 
}

void print_time(struct timespec *start, struct timespec *end)
{
    clock_gettime(CLOCK_MONOTONIC, end);
    double elapsed = (end->tv_sec+end->tv_nsec/1e9)-(start->tv_sec+start->tv_nsec/1e9);
    printf("Running time: %f secs\n", elapsed);
}

void compute_dist(double *room, int resolution, double epsilon)
{

    double max = epsilon+0.01;
    double *tmp = malloc(2*(resolution*sizeof(double))*(resolution*sizeof(double)));
    int iteration = 0;
    int curr = 0;
    int next = 1;
    
       
    memcpy(tmp, room, (resolution*sizeof(double))*(resolution*sizeof(double)));
    memcpy(&tmp[resolution*resolution], room, (resolution*sizeof(double))*(resolution*sizeof(double)));
    
    while (max >= epsilon) {
        max = 0;
        iteration++;
    
        
        #pragma omp parallel for reduction (max:max) schedule (dynamic)
        for (int i = 1; i < resolution - 1; i++) {
            
            for (int j = 1; j < resolution - 1; j++) {
                double old = tmp[(curr*resolution*resolution)+(i*resolution)+j];
                double new = (tmp[(curr*resolution*resolution)+((i-1)*resolution)+j]+tmp[(curr*resolution*resolution)+((i+1)*resolution)+j]+tmp[(curr*resolution*resolution)+(i*resolution)+j-1]+tmp[(curr*resolution*resolution)+(i*resolution)+j+1])/4;
                tmp[(next*resolution*resolution)+(i*resolution)+j] = new;
                float diff = fabs(new - old);
                if (diff > max) max = diff;
            }
        }
        
        curr = next;
        next = (curr+1) % 2;
        if ((iteration & (iteration - 1)) == 0)
            printf("%-7d %f\n", iteration, max);
    }
    
    memcpy(room, &tmp[curr*resolution*resolution], (resolution*sizeof(double))*(resolution*sizeof(double)));
    free(tmp);
    printf("%-7d %f\n", iteration, max);
    return;
}

int main(int argc, char **argv)
{
    int resolution;
    double fire_temp, wall_temp, epsilon;
    double *room;
    double mean_value;
    struct timespec *start = malloc(sizeof(struct timespec));
    struct timespec *end = malloc(sizeof(struct timespec));

    check_args(argc, argv, &resolution, &fire_temp, &wall_temp, &epsilon);

    room = malloc((resolution*sizeof(double))*(resolution*sizeof(double)));

    start_time(start);
    initialize_edges(room, resolution, fire_temp, wall_temp);
    mean_value = initialize_interior(room, resolution);
    print_time(start, end);

    printf("mean: %f\n", mean_value);
    printf("hmap: %u\n", matrix_checksum(resolution, room));

    start_time(start);
    compute_dist(room, resolution, epsilon);
    print_time(start, end);

    printf("hmap: %u\n", matrix_checksum(resolution, room));

    free(start);
    free(end);
    free(room);
    return 0;
}
