#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <inttypes.h>
#include <omp.h>
#include <mpi.h>

#include "gen_points.c"
#include "binaryTree.c"

struct parameters
{
    int n_dims;
    long n_my_points;
    int process_id;
    int n_processes;
    int* my_elements;
    int* my_team; //includes myself
    int n_iteration; //what iteration the whole program is expected to be at
    int max_parallel_iteration; 
};

double** points;
long n_points_total;

