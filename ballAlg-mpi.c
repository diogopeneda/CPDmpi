#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <inttypes.h>
#include <omp.h>
#include <mpi.h>

#include "gen_points.c"
#include "binaryTree.c"

//=================== STRUCTS =================
struct pos_projection 
{
    long point_id;
    double *projection;
    double distance_to_a;
};

struct branches
{
    struct nodeValue item;
    int left_id;
    int right_id;
};

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
    long n_nodes;
    struct branches* my_branch;
    int branch_size; 
    struct node* my_root;
    MPI_Comm my_split_comm;
    int split_process_id;
    int split_n_processes; 
};

double** points;
long n_points_total;
struct parameters* my_parameters;

//=============================================
//=================== GENERAL =================

int cmpfunc(const void *a, const void *b) 
{

    struct pos_projection *orderA = (struct pos_projection *) a;
    struct pos_projection *orderB = (struct pos_projection *) b;

    if ((double) orderA->distance_to_a > (double) orderB->distance_to_a)
        return 1;
    else if ((double) orderA->distance_to_a < (double) orderB->distance_to_a)
        return -1;
    else
        return 0;
}

//=============================================
//=================== OTHERS ==================
void printPoints(int n_dims, long n_points) 
{
    long i;
    int j;

    for (i = 0; i < n_points; i++) {
        for (j = 0; j < n_dims; j++)
            printf("[%f] ", points[i][j]);
        printf("\n");
    }
    printf("\n");
}

double relative_euclid_dist_by_point(double* point_a, double* point_b, int n_dims)
{
    double relative_dist = 0;
    double sub = 0;

    //only parallelize if over certain threshold
    //#pragma omp parallel for private(sub, relative_dist) reduction(+: relative_dist)
    for(int i = 0; i<n_dims; i++){
        sub = point_a[i] - point_b[i];
        relative_dist += sub*sub; 
    }

    return relative_dist;
}

double relative_euclid_dist_by_id(int id_a, int id_b, int n_dims)
{
    double relative_dist = 0;
    double sub = 0;

    //only parallelize if over certain threshold
    //#pragma omp parallel for private(sub, relative_dist) reduction(+: relative_dist)
    for(int i = 0; i<n_dims; i++){
        sub = points[id_a][i] - points[id_b][i];
        relative_dist += sub*sub; 
    }

    return relative_dist;
}

double *subtract2points(double *pointA, double *pointB, int n_dims) 
{
    double *result;
    result = (double *) malloc(n_dims * sizeof(double));

    for (int j = 0; j < n_dims; j++) {
        result[j] = pointA[j] - pointB[j];
    }

    return result;
}

double findInnerProduct(double *pointA, double *pointB, int n_dims) 
{
    double innerProduct = 0;

    for (int i = 0; i < n_dims; i++)
        innerProduct += pointA[i] * pointB[i];

    return innerProduct;
}

int find_furthest_others(int* indexes, int n_dims, long n_points, int start){
    int furthest;
    double new_distance;
    double longest_distance = 0;

    for(long i = 0; i<n_points; i++){
        if( i != start){
            new_distance = relative_euclid_dist_by_id(indexes[start], indexes[i], n_dims);
            if(new_distance > longest_distance){
                longest_distance = new_distance;
                furthest = i;
            }
        }
    }
    return furthest;
}

int* recursive_find_furthest_pair_others(int* indexes, int n_dims, long n_points, int* last_pair)
{
    int new_furthest_a;

    if(last_pair[1] < 0){//first iteration
        last_pair[1] = find_furthest_others(indexes, n_dims, n_points, last_pair[0]);
    }
    
    new_furthest_a = find_furthest_others(indexes, n_dims, n_points, last_pair[1]);

    if(new_furthest_a != last_pair[0]){ //not found, but know that last_pair[1] is furthest from new furthest, so need to check this connection
        last_pair[0] = last_pair[1];
        last_pair[1] = new_furthest_a;
        recursive_find_furthest_pair_others(indexes, n_dims, n_points, last_pair);
    }else{
        return last_pair;// if indexes is NULL gives points id, else gives indexes id
    }
}

int* get_furthest_pair_others(int n_dims, long n_points, int* indexes)
{
    int* furthest_pair;
    int* start_pair = (int*) malloc (2*sizeof(int));

    if(n_points > 2){
        start_pair[1] = -1;
        start_pair[0] = 0;
        furthest_pair = recursive_find_furthest_pair_others(indexes, n_dims, n_points, start_pair);
        //free(start_pair);
        return furthest_pair;
    }else if(n_points == 2){
        start_pair[1] = 1;
        start_pair[0] = 0;
            //free(furthest_pair);
        return start_pair;

    }

    printf("ERROR IN NUMBER OF POINTS!\n");
    exit(1);
}

struct pos_projection* get_ortogonal_projections_others(int* point_indexes, int n_dims, long n_points, int* furthest_pair)
{
    // Po = ((alfa . beta) / (beta . beta)) * beta + a
    struct pos_projection* point_projections = (struct pos_projection*) malloc (n_points*sizeof(struct pos_projection));
    
    double *beta;
    //printf("%ld furthest pair is: %d %d\n",n_points, furthest_pair[0], furthest_pair[1]);
    
    //printf("point ids is: %d %d\n", point_indexes[furthest_pair[0]], point_indexes[furthest_pair[1]]);
    beta = subtract2points(points[point_indexes[furthest_pair[0]]], points[point_indexes[furthest_pair[1]]], n_dims);

    double innerProduct2 = findInnerProduct(beta, beta, n_dims);
    
    for (long i = 0; i < n_points; i++) {
        
        point_projections[i].projection = (double*) malloc (n_dims*sizeof(double));
        double *alfa;
        
        alfa = subtract2points(points[point_indexes[i]], points[point_indexes[furthest_pair[0]]], n_dims);
        
        //printf("alpha is: %f\n", *alfa);
        ////printf("beta is: %f\n", *beta);
        //printf("furthest pair is: %d %d\n", furthest_pair[0], furthest_pair[1]);
        double innerProduct1 = findInnerProduct(alfa, beta, n_dims);
        
        double p = innerProduct1 / innerProduct2;
        for (int k = 0; k < n_dims; k++) {
            //build projection
            point_projections[i].projection[k] = beta[k] * p + points[point_indexes[furthest_pair[0]]][k];
        }
        //printf("point is: %f %f\n", point_projections[i].projection[0] , point_projections[i].projection[1]);
        
        //stands on the fact that p will be a value between -1 and 0, or 1 and 0, according to the
        //position in [AB]
        point_projections[i].distance_to_a = fabs(p);
        
        point_projections[i].point_id = point_indexes[i];
        
        free(alfa);
    }
    
    free(beta);
    return point_projections;
}

double get_radius_others(double* center, int* point_indexes, int n_dims, int n_points)
{
    double max_distance_squared = 0;
    double new_distance_squared = 0;

    for (int i = 0; i < n_points; i++){
        new_distance_squared = relative_euclid_dist_by_point(center, points[point_indexes[i]], n_dims);
        if(new_distance_squared > max_distance_squared){
            max_distance_squared = new_distance_squared;
        }
    }
    return sqrt(max_distance_squared);
}

struct node *ballTreeAlgo(int* point_indexes, int n_dims, long n_points, int id)
{
    struct node * root;
    struct nodeValue value;
    int n_left, n_right, *left_ids, *right_ids;
    int reversed = 0;

    value.node_id = id;
    value.n_dims = n_dims;
    value.center_coordinates = (double *) malloc(n_dims * sizeof(double));

    if(n_points >1){
        //find furthest pair
        int* furthest_pair = get_furthest_pair_others(n_dims, n_points, point_indexes);
        //find orthogonal projections
        struct pos_projection* ortogonal_projections = get_ortogonal_projections_others(point_indexes, n_dims, n_points, furthest_pair);
        //get orthogonal projections median
        /* for(int i = 0; i<n_points; i++){
            //printf("point is: %f %f\n", ortogonal_projections[i].projection[0] , ortogonal_projections[i].projection[1]);
        } */
        qsort(ortogonal_projections, n_points, sizeof(struct pos_projection), cmpfunc);
        if(points[furthest_pair[0]][0] > points[furthest_pair[1]][0]){
            reversed = 1;
        }
        //printf("YO\n");
        free(furthest_pair);
        //printf("YO\n");
        
        //median is center
        if(n_points%2 == 0){
            int firstCenter = (n_points / 2);
            int secondCenter = firstCenter - 1;
            for (int i = 0; i < n_dims; i++) {
                double firstCenterX = ortogonal_projections[firstCenter].projection[i];
                //printf("First center is: %f\n", ortogonal_projections[firstCenter].projection[i]);
                double secondCenterX = ortogonal_projections[secondCenter].projection[i];
                //printf("Second center is: %f\n", ortogonal_projections[secondCenter].projection[i]);
                value.center_coordinates[i] = (firstCenterX + secondCenterX) / 2;
            }
            //printf("center is: %f %f\n", value.center_coordinates[0], value.center_coordinates[1]);
            n_left = firstCenter;
            n_right = firstCenter;
        }else{
            int centerProjection = (n_points-1)/2;
            for (int i = 0; i < n_dims; i++) {
                value.center_coordinates[i] = ortogonal_projections[centerProjection].projection[i];
            }
            n_left = centerProjection;
            n_right = centerProjection + 1;
        }
        //find radius
        value.radius = get_radius_others(value.center_coordinates, point_indexes, n_dims, n_points);
        root = createNode(value);
        
        //create 2 index arrays for right and left points
        left_ids = (int*) malloc (n_left*sizeof(int));
        right_ids = (int*) malloc (n_right*sizeof(int));

        if(reversed == 0){
            //left to right
            for(long i=0; i<n_points; i++){
                if(i < n_left){//left
                    left_ids[i] = ortogonal_projections[i].point_id;
                }else{//right
                    right_ids[i-n_left] = ortogonal_projections[i].point_id;
                }
            }
        }else{
            //right to left
            for(long i=0; i<n_points; i++){
                if(i < n_right){//right
                    right_ids[i] = ortogonal_projections[i].point_id;
                }else{//left
                    left_ids[i-n_right] = ortogonal_projections[i].point_id;
                }
            }
        }
        
        //free orthogonal projections list
        for(long i = 0; i < n_points; i++){
            free(ortogonal_projections[i].projection);
        }
        free(ortogonal_projections);
        free(point_indexes);
        //send to recursive on left
        if(n_left > 0){
            root->left = ballTreeAlgo(left_ids, n_dims, n_left, ((id+1)*2)-1);
        }
        //send to recursive on right
        if(n_left > 0){
            root->right = ballTreeAlgo(right_ids, n_dims, n_right, ((id+1)*2));
        }
        return root;
    }

    //else jump all this, point is center itself
    for (int i = 0; i < n_dims; i++) {
        value.center_coordinates[i] = points[point_indexes[0]][i];
    }
    value.radius = 0;
    root = createNode(value);

    return root;
}

//=============================================
//==================== FIRST ==================
int find_furthest_first(int n_dims, long n_points, int start)
{
    int furthest;
    double new_distance;
    double longest_distance = 0;

    //for first iteration
    for(long i = 0; i<n_points; i++){
        if( i != start){
            new_distance = relative_euclid_dist_by_id(start, i, n_dims);
            if(new_distance > longest_distance){
                longest_distance = new_distance;
                furthest = i;
            }
        }
    }
    return furthest;
}

int* recursive_find_furthest_pair_first(int n_dims, long n_points, int* last_pair)
{
    int new_furthest_a;

    if(last_pair[1] < 0){//first iteration
        last_pair[1] = find_furthest_first(n_dims, n_points, last_pair[0]);
    }
    
    new_furthest_a = find_furthest_first(n_dims, n_points, last_pair[1]);

    if(new_furthest_a != last_pair[0]){ //not found, but know that last_pair[1] is furthest from new furthest, so need to check this connection
        last_pair[0] = last_pair[1];
        last_pair[1] = new_furthest_a;
        recursive_find_furthest_pair_first(n_dims, n_points, last_pair);
    }else{
        return last_pair;// if indexes is NULL gives points id, else gives indexes id
    }
}

int* get_furthest_pair_first(int n_dims, long n_points)
{
    int* furthest_pair;
    int* start_pair = (int*) malloc (2*sizeof(int));

    if(n_points > 2){
        start_pair[1] = -1;
        start_pair[0] = 0;
        furthest_pair = recursive_find_furthest_pair_first(n_dims, n_points, start_pair);
        //free(start_pair);
        return furthest_pair;
    }else if(n_points == 2){
        start_pair[1] = 1;
        start_pair[0] = 0;
            //free(furthest_pair);
        return start_pair;

    }

    printf("ERROR IN NUMBER OF POINTS!\n");
    exit(1);
}

struct pos_projection* get_ortogonal_projections_first(int n_dims, long n_points, int* furthest_pair)
{
    // Po = ((alfa . beta) / (beta . beta)) * beta + a
    struct pos_projection* point_projections = (struct pos_projection*) malloc (n_points*sizeof(struct pos_projection));
    
    double *beta;
    //printf("%ld furthest pair is: %d %d\n",n_points, furthest_pair[0], furthest_pair[1]);
    
    beta = subtract2points(points[furthest_pair[0]], points[furthest_pair[1]], n_dims);

    double innerProduct2 = findInnerProduct(beta, beta, n_dims);
    
    for (long i = 0; i < n_points; i++) {   
        point_projections[i].projection = (double*) malloc (n_dims*sizeof(double));
        double *alfa;
        
        alfa = subtract2points(points[i], points[furthest_pair[0]], n_dims);
        //printf("alpha is: %f\n", *alfa);
        ////printf("beta is: %f\n", *beta);
        //printf("furthest pair is: %d %d\n", furthest_pair[0], furthest_pair[1]);
        double innerProduct1 = findInnerProduct(alfa, beta, n_dims);
        
        double p = innerProduct1 / innerProduct2;
        for (int k = 0; k < n_dims; k++) {
            //build projection
            point_projections[i].projection[k] = beta[k] * p + points[furthest_pair[0]][k];
        }
        //printf("point is: %f %f\n", point_projections[i].projection[0] , point_projections[i].projection[1]);
        
        //stands on the fact that p will be a value between -1 and 0, or 1 and 0, according to the
        //position in [AB]
        point_projections[i].distance_to_a = fabs(p);
        
        point_projections[i].point_id = i;
        
        free(alfa);
    }
    
    free(beta);
    return point_projections;
}

double get_radius_first(double* center, int n_dims, int n_points)
{
    double max_distance_squared = 0;
    double new_distance_squared = 0;

    for (int i = 0; i < n_points; i++){
        new_distance_squared = relative_euclid_dist_by_point(center, points[i], n_dims);
        if(new_distance_squared > max_distance_squared){
            max_distance_squared = new_distance_squared;
        }
    }
    return sqrt(max_distance_squared);

}

void start_algo()
{
    my_parameters->n_iteration++;
    int n_left, n_right, *left_ids, *right_ids;
    int reversed = 0;

    
        struct nodeValue value;
        value.node_id = 0;
        value.n_dims = my_parameters->n_dims;
        value.center_coordinates = (double *) malloc(my_parameters->n_dims * sizeof(double));
        
        if(n_points_total>1){
            int* furthest_pair = get_furthest_pair_first(my_parameters->n_dims, n_points_total);
            struct pos_projection* ortogonal_projections = get_ortogonal_projections_first(my_parameters->n_dims, n_points_total, furthest_pair);
            qsort(ortogonal_projections, n_points_total, sizeof(struct pos_projection), cmpfunc);
            if(points[furthest_pair[0]][0] > points[furthest_pair[1]][0]){
                reversed = 1;
            }
            free(furthest_pair);
            if(n_points%2 == 0){
                int firstCenter = (n_points_total / 2);
                int secondCenter = firstCenter - 1;
                for (int i = 0; i < n_dims; i++) {
                    double firstCenterX = ortogonal_projections[firstCenter].projection[i];
                    double secondCenterX = ortogonal_projections[secondCenter].projection[i];
                    value.center_coordinates[i] = (firstCenterX + secondCenterX) / 2;
                }
                n_left = firstCenter;
                n_right = firstCenter;
            }else{
                int centerProjection = (n_points_total-1)/2;
                for (int i = 0; i < n_dims; i++) {
                    value.center_coordinates[i] = ortogonal_projections[centerProjection].projection[i];
                }
                n_left = centerProjection;
                n_right = centerProjection + 1;
            }
            value.radius = get_radius_first(value.center_coordinates, my_parameters->n_dims, n_points_total);
            if(my_parameters->n_processes == 1){//serial implementation
                my_parameters->my_root = createNode(value);
            }else{
                my_parameters->my_branch[0].item = value;
                my_parameters->my_branch[0].left_id = 1;
                my_parameters->my_branch[0].right_id = 2;
            }   
            left_ids = (int*) malloc (n_left*sizeof(int));
            right_ids = (int*) malloc (n_right*sizeof(int));
            if(reversed == 0){
                //left to right
                for(long i=0; i<n_points; i++){
                    if(i < n_left){//left
                        left_ids[i] = ortogonal_projections[i].point_id;
                    }else{//right
                        right_ids[i-n_left] = ortogonal_projections[i].point_id;
                    }
                }
            }else{
                //right to left
                for(long i=0; i<n_points; i++){
                    if(i < n_right){//right
                        right_ids[i] = ortogonal_projections[i].point_id;
                    }else{//left
                        left_ids[i-n_right] = ortogonal_projections[i].point_id;
                    }
                }
            }
            for(long i = 0; i < n_points; i++){
                free(ortogonal_projections[i].projection);
            }
            free(ortogonal_projections);
            if(my_parameters->n_processes == 1){
                if(n_left > 0){        
                    my_parameters->my_root->left = ballTreeAlgo(left_ids, my_parameters->n_dims, n_left, 1); 
                } 
                //send to recursive on right
                if(n_right > 0){
                    my_parameters->my_root->right = ballTreeAlgo(right_ids, my_parameters->n_dims, n_right, 2);
                }
            }else{
                //split comm
                if(my_parameters->process_id%2 == 0){
                    MPI_Comm_split(MPI_COMM_WORLD, 1, my_parameters->process_id, &my_parameters->my_split_comm);
                    MPI_Comm_size(MPI_COMM_WORLD, &my_parameters->split_n_processes);
                    MPI_Comm_rank(MPI_COMM_WORLD, &my_parameters->split_process_id);
                    if(n_left > 0){
                        if(my_parameters->split_n_processes == 1){
                            my_parameters->my_root = ballTreeAlgo(left_ids, my_parameters->n_dims, n_left, 1);
                        }else{
                            //ballTreeAlgoBranch(left_ids, my_parameters->n_dims, n_left, 1);
                        }        
                         
                    } 
                }else{
                    MPI_Comm_split(MPI_COMM_WORLD, 2, my_parameters->process_id, &my_parameters->my_split_comm);
                    MPI_Comm_size(MPI_COMM_WORLD, &my_parameters->split_n_processes);
                    MPI_Comm_rank(MPI_COMM_WORLD, &my_parameters->split_process_id);
                    if(n_right > 0){
                        if(my_parameters->split_n_processes == 1){
                            my_parameters->my_root = ballTreeAlgo(right_ids, my_parameters->n_dims, n_right, 2);
                        }else{
                            //ballTreeAlgoBranch(right_ids, my_parameters->n_dims, n_right, 2);
                        } 
                        
                    }
                }
            }
            return;
        }
        //else jump all this, point is center itself
        for (int i = 0; i < my_parameters->n_dims; i++) {
            value.center_coordinates[i] = points[0][i];
        }
        value.radius = 0;
        my_parameters->my_root = createNode(value);
        

        return;
}


int main(int argc, char *argv[]) {
    // TODO error handling initial parameters
    //int n_dims;
    //long n_points;
    //int seed;
    struct node *root;
    //double **points;

    //start timer
    double exec_time;
    //exec_time = -omp_get_wtime();
    

    MPI_Init(&argc, &argv);
    double t1 = MPI_Wtime();

    my_parameters = (struct parameters*) malloc (sizeof(struct parameters));

    //get number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &my_parameters->n_processes);
    //get process rank
    MPI_Comm_rank(MPI_COMM_WORLD, &my_parameters->process_id);

    //get points
    points = get_points(argc, argv, &my_parameters->n_dims, &n_points_total);
    //calculate number of points
    long n_nodes_aux = 0;
    long max_pow = ceil(log2(n_points));
    for(int i = 0; i<ceil(log2(n_points_total)); i++){
        n_nodes_aux = n_nodes_aux + pow(2, i);
    }
    my_parameters->n_nodes = n_nodes_aux + (n_points_total - pow(2,max_pow-1))*2;
    my_parameters->my_elements = NULL;
    my_parameters->n_my_points = n_points_total;
    my_parameters->n_iteration = 0;
    my_parameters->my_team = (int*) malloc(my_parameters->n_processes*sizeof(int));
    for(int i = 0; i<my_parameters->n_processes; i++){
        my_parameters->my_team[i] = i;
    }
    my_parameters->branch_size = (int)ceil(log2(my_parameters->n_processes));
    my_parameters->my_branch = (struct branches*) malloc (my_parameters->branch_size*sizeof(struct branches));


    //do first call
    //check every time if left or right
    //change how nodes are created
    start_algo();

    double t2 = MPI_Wtime();
    exec_time = t2 - t1;

    

    for(int i = 0; i < my_parameters->n_processes; i++){
        if(i == my_parameters->process_id){
            if(i == 0){
                //Print first line
                printf("%d %ld \n", my_parameters->n_dims, my_parameters->n_nodes);
            }
            for(int j = 0; j<my_parameters->branch_size; j++){
                if(j == my_parameters->process_id){
                    printf("%d %d %d %f", my_parameters->my_branch[j].item.node_id, my_parameters->my_branch[j].left_id, my_parameters->my_branch[j].right_id, my_parameters->my_branch[j].item.radius);
                    for(int k = 0; k < my_parameters->n_dims; k++){
                        printf(" %f", my_parameters->my_branch[j].item.center_coordinates[j]);
                    }
                }
            }
            postorderTraversal(my_parameters->my_root);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    

    

    //stop timer
    if(my_parameters->process_id == 0){
        fprintf(stderr, "%.1lf\n", exec_time);
    }
    

    printf("\n");
    clean(root);
    free(points[0]);
    free(points);

    //should happen only on rank 0
    MPI_Finalize();
    return 0;
}