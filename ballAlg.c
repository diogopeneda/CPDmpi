#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <inttypes.h>
#include <omp.h>

#include "gen_points.c"
#include "binaryTree.c"

_Atomic int global = 0;
_Atomic int nodeLevel = 0;

struct pos_distance {
    int pos;
    double distance;
};

struct pairPos {
    int posA;
    int posB;
};

struct pos_projection {
    int pos;
    double *projection;
    double distance_to_a;
};

int printPosProjections(struct pos_projection *projections, int n_dims, long n_points) {
    long i;
    int j;

    for (i = 0; i < n_points; i++) {
        printf("{%d} ", projections[i].pos);
        for (j = 0; j < n_dims; j++) {
            printf("[%f] ", projections[i].projection[j]);
        }

        printf("\n");
    }
    printf("\n");

    return 0;
}

int printPoints(double **points, int n_dims, long n_points) {
    long i;
    int j;

    for (i = 0; i < n_points; i++) {
        for (j = 0; j < n_dims; j++)
            printf("[%f] ", points[i][j]);
        printf("\n");
    }
    printf("\n");

    return 0;
}

//int print_point(double *point, int n_dims) {
//    int j;
//
//    for (j = 0; j < n_dims; j++)
//        printf("[%f] ", point[j]);
//
//    printf("\n");
//    printf("\n");
//
//    return 0;
//}

double computeDistance(int n_dims, double *point_1, double *point_2) {
    int j;
    double sum = 0;

    for (j = 0; j < n_dims; j++) {
        double sub = point_2[j] - point_1[j];
        sum += (sub * sub);
    }

    return sqrt(sum);
}

struct pos_distance *findMostDistanceByPoint(double **points, int n_dims, long n_points, double *point_1,  int excludePos) {
    int i;
    struct pos_distance *posDist = (struct pos_distance *) malloc(sizeof(struct pos_distance));
    posDist->pos = 0;
    posDist->distance = computeDistance(n_dims, point_1, points[0]);

    for(i = 1; i < n_points; i++){
        double distance = 0;
        if(i != excludePos){
            distance = computeDistance(n_dims, point_1, points[i]);
        }

        if(distance >= posDist->distance){
            posDist->pos = i;
            posDist->distance = distance;
        }
    }

    return posDist;
}

struct pairPos findMostDistance(double **points, int n_dims, long n_points, struct pairPos oldPair) {
    struct pos_distance *posDistA = findMostDistanceByPoint(points, n_dims, n_points, points[oldPair.posB], oldPair.posB);
    struct pos_distance *posDistB = findMostDistanceByPoint(points, n_dims, n_points, points[posDistA->pos], posDistA->pos);

    struct pairPos pairPos;
    pairPos.posA = posDistA->pos;
    pairPos.posB = posDistB->pos;

    free(posDistA);
    free(posDistB);

    if (oldPair.posA == pairPos.posA && oldPair.posB == pairPos.posB) {
        return pairPos;
    }

    return findMostDistance(points, n_dims, n_points, pairPos);
}

struct pairPos findPairMostDistance(double **points, int n_dims, long n_points) {
    struct pairPos pairPos;
    pairPos.posB = 0;

    if (n_points == 2) {
        pairPos.posA = 0;
        pairPos.posB = 1;
    } else {
        pairPos = findMostDistance(points, n_dims, n_points, pairPos);

    }
    return pairPos;
}

double findInnerProduct(double *pointA, double *pointB, int n_dims) {
    double innerProduct = 0;

    for (int i = 0; i < n_dims; i++)
        innerProduct += pointA[i] * pointB[i];

    return innerProduct;
}

double *subtract2points(double *pointA, double *pointB, int n_dims) {
    double *result;
    result = (double *) malloc(n_dims * sizeof(double *));

    for (int j = 0; j < n_dims; j++) {
        result[j] = pointA[j] - pointB[j];
    }

    return result;
}


struct pos_projection *findOrthogonalProjection(double **points, int n_dims, int n_points, struct pairPos AB) {
    // Po = ((alfa . beta) / (beta . beta)) * beta + a
    double *alfa;
    struct node *root;

    struct pos_projection *pointsWithProjection = (struct pos_projection *) malloc(
            n_points * sizeof(struct pos_projection));

    double *beta = subtract2points(points[AB.posB], points[AB.posA], n_dims);
    double innerProduct2 = findInnerProduct(beta, beta, n_dims);

    // run each point
    for (int i = 0; i < n_points; i++) {

        pointsWithProjection[i].pos = i;
        pointsWithProjection[i].projection = (double *) malloc(n_dims * sizeof(double));

        alfa = subtract2points(points[i], points[AB.posA], n_dims);
        double innerProduct1 = findInnerProduct(alfa, beta, n_dims);
        double pi = innerProduct1 / innerProduct2;

        for (int k = 0; k < n_dims; k++) {
            //build projection
            pointsWithProjection[i].projection[k] = beta[k] * pi + points[AB.posA][k];
        }
        pointsWithProjection[i].distance_to_a = fabs(pi);
        free(alfa);
    }
    free(beta);
    return pointsWithProjection;
}

int cmpfunc(const void *a, const void *b) {

    struct pos_projection *orderA = (struct pos_projection *) a;
    struct pos_projection *orderB = (struct pos_projection *) b;

    if ((double) orderA->distance_to_a > (double) orderB->distance_to_a)
        return 1;
    else if ((double) orderA->distance_to_a < (double) orderB->distance_to_a)
        return -1;
    else
        return 0;
}

struct node *ballTreeAlgo(double **points, int n_dims, long n_points, int id) {
    nodeLevel++;

    struct node *root;
    struct nodeValue value;
    struct pos_projection *pointsWithProjection;

    long n_left, n_right;
    int reversed = 0;

    value.node_id = id;
    value.n_dims = n_dims;
    value.center_coordinates = (double *) malloc(n_dims * sizeof(double));

    if (n_points > 1) {
        struct pairPos pairPos;
        pairPos = findPairMostDistance(points, n_dims, n_points);
        pointsWithProjection = findOrthogonalProjection(points, n_dims, n_points, pairPos);

        qsort(pointsWithProjection, n_points, sizeof(struct pos_projection), cmpfunc);
        if(pointsWithProjection[pairPos.posA][0] > pointsWithProjection[pairPos.posB][0]){
            reversed = 1;
        }

        if (n_points % 2 == 0) {
            // Compute average
            int firstCenter;
            int secondCenter;


            firstCenter = (n_points / 2);
            secondCenter = firstCenter - 1;
            int i;

            struct pos_projection posProjectionFirstCenter = pointsWithProjection[firstCenter];
            struct pos_projection posProjectionSecondCenter = pointsWithProjection[secondCenter];
            for (i = 0; i < n_dims; i++) {
                double firstCenterX = posProjectionFirstCenter.projection[i];
                double secondCenterX = posProjectionSecondCenter.projection[i];
                value.center_coordinates[i] = (firstCenterX + secondCenterX) / 2;
            }
            n_left = firstCenter;
            n_right = firstCenter;
        } else {
            int centerProjection = (n_points-1)/1;
            for (int i = 0; i < n_dims; i++) {
                value.center_coordinates[i] = pointsWithProjection[centerProjection].projection[i];
            }
            n_left = centerProjection;
            n_right = centerProjection + 1;
        }
    } else {
        for (int i = 0; i < n_dims; i++) {
            value.center_coordinates[i] = points[0][i];
        }
    }
    struct pos_distance *posRadios = findMostDistanceByPoint(points, n_dims, n_points, value.center_coordinates, -1);
    value.radius = posRadios->distance;
    free(posRadios);
    root = createNode(value);

    int j;
    //int leftPos = 0;
    //int rightPos = 0;

    if (n_points > 1) {
        double **left = create_array_pts(n_left, n_points);
        double **right = create_array_pts(n_right, n_points);

        //leftPos = 0;
        //rightPos = 0;
        /* for (j = 0; j < n_points && pointsWithProjection != NULL; j++) {
            if (pointsWithProjection[j].projection[0] < value.center_coordinates[0]) {
                int z;
                for (z = 0; z < n_dims; z++) {
                    left[leftPos][z] = points[pointsWithProjection[j].pos][z];
                }
                leftPos++;
            } else {
                int z;
                for (z = 0; z < n_dims; z++) {
                    right[rightPos][z] = points[pointsWithProjection[j].pos][z];
                }
                rightPos++;
            }
        } */
        if(reversed == 0){
            for(long i = 0; i<n_points, i++){
                if(i < n_left){
                    for (int z = 0; z < n_dims; z++) {
                        left[i][z] = points[pointsWithProjection[i].pos][z];
                    }
                }else{
                    for (int z = 0; z < n_dims; z++) {
                        right[i-n_left][z] = points[pointsWithProjection[i].pos][z];
                    }
                }
            }
        }else{
            for(long i = 0; i<n_points, i++){
                if(i < n_right){
                    for (int z = 0; z < n_dims; z++) {
                        right[i][z] = points[pointsWithProjection[i].pos][z];
                    }
                }else{
                    for (int z = 0; z < n_dims; z++) {
                        left[i-n_right][z] = points[pointsWithProjection[i].pos][z];
                    }
                }
            }
        }

        //free points with projection
        for (int i = 0; i < n_points; i++) {
            free(pointsWithProjection[i].projection);
        }
        free(pointsWithProjection);

        if (n_left > 0) {
            global++;
            root->left = ballTreeAlgo(left, n_dims, leftPos, global);

            free(left[0]);
            free(left);
        }

        if (n_right > 0) {
            global++;
            root->right = ballTreeAlgo(right, n_dims, rightPos, global);

            free(right[0]);
            free(right);
        }
    }
    return root;
}

int main(int argc, char *argv[]) {
    // TODO error handling initial parameters
    int n_dims;
    long n_points;
    int seed;
    struct node *root;
    double **points;

    //start timer
    double exec_time;
    exec_time = -omp_get_wtime();

    //get points
    points = get_points(argc, argv, &n_dims, &n_points);

    //printPoints(points, n_dims, n_points);
    root = ballTreeAlgo(points, n_dims, n_points, 0);
    exec_time += omp_get_wtime();

    //Print first line
    printf("%d %d \n", n_dims, nodeLevel);
    postorderTraversal(root);

    //stop timer
    fprintf(stderr, "%.1lf\n", exec_time);

    printf("\n");
    clean(root);
    free(points[0]);
    free(points);
    return -1;
}