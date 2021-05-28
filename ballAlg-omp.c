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
    double sum = 0;
    for (int j = 0; j < n_dims; j++) {
        double sub = point_2[j] - point_1[j];
        sum += (sub * sub);
    }

    return sqrt(sum);
}

#pragma omp declare reduction(maximum : struct pos_distance : omp_out = omp_in.distance > omp_out.distance ? omp_in : omp_out)
struct pos_distance findMostDistanceByPoint(double **points, int n_dims, long n_points, double *point_1,  int excludePos) {
    struct pos_distance posDist;
    posDist.pos = excludePos;
    posDist.distance = 0;

    if(n_points >= 1000){
        #pragma omp parallel for reduction (maximum:posDist)
        for(int i = 0; i < n_points; i++){
            if(i != excludePos){
                double distance = computeDistance(n_dims, point_1, points[i]);
                posDist.pos = distance >= posDist.distance ? i : posDist.pos ;
                posDist.distance = distance >= posDist.distance ? distance : posDist.distance;
            }
        }
    }else{
        for(int i = 0; i < n_points; i++){
            if(i != excludePos){
                double distance = computeDistance(n_dims, point_1, points[i]);
                posDist.pos = distance >= posDist.distance ? i : posDist.pos ;
                posDist.distance = distance >= posDist.distance ? distance : posDist.distance;
            }
        }
    }

    return posDist;
}

struct pairPos findMostDistance(double **points, int n_dims, long n_points, struct pairPos oldPair) {
    struct pos_distance posDistA = findMostDistanceByPoint(points, n_dims, n_points, points[oldPair.posB], oldPair.posB);
    struct pos_distance posDistB = findMostDistanceByPoint(points, n_dims, n_points, points[posDistA.pos], posDistA.pos);

    struct pairPos pairPos;
    pairPos.posA = posDistA.pos;
    pairPos.posB = posDistB.pos;

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

double *add2points(double *pointA, double *pointB, int n_dims) {
    double *result;
    result = (double *) malloc(n_dims * sizeof(double *));

    for (int j = 0; j < n_dims; j++) {
        result[j] = pointA[j] + pointB[j];
    }

    return result;
}

double *scalarMulti(double *pointA, double scalar,  int n_dims) {
    double *result;
    result = (double *) malloc(n_dims * sizeof(double *));

    for (int j = 0; j < n_dims; j++) {
        result[j] = pointA[j] * scalar;
    }

    return result;
}

struct pos_projection *findOrthogonalProjection(double **points, int n_dims, int n_points, struct pairPos AB) {
    // Po = ((alfa . beta) / (beta . beta)) * beta + a
    struct pos_projection *pointsWithProjection =
            (struct pos_projection *) malloc(n_points * sizeof(struct pos_projection));

    double *beta = subtract2points(points[AB.posB], points[AB.posA], n_dims);
    double innerProduct2 = findInnerProduct(beta, beta, n_dims);

    if(n_points > 1000){
        #pragma omp parallel for
        for (int i = 0; i < n_points; i++) {
            double *alfa = subtract2points(points[i], points[AB.posA], n_dims);
            double innerProduct1 = findInnerProduct(alfa, beta, n_dims);
            free(alfa);
            double pi = innerProduct1 / innerProduct2;
            double * scalar = scalarMulti(beta, pi, n_dims);
            double *  result = add2points(scalar, points[AB.posA], n_dims );
            free(scalar);

            pointsWithProjection[i].pos = i;
            pointsWithProjection[i].projection = result;
        }
    }else{
        for (int i = 0; i < n_points; i++) {
            double *alfa = subtract2points(points[i], points[AB.posA], n_dims);
            double innerProduct1 = findInnerProduct(alfa, beta, n_dims);
            free(alfa);
            double pi = innerProduct1 / innerProduct2;
            double * scalar = scalarMulti(beta, pi, n_dims);
            double *  result = add2points(scalar, points[AB.posA], n_dims );
            free(scalar);

            pointsWithProjection[i].pos = i;
            pointsWithProjection[i].projection = result;
        }
    }


    free(beta);
    return pointsWithProjection;
}

int cmpfunc(const void *a, const void *b) {

    struct pos_projection *orderA = (struct pos_projection *) a;
    struct pos_projection *orderB = (struct pos_projection *) b;

    if ((double) orderA->projection[0] > (double) orderB->projection[0])
        return 1;
    else if ((double) orderA->projection[0] < (double) orderB->projection[0])
        return -1;
    else
        return 0;
}

void copyPoint(double *point, double *point2, int n_dims){
    for (int i = 0; i < n_dims; i++) {
        point2[i] = point[i];
    }
}
struct node *ballTreeAlgo(double **points, int n_dims, long n_points, int id) {
    nodeLevel++;
    struct node *root;
    struct nodeValue value;
    struct pos_projection *pointsWithProjection;

    value.node_id = id;
    value.n_dims = n_dims;
    value.center_coordinates = (double *) malloc(n_dims * sizeof(double));

    if (n_points > 1) {

        struct pairPos pairPos = findPairMostDistance(points, n_dims, n_points);
        pointsWithProjection = findOrthogonalProjection(points, n_dims, n_points, pairPos);
        qsort(pointsWithProjection, n_points, sizeof(struct pos_projection), cmpfunc);

        if (n_points % 2 == 0) {
            // Compute average
            int firstCenter;
            int secondCenter;


            firstCenter = (n_points / 2);
            secondCenter = firstCenter + 1;
            int i;

            struct pos_projection posProjectionFirstCenter = pointsWithProjection[firstCenter - 1];
            struct pos_projection posProjectionSecondCenter = pointsWithProjection[secondCenter - 1];
            for (i = 0; i < n_dims; i++) {
                double firstCenterX = posProjectionFirstCenter.projection[i];
                double secondCenterX = posProjectionSecondCenter.projection[i];
                value.center_coordinates[i] = (firstCenterX + secondCenterX) / 2;
            }
        } else {
            int centerProjection = (n_points / 2) + 1;
            copyPoint( pointsWithProjection[centerProjection - 1].projection, value.center_coordinates, n_dims);
        }
    } else {
        copyPoint( points[0], value.center_coordinates, n_dims);
    }
    struct pos_distance posRadios = findMostDistanceByPoint(points, n_dims, n_points, value.center_coordinates, -1);
    value.radius = posRadios.distance;
    root = createNode(value);

    if (n_points > 1) {
        double **left = create_array_pts(n_dims, n_points);
        double **right = create_array_pts(n_dims, n_points);

        int leftPos = 0;
        int rightPos = 0;

        if(pointsWithProjection != NULL){
            for (int j = 0; j < n_points; j++) {
                if (pointsWithProjection[j].projection[0] < value.center_coordinates[0]) {
                    copyPoint( points[pointsWithProjection[j].pos], left[leftPos], n_dims);
                    leftPos++;
                } else {
                    copyPoint( points[pointsWithProjection[j].pos],  right[rightPos], n_dims);
                    rightPos++;
                }
            }
        }

        //free points with projection
        for (int i = 0; i < n_points; i++) {
            free(pointsWithProjection[i].projection);
        }
        free(pointsWithProjection);

        #pragma omp parallel
        {
            #pragma omp sections nowait
            {
                #pragma omp section
                {
                    if (leftPos > 0) {
                        global++;
                        root->left = ballTreeAlgo(left, n_dims, leftPos, global);

                        free(left[0]);
                        free(left);
                    }
                }
                #pragma omp section
                {
                    if (rightPos > 0) {
                        global++;
                        root->right = ballTreeAlgo(right, n_dims, rightPos, global);

                        free(right[0]);
                        free(right);
                    }
                }

            }
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

//    omp_set_num_threads(8);
//    int a = omp_get_max_threads();
//    printf("\nTHREADS:%d\n", a);

    //start timer
    double exec_time;
    exec_time = -omp_get_wtime();

    //get points
    points = get_points(argc, argv, &n_dims, &n_points);
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