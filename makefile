all: program1 program2
program1: ballAlg.c
		gcc -fopenmp -O3 -o ballAlg ballAlg.c -lm

program2: ballAlg-omp.c
		gcc -fopenmp -O3 -o ballAlg-omp ballAlg-omp.c -lm

program3: ballAlg-mpi.c
		mpicc -fopenmp -g -o ballAlg-mpi ballAlg-mpi.c -lm