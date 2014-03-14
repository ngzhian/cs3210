#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <malloc.h>

#define N 4

void populate(double [][N], int, int);
int rank_above_me(int, int);
int rank_below_me(int, int);

int main(int argc, char* argv[]) {
  double start, end, elapsed;
  int my_rank, p;

  MPI_Status status;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  srand(my_rank + time(NULL)); // to ensure each process gets a different seed

  int num_per_proc = N/p;
  double my_rows[num_per_proc][N];
  double my_cols[num_per_proc][N];
  double my_result[num_per_proc][N];

  populate(my_rows, num_per_proc, N);
  populate(my_cols, num_per_proc, N);

  start = MPI_Wtime();

  end = MPI_Wtime();
  elapsed = end - start;
  printf("%d took %f\n", my_rank, elapsed);

  MPI_Finalize();
}

void populate(double matrix[][N], int num, int count) {
  double r;
  int i, j;
  for (i = 0; i < num; i++) {
    for (j = 0; j < count; j++) {
      r = (double)rand()/(double)RAND_MAX;
      matrix[i][j] = r;
    }
  }
}

int rank_above_me(int p, int my_rank) {
  if (my_rank == 0) {
    return p-1;
  } else {
    return my_rank - 1;
  }
}

int rank_below_me(int p, int my_rank) {
  if (my_rank == p-1) {
    return 0;
  } else {
    return my_rank + 1;
  }
}
