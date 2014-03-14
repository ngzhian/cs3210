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

  int round, r, c, i;
  for (round = 0; round < p; round++) {
    double result = 0;
    // compute locally
    for (r = 0; r < num_per_proc; r++) {
      for (c = 0; c < num_per_proc; c++) {
        for (i = 0; i < N; i++) {
          result += my_rows[r][i] * my_cols[c][i];
        }
      my_result[c][round*num_per_proc + r] = result;
      }
    }
    if (round+1 == p) break; // break early to avoid unnecessary communication

    if (my_rank % 2 == 0) {
      int send_to = rank_above_me(p, my_rank);
      MPI_Send(my_rows, num_per_proc*N, MPI_DOUBLE, send_to, my_rank, MPI_COMM_WORLD);
      int receive_from = rank_below_me(p, my_rank);
      MPI_Recv(my_rows, num_per_proc*N, MPI_DOUBLE, receive_from, receive_from, MPI_COMM_WORLD, &status);
    } else {
      int receive_from = rank_below_me(p, my_rank);
      MPI_Status status;
      double buffer[num_per_proc][N]; // need to store in buffer else it will overwrite what we send
      MPI_Recv(buffer, N*num_per_proc, MPI_DOUBLE, receive_from, receive_from, MPI_COMM_WORLD, &status);
      int send_to = rank_above_me(p, my_rank);
      MPI_Send(my_rows, N*num_per_proc, MPI_DOUBLE, send_to, my_rank, MPI_COMM_WORLD);
      // copy numbers from buffer to row
      int i, j;
      for (i = 0; i < num_per_proc; i++) {
        for (j = 0; j < N; j++) {
          my_rows[i][j] = buffer[i][j];
        }
      }
    }
  }

  // gather results
  double results[N][N];
  MPI_Gather(my_result,N*num_per_proc, MPI_DOUBLE, results, N*num_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);

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
