#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <malloc.h>

#define N 416

void populate(double [N][N]);
void transpose(double matrix[N][N], double transposed[N][N]);
int rank_above_me(int, int);
int rank_below_me(int, int);

int main(int argc, char* argv[]) {
  double start, end, elapsed;
  int my_rank;
  int p;

  MPI_Status status;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  int num_per_proc = N/p;

  double matrixA[N][N];
  double matrixB[N][N];
  double matrixT[N][N];

  if (my_rank == 0) {
    populate(matrixA);
    populate(matrixB);
    transpose(matrixB, matrixT);
  }

  start = MPI_Wtime();

  double row[num_per_proc][N];
  double col[num_per_proc][N];

  MPI_Scatter(matrixA, num_per_proc*N, MPI_DOUBLE, row, num_per_proc*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatter(matrixT, num_per_proc*N, MPI_DOUBLE, col, num_per_proc*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  int i;

  //compute
  double result = 0;
  int round = 0;
  double my_c_row[N][num_per_proc];
  int r, c;
  for (r = 0; r < num_per_proc; r++) {
    for (c = 0; c < num_per_proc; c++) {
      for (i = 0; i < N; i++) {
        result += row[r][i]*col[c][i];
      }
      my_c_row[r][c] = result; 
    }
    round++;
  }
  // send up the row
  while (round < N) {
    if (my_rank % 2 == 0) {
      int send_to = rank_above_me(p, my_rank);
      MPI_Send(row, num_per_proc*N, MPI_DOUBLE, send_to, my_rank, MPI_COMM_WORLD);
      int receive_from = rank_below_me(p, my_rank);
      MPI_Recv(row, num_per_proc*N, MPI_DOUBLE, receive_from, receive_from, MPI_COMM_WORLD, &status);
    } else {
      int receive_from = rank_below_me(p, my_rank);
      MPI_Status status;
      MPI_Recv(row, num_per_proc*N, MPI_DOUBLE, receive_from, receive_from, MPI_COMM_WORLD, &status);
      int send_to = rank_above_me(p, my_rank);
      MPI_Send(row, num_per_proc*N, MPI_DOUBLE, send_to, my_rank, MPI_COMM_WORLD);
    }
    for (r = 0; r < num_per_proc; r++) {
      for (c = 0; c < num_per_proc; c++) {
        for (i = 0; i < N; i++) {
          result += row[r][i]*col[c][i];
        }
        my_c_row[round+r][c] = result;
      }
    }
    round += num_per_proc;
  }

  double result_cols[N][N];
  // gather results
  MPI_Gather(my_c_row, num_per_proc*N, MPI_DOUBLE, result_cols, num_per_proc*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  end = MPI_Wtime();
  elapsed = end - start;
  //if (my_rank == 0) {
    printf("%d : Took %f\n", my_rank, elapsed);
  //}

  MPI_Finalize();
}

void populate(double matrix[N][N]) {
  double r;
  srand(time(NULL));
  int i, j;
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      r = (double)rand()/(double)RAND_MAX;
      matrix[i][j] = r;
    }
  }
}

void transpose(double matrix[N][N], double transposed[N][N]) {
  int i, j;
  for (i = 0; i < N; i++) {
    for (j = i; j < N; j++) {
      transposed[i][j] = matrix[j][i];
      transposed[j][i] = matrix[i][j];
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
