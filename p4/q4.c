#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define EPS 0.10

double ras(double, double);
double random_eps();
void get_start_end_range(double *, double *);
double prob_accept(double new_state, double old_state);
int is_valid_point(double next_1, double next_2, double start, double end);
void pick_random_point_from_cur(double, double, double *, double *);

int main(int argc, char *argv[]) {
  int n = 2,
      my_rank,
      num_procs;

  // time keeping variables
  double start,
         end,
  // ranges that each processor takes care of
         range_start,
         range_end,
  // the domain of Rastrigin function
         lower_bound = -5.12,
         upper_bound = 5.12;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  srand(time(NULL) + my_rank);

  get_start_end_range(&range_start, &range_end);

  int tries = 0;
  double total_range = upper_bound - lower_bound;
  double my_range = total_range / num_procs;
  double cur_x1 = range_start,
         cur_x2 = lower_bound,
         cur_ras = ras(cur_x1, cur_x2),
         next_x1,
         next_x2,
         next_ras,
         low_x1,
         low_x2,
         low_ras = 1000;

  double h = -10, i = -10;

  while (tries < (total_range/EPS * 9 )) {
  //while (tries < 1000) {
    if (cur_ras < low_ras) {
      low_ras = cur_ras;
      low_x1 = cur_x1;
      low_x2 = cur_x2;
    }

    do {
      pick_random_point_from_cur(cur_x1, cur_x2, &next_x1, &next_x2);
    } while (!is_valid_point(next_x1, next_x2, range_start, range_end));
    if (next_x1 > h) {
      h = next_x1;
    }
    if (next_x2 > i) {
      i = next_x2;
    }

    double next_ras = ras(next_x1, next_x2);
    double prob = prob_accept(next_ras, cur_ras);
    //printf("probability is %f, (%.3f, %.3f) -> (%.3f, %.3f)\n",
        //prob, cur_x1, cur_x2, next_x1, next_x2);
    double ran = (double) rand() / (double) RAND_MAX;
    if (ran > prob) {
    } else {
      cur_x1 = next_x1;
      cur_x2 = next_x2;
      cur_ras = next_ras;
    }
    tries++;
  }
  printf("highest x1: %.3f, highest x2: %.3f\n", h, i);
  //start = MPI_Wtime();

  printf("lows: %f at (%f,%f)\n", low_ras, low_x1, low_x2);
  MPI_Finalize();
}

void get_start_end_range(double *start, double *end) {
  int my_rank, num_procs;
  double lower_bound = -5.12,
    upper_bound = 5.12;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  double total_range = upper_bound - lower_bound;
  double my_range = total_range / num_procs;
  double range_start = my_rank * my_range - 5.12;
  double range_end = my_rank * my_range + my_range - 5.12;
  *start = range_start;
  *end = range_end;
  printf("start %f, end %f\n", range_start, range_end);
}


/*
 * Randomly pick a point that is eps distance away from
 * the curent point, in the x1 and/or x2 axis.
 */
void pick_random_point_from_cur(double cur_x1,
    double cur_x2,
    double *next_x1,
    double *next_x2) {
  double eps_x1, eps_x2;
    //eps_x1 = random_eps();
    //eps_x2 = random_eps();
    eps_x1 = 0.01;
    eps_x2 = 0.01;
    *next_x1 = cur_x1 + eps_x1;
    *next_x2 = cur_x2 + eps_x2;
}

double random_eps() {
  int r = rand() % 3;
  switch (r) {
    case 0:
      return -EPS;
    case 1:
      return 0.00;
    case 2:
      return EPS;
  }
}

int is_valid_point(double next_1, double next_2, double start, double end) {
  return next_1 >= start &&
    next_1 <= end &&
    next_2 >= -5.12 &&
    next_2 <= 5.12;
}

double ras(double x1, double x2) {
  const n = 2;
  const double pi = 4. * atan(1.);
  int an = 10 * n;
  double summ = an;
  summ += ((x1*x1) - (10 * cos(2 * pi * x1)));
  summ += ((x2*x2) - (10 * cos(2 * pi * x2)));
  return summ;
}

double prob_accept(double next_ras, double cur_ras) {
  const double t = 0.8;
  double prob = fmin(1.0, exp(((next_ras - cur_ras) / -t)));
  return prob;
}
