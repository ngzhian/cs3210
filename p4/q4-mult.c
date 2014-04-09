#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

double ras(double, double);
double gen_x1(double, double, double, double);
double gen_x2(double, double);
void get_range(int, int, double *, double *, double, double);
double prob_accept(double new_state, double old_state);
void do_single_run(double *m_ras, double *m_x1, double *m_x2, double *low, double *upp, double);

int main(int argc, char *argv[]) {
  int n = 2,
      my_rank,
      num_procs;

  double start, end,
         range_start, range_end,
         lower_bound = -5.12, upper_bound = 5.12;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  srand(time(NULL) + my_rank);

  double global_min_ras = 100.0,
         global_min_x1, global_min_x2;
  double run_min_ras = 100.0,
         run_min_x1, run_min_x2;

  int runs, root = 0;
  double eps = 1.0;
  for (runs = 0; runs < 100; runs++) {
    // each processor finds its local minimum
    do_single_run(&run_min_ras, &run_min_x1, &run_min_x2, &lower_bound, &upper_bound, eps);
    // reduce and get the lowest minimum out of all processors
    struct {
      double val;
      int rank;
    } min_ras_and_rank = {run_min_ras, my_rank}, run_min;
    // get the best value for a run
    MPI_Reduce(&min_ras_and_rank, &run_min, 1, MPI_DOUBLE_INT, MPI_MINLOC, 0, MPI_COMM_WORLD);
    
    root = run_min.rank;
    MPI_Bcast(&run_min_x1, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);

    lower_bound = (run_min_x1 - eps < lower_bound) ? lower_bound : run_min_x1 - eps;
    upper_bound = (run_min_x1 + eps > upper_bound) ? upper_bound : run_min_x1 + eps;
    eps /= 10;
  }
  if (my_rank == 0) {
    printf("%f\n", run_min_ras);
  }



  MPI_Finalize();
}

/*
 * Runs a single attempt to search for the minimum point.
 * m_ras is where the min value is stored, m_x1 and m_x2 are the coordiates for that value
 * low is the lower bound of x1, and upp the upper bound.
 * eps is the max allowed deviation from the current point, used for finding next point
 */
void do_single_run(double *m_ras, double *m_x1, double *m_x2, double *low, double *upp, double eps) {
  int my_rank, num_procs;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

  double x1_start,  // start x1 value for this processor
         x1_end,    // max x1 value for this processor
         x2_min = -5.12, x2_max = 5.12;

  // partition the problem space to processors
  get_range(my_rank, num_procs, &x1_start, &x1_end, *low, *upp);

  //printf("%d: x1 start %f, x1 end %f\n", my_rank, x1_start, x1_end);
  int tries = 0;
  double cur_x1 = x1_start,
         cur_x2 = *low,
         cur_ras = ras(cur_x1, cur_x2),
         next_x1,
         next_x2,
         next_ras,
         low_x1 = cur_x1,
         low_x2 = cur_x2,
         low_ras = cur_ras;
  while (tries++ < 1000) {
      if (fabs(low_ras < 0.00001)) break;
      if (cur_ras < low_ras) {
        low_ras = cur_ras;
        low_x1 = cur_x1;
        low_x2 = cur_x2;
      }

      next_x1 = gen_x1(x1_start, x1_end, cur_x1, eps);
      next_x2 = gen_x2(cur_x2, eps);

      double next_ras = ras(next_x1, next_x2);
      double prob = prob_accept(next_ras, cur_ras);
      double ran = (double) rand() / (double) RAND_MAX;
      if (ran > prob) {
      } else {
        cur_x1 = next_x1;
        cur_x2 = next_x2;
        cur_ras = next_ras;
      }
      tries++;
  }
  *m_ras = low_ras;
  *m_x1 = low_x1;
  *m_x2 = low_x2;
}

/*
 * Places in s and e the start and end x1 values for processor of rank
 * l and u are the lower and upper bounds of x1 respectively
 */
void get_range(int rank, int num_procs, double *s, double *e, double l, double u) {
  double total_range = u - l;
  double my_range = total_range / num_procs;
  // have an overlap area for the range_start, unless it is outside of domain
  double range_start = l + (rank * my_range);
  // have an overlap area for the range_end, unless it is outside of domain
  double range_end = range_start + my_range;
  *s = range_start;
  *e = range_end;
}

/*
 * Make a new_x1 such that
 * cur_x1 - eps <= new_x1 <= cur_x1 + eps
 * and
 * s <= new_x1 <= e
 */
double gen_x1(double s, double e, double cur_x1, double eps) {
  s = (cur_x1 - eps < s) ? s : cur_x1 - eps;
  e = (cur_x1 + eps > e) ? e : cur_x1 + eps;
  double r = (double) rand() / (double) RAND_MAX;
  double range = e - s;
  double random_double_in_range = r * range;
  return s + random_double_in_range;
}

/*
 * Make a new_x2 such that
 * cur_x2 - eps <= new_x2 <= cur_x2 + eps
 * and
 * -5.12 <= new_x2 <= 5.12
 */
double gen_x2(double cur_x2, double eps) {
  return gen_x1(-5.12, 5.12, cur_x2, eps);
}

/*
 * Compute the value of the Rastigrin function at (x1, x2)
 */
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
