#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <malloc.h>

//#define N 8

struct Vector {
  double x;
  double y;
  double z;
};

struct Particle {
  struct Vector position;
  int charge;
};

int rank_above_me(int, int);
int rank_below_me(int, int);
struct Vector subtract_vectors(struct Vector, struct Vector);
struct Vector add_vectors(struct Vector, struct Vector);
struct Vector calc_force_vector(struct Particle, struct Particle);
double calc_magnitude(struct Vector);
struct Vector scalar_multiply(struct Vector, double);
double gen_coordinate();
int gen_charge();
struct Vector calc_total_force(struct Particle *, int, int);
int main(int argc, char* argv[]) {
  int N = atoi(argv[1]);
  if (N == 0) return;
  //printf("Placing %d particles.\n", N);

  double start, end, elapsed;
  int my_rank, p;

  int charges[6] = {-3, -2, -1, 1, 2, 3};

  MPI_Status status;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  srand(my_rank + time(NULL)); // to ensure each process gets a different seed

  const int nitems = 3;
  int blocklengths[3] = {1,1,1};
  MPI_Datatype types[3] = { MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE };
  MPI_Datatype mpi_vector_type;
  MPI_Aint offsets[3];
  offsets[0] = offsetof(struct Vector, x);
  offsets[1] = offsetof(struct Vector, y);
  offsets[2] = offsetof(struct Vector, z);
  MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_vector_type);
  MPI_Type_commit(&mpi_vector_type);

  int num_particles = N/p;

  // 10 x 10 x 10 box, generate x,y,z (0,10]
  // each x,y,z generate Q {-3,-2,-1,1,2,3}
  struct Particle *particles = malloc(num_particles * sizeof(struct Particle));
  int i;
  for (i = 0; i < num_particles; i++) {
    struct Vector point = { gen_coordinate(),
      gen_coordinate(),
      gen_coordinate() };
    struct Particle particle = {
      point,
      gen_charge()
    };
    particles[i] = particle;
  }

  /*
  // print the particles this process has
  for (i = 0; i < num_particles; i++) {
    printf("%f %f %f, %d\n", particles[i].position.x,
        particles[i].position.y,
        particles[i].position.z,
        particles[i].charge);
  }
  */

  start = MPI_Wtime();

  struct Vector *forces = malloc(num_particles * sizeof(struct Vector));
  for (i = 0; i < num_particles; i++) {
    struct Vector total = calc_total_force(particles, num_particles, i);
    forces[i] = total;
  }

  struct Vector results[N];
  MPI_Gather(forces, num_particles, mpi_vector_type, results, num_particles, mpi_vector_type, 0, MPI_COMM_WORLD);

  /*
  // print all results
  if (my_rank == 0) {
    for (i = 0; i < N; i++) {
      printf("particle %d: %f %f %f\n", i, results[i].x, results[i].y, results[i].z);
    }
  }
  */
  
  MPI_Type_free(&mpi_vector_type);

  end = MPI_Wtime();

  MPI_Finalize();

  elapsed = end - start;
  printf("%d took %f\n", my_rank, elapsed);
}

double gen_coordinate() {
  return (double)rand()/(double)RAND_MAX * 10.0;
}

int gen_charge() {
  int charges[6] = {-3, -2, -1, 1, 2, 3};
  int index = (int)((double)rand()/(double)RAND_MAX * 6)%6;
  return charges[index];
}

struct Vector add_vectors(struct Vector i, struct Vector j) {
  struct Vector add = { i.x + j.x, i.y + j.y, i.z + j.z };
  return add;
}

struct Vector subtract_vectors(struct Vector u, struct Vector v) {
 struct Vector difference = { u.x - v.x, u.y - v.y, u.z - v.z };
 return difference;
}

double calc_magnitude(struct Vector point) {
  return sqrt(point.x*point.x + point.y*point.y + point.z*point.z);
}

struct Vector scalar_multiply(struct Vector v, double k) {
  struct Vector result = { v.x*k, v.y*k, v.z*k};
  return result;
}

struct Vector calc_force_vector(struct Particle i, struct Particle j) {
  struct Vector difference = subtract_vectors(i.position, j.position);
  double multiplier = (double)(i.charge * j.charge) /
    pow(calc_magnitude(difference), 3);
  return scalar_multiply(difference, multiplier);
}

struct Vector calc_total_force(struct Particle *particles,
    int num_particles,
    int index) {
  struct Vector result = { 0.0, 0.0, 0.0 };
  int j;
  for (j = 0; j < num_particles; j++) {
    if (j == index) continue;
    struct Vector force = calc_force_vector(particles[index], particles[j]);
    result = add_vectors(result, force);
  }
  return result;
}

