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
  //double x, y, z;
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
struct Vector calc_total_force(struct Particle *, int, struct Particle);
void array_copy(struct Particle *, struct Particle *, int);

int main(int argc, char* argv[]) {
  int N = atoi(argv[1]);
  if (N == 0) return 0;
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
  MPI_Type_struct(nitems, blocklengths, offsets, types, &mpi_vector_type);
  MPI_Type_commit(&mpi_vector_type);

  const int oitems = 2;
  int oblocklengths[2] = {1, 1};
  MPI_Datatype otypes[2] = {mpi_vector_type, MPI_INT};
  MPI_Datatype mpi_particle_type;
  MPI_Aint ooffsets[2] = { 0, 24 };
  offsets[0] = offsetof(struct Particle, position);
  offsets[1] = offsetof(struct Particle, charge);
  MPI_Type_struct(oitems, oblocklengths, ooffsets, otypes, &mpi_particle_type);
  MPI_Type_commit(&mpi_particle_type);
  
  /*
  const int oitems = 2;
  int oblocklengths[2] = {3, 1};
  MPI_Datatype otypes[2] = {MPI_DOUBLE, MPI_INT};
  MPI_Datatype mpi_particle_type;
  MPI_Aint ooffsets[2];
  offsets[0] = offsetof(struct Particle, position);
  offsets[1] = offsetof(struct Particle, charge);
  MPI_Type_struct(oitems, oblocklengths, ooffsets, otypes, &mpi_particle_type);
  MPI_Type_commit(&mpi_particle_type);
  */

  int num_particles = N/p;

  // 10 x 10 x 10 box, generate x,y,z (0,10]
  // each x,y,z generate Q {-3,-2,-1,1,2,3}
  struct Particle *my_particles = malloc(num_particles * sizeof(struct Particle));
  struct Particle *other_particles = malloc(num_particles * sizeof(struct Particle));

  int i;
  for (i = 0; i < num_particles; i++) {
    struct Vector point = { gen_coordinate(),
      gen_coordinate(),
      gen_coordinate() };
    struct Particle particle = {
      point,
      gen_charge()
    };
    my_particles[i] = particle;
  }

  array_copy(my_particles, other_particles, num_particles);

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

  // for all the N on all processors, sum to the total;

  // array of forces on the particle that I am taking care of
  // this is slowly summed up as I retrieve more elements
  struct Vector *forces = malloc(num_particles * sizeof(struct Vector));
  int round;
  for (round = 0; round < p; round++) {
    int i;
    for (i = 0; i < num_particles; i++) {
      struct Vector total = calc_total_force(other_particles, num_particles, my_particles[i]);
      forces[i] = add_vectors(forces[i], total);
    }
    if (round == p-1) break;

    // send and receive
    if (my_rank % 2 == 0) {
      int send_to = rank_above_me(p, my_rank);
      MPI_Send(other_particles, num_particles, mpi_particle_type, send_to, my_rank, MPI_COMM_WORLD);
      int receive_from = rank_below_me(p, my_rank);
      MPI_Recv(other_particles, num_particles, mpi_particle_type, receive_from, receive_from, MPI_COMM_WORLD, &status);
    } else {
      int receive_from = rank_below_me(p, my_rank);
      MPI_Status status;
      struct Particle buffer[num_particles]; // need to store in buffer else it will overwrite what we send
      MPI_Recv(buffer, num_particles, mpi_particle_type, receive_from, receive_from, MPI_COMM_WORLD, &status);
      int x;
      int send_to = rank_above_me(p, my_rank);
      MPI_Send(other_particles, num_particles, mpi_particle_type, send_to, my_rank, MPI_COMM_WORLD);
      // copy numbers from buffer to row
      int i;
      for (i = 0; i < num_particles; i++) {
        other_particles[i] = buffer[i];
      }
    }
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
  MPI_Type_free(&mpi_particle_type);

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
  struct Vector difference;
  if (i.position.x == j.position.x &&
      i.position.y == j.position.y &&
      i.position.z == j.position.z) return difference;
  difference = subtract_vectors(i.position, j.position);
  double multiplier = (double)(i.charge * j.charge) /
    pow(calc_magnitude(difference), 3);
  return scalar_multiply(difference, multiplier);
}

/*
 * Calculate the force that a group of particles exert
 * one one particular particle.
 */
struct Vector calc_total_force(struct Particle *particles,
    int num_particles, struct Particle particle) {
  struct Vector result = { 0.0, 0.0, 0.0 };
  int j;
  for (j = 0; j < num_particles; j++) {
    struct Vector force = calc_force_vector(particle, particles[j]);
    result = add_vectors(result, force);
  }
  return result;
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

void array_copy(struct Particle *from, struct Particle *to, int size) {
  int i;
  for (i = 0 ; i < size; i++) {
    struct Vector v = { from[i].position.x,
      from[i].position.y,
      from[i].position.z };
    struct Particle p = { v, from[i].charge };
    to[i] = p;
  }
}

