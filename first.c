#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

int main(int argc, char* argv[]) {
	int my_rank;
	int p;
	double start_time, end_time, elapse_time;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	if (my_rank == 0) {
		int start = 0;
		printf("hi from %d\n", my_rank);
		MPI_Send(&start, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
	} else {
		int received_from = my_rank - 1;
		MPI_Status status;
		MPI_Recv(&received_from, 1, MPI_INT, my_rank - 1, 0, MPI_COMM_WORLD, &status);
		printf("hi from %d\n", my_rank);
		if (my_rank + 1 == p) {
		} else {
			int send_to = my_rank + 1;
			MPI_Send(&send_to, 1, MPI_INT, send_to, 0, MPI_COMM_WORLD);
		}
	}
	MPI_Finalize();
}

