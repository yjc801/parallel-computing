#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

int main(int argc, char const *argv[])
{
	int	rank, np;
	// MPI_Comm MPI_COMM_WORLD;
	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&np);

	int buf;

	buf = rank;
	printf("rank %d: %d\n", rank,buf);
	MPI_Sendrecv_replace(&buf,np,MPI_INT,(rank+np-1)%np,0,(rank+1)%np,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	// MPI_Barrier(MPI_COMM_WORLD);
	printf("rank %d: %d\n", rank,buf);
	MPI_Finalize();
	return 0;
}
