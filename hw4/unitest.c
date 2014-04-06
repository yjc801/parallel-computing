#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ITS 1

void MPItest(int rank, int nproc, int size, double *buf, double *buf2){
/*
  Implement the send-recv test, uni/bi-direction ring, MPI_scan in this function.
  Submit File Names:
  send-recv		:	sendtest.c
  uni-direction ring	:	unitest.c
  bi-direction ring	: 	bitest.c
  MPI_Scan		:	scantest.c
*/
    if (rank == 0) {
      MPI_Send(buf2,size,MPI_DOUBLE,rank+1,0,MPI_COMM_WORLD);
      MPI_Recv(buf,size,MPI_DOUBLE,nproc-1,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      printf("%s\n", buf[1]);   
    }
    else if (rank == nproc - 1) {
      MPI_Recv(buf,size,MPI_DOUBLE,rank-1,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      printf("%s\n", buf[1]); 
      MPI_Send(buf,size,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
    }else{
      MPI_Recv(buf,size,MPI_DOUBLE,rank-1,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      printf("%s\n", buf[1]);
      MPI_Send(buf,size,MPI_DOUBLE,rank+1,0,MPI_COMM_WORLD); 
  }
}


int main(int argc, char **argv)
{
  int    i, rank, nproc;
  double startTime;
  double t[ITS];
  char   procName[200];
  int	 procNameLen;
  int    size;
  double *buf=NULL;
  double *buf2=NULL;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&nproc);
  MPI_Get_processor_name(procName,&procNameLen);
  printf("Proc %d on %s\n",rank,procName);
  if (argc!=2){
    printf("usage: ./MPItest size\n");
  }
  else{
    size=4096; //2^12
  }

  size=atoi(argv[1]);   //message length, in size of doubles

  //allocate two buffers
  //In some cases, you may only use one of them
  buf=(double*)malloc(sizeof(double)*size);
  buf2=(double*)malloc(sizeof(double)*size); 
  for(i=0;i<size;i++){
    buf[i]=i; buf2[i]=-i;
  }

  //Run the test ITS times
  for (i=0; i<ITS; i++) {
    MPI_Barrier(MPI_COMM_WORLD);
    startTime = MPI_Wtime();
    MPItest(rank, nproc, size, buf, buf2);
    MPI_Barrier(MPI_COMM_WORLD);
    t[i] = MPI_Wtime() - startTime;
  }
  if (buf!=NULL) {free(buf);buf=NULL;}
  if (buf2!=NULL) {free(buf2);buf2=NULL;}

  //Compute the average time and the standard deviation
  if (rank == 0) {
    double avg=0;
    double std=0;
    for (i=0; i<ITS; i++) {
      avg += t[i];
    }
    avg=avg/ITS;
    printf("conducted %i trials, numprocs = %i \n", ITS, nproc);
    for (i=0;i<ITS;i++){
      std+=(t[i]-avg)*(t[i]-avg);
    }
    std=sqrt(std/ITS);
    printf("average time = %e\n", avg);
    printf("standard deviation = %e\n", std);
  }

  MPI_Finalize();
  return(0);
}
