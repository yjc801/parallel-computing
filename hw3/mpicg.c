#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define CONST_TOL	0.00001
#define MPI_RANK_ROOT	0
#define MAX_ITER	1000

void readFile(char *fileName, int *n, double **value, int **colind, int **rbegin, double **b){
	FILE *fp=fopen(fileName,"r");
	char buf[200];
	int i,j,k,l;
	int p,q;
	double w;
	int m;
	if ((fp=fopen(fileName,"r"))==NULL){
	  fprintf(stderr,"Open file errer, check filename please.\n");
	}
	fgets(buf,200,fp);
	*n=atoi(buf);
	l=0;
	while(buf[l++]!=' ');
	m=atoi(&buf[l]);
	printf("Matrix size: %d, #Non-zeros: %d\n",*n,m);
	(*value)=(double*)malloc(sizeof(double)*m);
	(*colind)=(int*)malloc(sizeof(int)*m);
	(*rbegin)=(int*)malloc(sizeof(int)*((*n)+1));
	(*b)=(double*)malloc(sizeof(double)*(*n));
	for(i=0;i<(*n);i++){(*b)[i]=1.0;}
	k=-1;
	for(i=0;i<m;i++){
	  fgets(buf,200,fp);
	  l=0;p=atoi(&buf[l]);
	  while(buf[l++]!=' '); q=atoi(&buf[l]);
	  while(buf[l++]!=' '); w=atof(&buf[l]);
	  (*value)[i]=w;
	  (*colind)[i]=q;
	  if(p!=k){
	    k=p;
	    (*rbegin)[p]=i;
	  }
	}
	(*rbegin)[*n]=m;
	fclose(fp);
}
void scatterData(int *n, int *m, double **value, int **colind, int **rbegin, double **b, int rank, int nproc){
	int np;
	int *sendcnts;
	int *displs;
	double *gvalue;
	int *gcolind;
	int *grbegin;
	double *gb;
	int i,j;
	if (rank==MPI_RANK_ROOT){
	  sendcnts=(int*)malloc(sizeof(int)*(nproc));
	  displs=(int*)malloc(sizeof(int)*(nproc));
	  np=(*n)/nproc;
	  gvalue=(*value);(*value)=NULL;
	  gcolind=(*colind);(*colind)=NULL;
	  grbegin=(*rbegin);(*rbegin)=NULL;
	  gb=(*b);(*b)=NULL;
	  for(i=0;i<nproc;i++){
	    sendcnts[i]=0;
	    for(j=i*np;j<(i+1)*np;j++){
	      sendcnts[i]+=grbegin[j+1]-grbegin[j];
	    }//could be simplified
	  }
	  displs[0]=0;
	  for(i=1;i<nproc;i++){
	    displs[i]=displs[i-1]+sendcnts[i-1];
	  }
	}
	fflush(stdout);
	MPI_Bcast(n,1,MPI_INT,MPI_RANK_ROOT,MPI_COMM_WORLD);
	MPI_Scatter(sendcnts,1,MPI_INT,m,1,MPI_INT,MPI_RANK_ROOT,MPI_COMM_WORLD);
	np=(*n)/nproc;
	printf("Process %d get n=%d, m=%d\n",rank,*n,*m);
	fflush(stdout);
	(*value)=(double*)malloc(sizeof(double)*(*m));
	(*colind)=(int*)malloc(sizeof(int)*(*m));
	(*rbegin)=(int*)malloc(sizeof(int)*(np+1));
	(*b)=(double*)malloc(sizeof(double)*np);
	MPI_Scatterv(gvalue, sendcnts, displs, MPI_DOUBLE, (*value), (*m),MPI_DOUBLE, MPI_RANK_ROOT, MPI_COMM_WORLD);
	MPI_Scatterv(gcolind, sendcnts, displs, MPI_INT, (*colind), (*m),MPI_INT, MPI_RANK_ROOT, MPI_COMM_WORLD);
	
	MPI_Scatter(grbegin, np, MPI_INT, (*rbegin), np, MPI_INT, MPI_RANK_ROOT, MPI_COMM_WORLD);

	MPI_Scatter(gb, np, MPI_DOUBLE, (*b), np, MPI_DOUBLE, MPI_RANK_ROOT, MPI_COMM_WORLD);

	int offset=(*rbegin)[0];
	for(i=0;i<np;i++){
	  (*rbegin)[i]-=offset;
	}
	(*rbegin)[np]=(*m);
	if (rank==MPI_RANK_ROOT){
	  free(gvalue);
	  free(gcolind);
	  free(grbegin);
	  free(gb);
	  free(sendcnts);
	  free(displs);
	}

}
void writeFile(int n,double *answer){
	FILE *fp=fopen("output.dat","w");
	int i;
	for(i=0;i<n;i++){
	  fprintf(fp,"%.10f\n",answer[i]);
	}
	fclose(fp);
}
double* cg(int n, double *value, int* colind, int* rbegin, double *b, int rank, int nproc){

	int i; 				// counter
	int k, k1, k2, j; 	// index
	int iter;
	double alpha, beta;
	double tmp; 		// temporary var
	int local_n;
	double local_sum_r,local_sum_pq, local_sum_r_next;	// local sum
	double sum_r, sum_pq, sum_r_next; // global sum
	double *local_r;
	double *local_r_next;
	double *local_p;
	double *local_x;
	double *local_q;
	double *p;
	double *x;

	// Initialization
	local_n = n/nproc;
	
	// local arrays
	local_r = 		(double *) malloc(local_n*sizeof(double));
	local_r_next = 	(double *) malloc(local_n*sizeof(double));
	local_p = 		(double *) malloc(local_n*sizeof(double));
	local_q = 		(double *) malloc(local_n*sizeof(double));
	local_x = 		(double *) malloc(local_n*sizeof(double));

	// global arrays
	p = 			(double *) malloc(n*sizeof(double));
	x = 			(double *) malloc(n*sizeof(double));

	iter = 0;

	for (i = 0; i < local_n; i++){
		local_x[i] = 0;
	}
	for (i = 0; i < local_n; i++){
		local_p[i] = b[i];
		local_r[i] = b[i];
	}
		
	MPI_Allgather(local_p,local_n,MPI_DOUBLE,p,local_n,MPI_DOUBLE,MPI_COMM_WORLD);

	// Iteration
	while (iter < MAX_ITER){

		// Sparse matrix multiplication
		for (i = 0; i < local_n; i++){
	  		tmp = 0;
	  		k1 = rbegin[i];
	  		k2 = rbegin[i+1]-1;
	  		if (k2 < k1){
	  		 	continue;
	  		}
	  		for (k = k1; k <= k2; k++){
	  			j = colind[k];
	  			tmp += value[k]*p[j];
		  	}
	  		local_q[i] = tmp;
		}
  		local_sum_r = 0;
  		local_sum_pq = 0;

  		for (i = 0; i < local_n; i++){
  			local_sum_r += pow(local_r[i],2);
  			local_sum_pq += local_p[i]*local_q[i];
  	}

  		MPI_Allreduce(&local_sum_r, &sum_r, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
  		MPI_Allreduce(&local_sum_pq, &sum_pq, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 

  		alpha = sum_r/sum_pq;

  if(rank==0)		printf("Iteration %d: sum_r is %e\n", iter,sum_r);

  		for (i = 0; i < local_n; i++){
  			local_x[i]+=alpha*local_p[i];
  			local_r_next[i] = local_r[i] - alpha*local_q[i];
  			local_r[i] = local_r_next[i];
  		}
  		
  		local_sum_r_next = 0;
  		 for (i = 0; i < local_n; i++){
  			local_sum_r_next += pow(local_r_next[i],2);
 		}
  		
  		MPI_Allreduce(&local_sum_r_next, &sum_r_next, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
		if (sum_r_next < CONST_TOL) break;
  		beta = sum_r_next/sum_r;

   		for (i = 0; i < local_n; i++){
   			local_p[i] = local_r[i]+beta*local_p[i];
   		}


  		MPI_Allgather(local_p,local_n,MPI_DOUBLE,p,local_n,MPI_DOUBLE,MPI_COMM_WORLD);
  		iter = iter + 1;
	} // end while

  	 MPI_Allgather(local_x,local_n,MPI_DOUBLE,x,local_n,MPI_DOUBLE,MPI_COMM_WORLD);
	
	if (rank == 0) printf("Total # of iteration = %d\n",iter);
	
	free(local_p);
	free(local_q);
	free(local_r);
	free(local_r_next);
	free(local_x);
	free(p);
	
	return x;
}


int main(int argc, char* argv[]){
  int n;
  int m;
  double *value=NULL;
  int *colind=NULL;
  int *rbegin=NULL;
  double *answer=NULL;
  double *b=NULL;

  int nproc,rank,namelen;
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&nproc);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Get_processor_name(processor_name,&namelen);
  printf("Process %d on %s out of %d\n",rank, processor_name, nproc);
  fflush(stdout);
  if (rank==MPI_RANK_ROOT){
    readFile(argv[1],&n,&value,&colind,&rbegin,&b);
    scatterData(&n,&m,&value,&colind,&rbegin,&b,rank,nproc);
  }
  else{
    scatterData(&n,&m,&value,&colind,&rbegin,&b,rank,nproc);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  double tv1,tv2;
  tv1=MPI_Wtime();
  answer=cg(n,value,colind,rbegin,b,rank,nproc);
  tv2=MPI_Wtime();
  if (rank==MPI_RANK_ROOT){
    printf("Process %d takes %.10f seconds\n",rank,tv2-tv1);
  }
  if (value!=NULL) {free(value);}
  if (colind!=NULL) {free(colind);}
  if (rbegin!=NULL) {free(rbegin);}
  if (b!=NULL) {free(b);}
  if (rank==MPI_RANK_ROOT) {
    writeFile(n,answer);
  }
  if (answer!=NULL) {free(answer);}
  MPI_Finalize();
  return 1;  
}
