#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#define CONST_ALPHA	0.85
#define CONST_TOL	0.00001

void readfile(char *fileName, int *n, double **value, int **colind, int **rbegin){
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
	close(fp);
}


double* pageRank(int n, double *value, int* colind, int *rbegin){

  // Initialization
  double* x;		//pagerank vector
  int i; 			// counter
  int iter;			// iteration
  double y[n];		// temp array to store pagerank vector
  double a;			// alpha
  double maxerr;	// max error
  double total;     // summation of vector y
  double gamma;     // coefficient

	x = (double *)malloc(sizeof y);

	for (i = 0; i < n; i++){
		x[i] = 1.0/n;
	}

	iter = 0;
	a = CONST_ALPHA;
	maxerr = 1.0;
	double thrhd = CONST_TOL;

	while (maxerr > thrhd) {   //termination condition
		maxerr = 0;
		total = 0;
		iter+=1;

#pragma omp parallel shared(n,a,maxerr) reduction(+:total)
{

		int k; 			//counteer
		int k1,k2,j; 	//index
		double prev;// x[i] before updating
		double tmp;	// updated x[i]
	 	double error;	// error between updated value and previous value

#pragma omp for //schedule(static,100)
	  for (i = 0; i < n; i++){
	  		tmp = 0;
	  		k1 = rbegin[i];
	  		k2 = rbegin[i+1]-1;
	  		if (k2 < k1){
	  			continue;
	  		}
	  		for (k = k1; k <= k2; k++){
	  			j = colind[k];
	  			tmp += value[k]*x[j];
	  		}
	  		y[i] = a*tmp;
	  		total += a*tmp;
		}
}

	gamma = 1 - total;

#pragma omp parallel shared(n,gamma) reduction(max:maxerr)
{
	double error, tmp, prev;

	#pragma omp for
		for (i = 0; i < n; i++){ 
	  		tmp = y[i];
	  		prev = x[i];
	  		tmp = tmp + gamma/n;
	  		error = fabs(tmp - prev);
	  		if(error > maxerr){
	  			maxerr = error;
	  		}
	  		x[i] = tmp;
		}// end for loop

}// end parallel

// double sum = 0;
		// for (i = 0; i < n; i++){
			// sum += x[i];
		// }
		// printf("The sum of vector is %e. Max error is %e\n", sum, maxerr);

}// end while loop

	printf("Iteration is %d\n", iter);
	return x;
}

void output(int n, double *answer){
	FILE *fp=fopen("output.dat","w");
	int i;
	for(i=0;i<n;i++){
	  fprintf(fp,"%.16f\n",answer[i]);
	}
	close(fp);
}
int main(int argc, char* argv[]){
	double *value=NULL;
	int *colind=NULL;
	int *rbegin=NULL;
	double *answer=NULL;
	int n;
	struct timeval tv1,tv2;
	readfile(argv[1],&n,&value,&colind,&rbegin);
	gettimeofday(&tv1,NULL);
	answer=pageRank(n,value,colind,rbegin);
	gettimeofday(&tv2,NULL);
	printf("#Threads=%d, Time=%10.5f\n",omp_get_max_threads(),((double)(tv2.tv_sec-tv1.tv_sec)+(double)(tv2.tv_usec-tv1.tv_usec)/1000000));
	output(n,answer);
	if (value!=NULL) {free(value);value=NULL;}
	if (colind!=NULL) {free(colind);colind=NULL;}
	if (rbegin!=NULL) {free(rbegin);rbegin=NULL;}
	if (answer!=NULL) {free(answer);answer=NULL;}
}
