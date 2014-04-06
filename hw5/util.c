#include "util.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void printMatrix(int n, double **M){
  int i,j;
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      printf("<%10.5f>",M[i][j]);
    }
    printf("\n");
  }
}
void allocateMatrix(int n, double ***M){
  int i,j;
  double *tmp;
  tmp=(double*)malloc(sizeof(double)*n*n);
  (*M)=(double**)malloc(sizeof(double)*n);
  for(i=0;i<n;i++){
    (*M)[i]=&tmp[i*n];
    for(j=0;j<n;j++){(*M)[i][j]=0.0;}
  }
}

void freeMatrix(int n,double **M){
  int i;
  double *tmp;
  tmp=&(M[0][0]);
  free(M);
  free(tmp);
}

void generateLocalA(int n, unsigned int seed, double **mA){
  int i,j;
  int rand_next;
  //Generate A:
  srand(seed);
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      rand_next=rand();
      mA[i][j]=((double)rand_next)/RAND_MAX;
    }
  }
}

void generateLocalB(int n, unsigned int seed, double **mB){
  int i,j;
  int rand_next;
  //Generate B:
  srand(seed);
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      rand_next=rand();
      mB[j][i]=((double)rand_next)/RAND_MAX;
    }
  }
}

void generateLocalMatrix(int n, int row_rank, int col_rank, int nproc_row,  double **mA, double **mB){
  int i,j;
  unsigned int seed;
  //Generate A;
  seed=12345+row_rank*nproc_row+col_rank;
  generateLocalA(n,seed,mA);
  //Generate B;
  seed=54321+row_rank*nproc_row+col_rank;
  srand(seed);
  generateLocalB(n,seed,mB);
}

//Verify C[0][0] for each block;
void verifyAnswerFast(int n, int rank, int row_rank, int col_rank, int nproc_row, double **mC){
  int k;
  unsigned int seed;
  int rand_next;
  int i,j;
  double *a=(double*)malloc(sizeof(double)*n);
  double *b=(double*)malloc(sizeof(double)*n);
  double s=0;
  for(k=0;k<nproc_row;k++){
    //generate a block of A and B;
    seed=12345+k*nproc_row+col_rank;
    srand(seed);
    for(i=0;i<n;i++){
      rand_next=rand();
      a[i]=((double)rand_next)/RAND_MAX;
    }
    seed=54321+row_rank*nproc_row+k;
    srand(seed);
    for(i=0;i<n;i++){
      rand_next=rand();
      b[i]=((double)rand_next)/RAND_MAX;
    }
    for(i=0;i<n;i++){
      s=s+a[i]*b[i];
    }
  }
  if (fabs(s-mC[0][0])<CONST_TOL){
    fprintf(stderr,"Proc %d: Correct! (Incomplete verify)\n",rank);
  }
  else{
    fprintf(stderr,"Proc %d: Wrong Answer! (Incomplete) verify)\n",rank);
  }
  free(a);
  free(b);
}

//Verify the answer C completely
void verifyAnswerComplete(int n, int rank, int row_rank, int col_rank, int nproc_row,  double **mC){
  double **A=NULL;
  double **B=NULL;
  double **C=NULL;
  int k;
  unsigned int seed;
  int i,j,l;
  allocateMatrix(n,&A);
  allocateMatrix(n,&B);
  allocateMatrix(n,&C);
  for(k=0;k<nproc_row;k++){
    seed=12345+k*nproc_row+col_rank;
    generateLocalA(n,seed,A);
    seed=54321+row_rank*nproc_row+k;
    generateLocalB(n,seed,B);
    for(i=0;i<n;i++){
      for(j=0;j<n;j++){
        for(l=0;l<n;l++){
          C[i][j]+=A[i][l]*B[l][j];
//if (row_rank+col_rank==0){
//printf("A=%10.5f, B=%10.5f,C=%10.5f\n",A[i][j],B[i][j],C[i][j]);
//}
        }
      }
    }
  }
  int flag=1;
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      if (fabs(mC[i][j]-C[i][j])>CONST_TOL){flag=0;}
    }
  }
  if (flag==1){
    fprintf(stderr,"Proc %d: Correct! (Complete verify)\n",rank);
  }
  else{
    fprintf(stderr,"Proc %d: Wrong Answer! (Complete) verify)\n",rank);
  }
  freeMatrix(n,A);
  freeMatrix(n,B);
  freeMatrix(n,C);
}


