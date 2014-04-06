#include <stdio.h>
#include "util.h"

int main(int argc, char const *argv[])
{
	double **A, **B, **C;
	int i,j,k;
	int np = 3;

	allocateMatrix(np,&A);
	allocateMatrix(np,&B);
	allocateMatrix(np,&C);

	for (i = 0; i < np; i++){
	    for (j = 0; j < np; j++){
	    	A[i][j] = 1;
	    }
	}

	for (i = 0; i < np; i++){
	    for (j = 0; j < np; j++){
	    	B[i][j] = 1;
	    }
	}

	printMatrix(np,A);
	printMatrix(np,B);

	for (i = 0; i < np; i++){
	    for (j = 0; j < np; j++){
	      for (k = 0; k < np; k++)
	        C[i][j] += A[i][k]*B[k][j];
	    }
	  }
	
	printMatrix(np,C);

	freeMatrix(np,A);
	freeMatrix(np,B);
	freeMatrix(np,C);

	return 0;
}