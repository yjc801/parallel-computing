#ifndef _HW5_UTIL_
#define _HW5_UTIL_

#define CONST_TOL 0.000001

void generateLocalMatrix(int n, int row_rank, int col_rank, int nproc_row,  double **mA, double **mB);


void verifyAnswerFast(int m, int rank, int row_rank, int col_rank, int nproc_row,double **mC);

void verifyAnswerComplete(int m, int rank, int row_rank, int col_rank, int nproc_row,double **mC);

//Free the allocated matrix
void freeMatrix(int n,double **M);

//Allocate memory for a matrix, guaranteed to be  contiguous. 
void allocateMatrix(int n, double ***M);


//Print the matrix
void printMatrix(int n, double **M);
#endif
