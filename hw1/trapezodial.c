/*
Implementation for Reducation Operation
*/

#include <omp.h>
#include <stdio.h>

main() {
    double  integral;   /* Store result in integral   */
    double  a, b;       /* Left and right endpoints   */
    int    n;      /* Number of intervals       */
    double  h;          /* interval  width       */
    double  x;
    int    i;
    int nthreads; 
    double start, end, compute_time; 
    double PI16D = 3.141592653589793;  
    double f(double x);  /* Function we're integrating */

//    int input_thread;
//    printf("Enter the number of threads:\n");
//    scanf("%d",&input_thread);
//    omp_set_num_threads(input_thread);


    a= 0.0; 
    b= 1.0; 
    printf("Enter the number of intervals,  n. (n-1) must be divisible by no. of threads. \n");
    scanf("%d",  &n);

    h = (b-a)/n;
    integral = (f(a) + f(b))/2.0;
    
    nthreads = omp_get_max_threads(); 
    printf("No. of threads = %d \n", nthreads); 
    start = omp_get_wtime(); 

    double local_sum[nthreads]; /* array which is used to store local sum of each thread */
    int pid; /* the thread id */

    for (i = 0; i < nthreads; i++) /* initialize the array*/
    {
    	local_sum[i] = 0;
    }

#pragma omp parallel shared(h,a,n) private(x,pid)
{
    pid = omp_get_thread_num();
#pragma omp for  
    for (i = 1; i < n; i++) { /* (n-1) iterations */
        x = a  + i* h;
        local_sum[pid] += f(x);
    }

#pragma omp master
    for (i = 0; i < nthreads; i++) {
    	// printf("%25.16e\n", local_sum[i]);
		integral+=local_sum[i];
    }
}
    integral = integral* h * 4.0;  /* we calculate pi/4 */

    end = omp_get_wtime();
    compute_time = end - start; 

    printf("With nthreads = %d threads, and n = %d intervals, the error in PI \n",
	   nthreads, n);
    printf(" = %25.16e\n",  PI16D - integral);
    printf("Time taken with implemented reduction operator is %15.2e\n",  compute_time);  

} /* main */


double f(double x) {
    double return_val;
    /* Calculate f(x).  Store calculation in return_val. */

    return_val = 1.0/(1.0 + x*x);
    return return_val;
}/* f */
