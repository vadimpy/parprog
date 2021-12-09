/******************************************************************************
* ЗАДАНИЕ: bugged2.c
* ОПИСАНИЕ:
*   Еще одна программа на OpenMP с багом. 
******************************************************************************/

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
    int nthreads, tid;
    float total = 0.0;
    int i;

    #pragma omp parallel private(tid) shared(total)
    {
        tid = omp_get_thread_num();
        if (tid == 0)
        {
            nthreads = omp_get_num_threads();
            printf("Number of threads = %d\n", nthreads);
        }

        #pragma omp barrier
        printf("Thread %d is starting...\n", tid);
        #pragma omp for reduction(+:total) private(i) schedule(static, 10)
        for (i = 0; i < 10000000; i++)
            total += i*1.0;

        #pragma omp barrier
        printf ("Thread %d is done! Total= %e\n", tid, total);
    }
}
