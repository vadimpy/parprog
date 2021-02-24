#include "mpi.h"
#include "euler.c"
#include "stdio.h"

#define N 4

double e_foo(int from, int to)
{
    long long fact = 1;

    for (int i = 1; i <= from; ++i)
        fact *= i;

    int i;
    double exp = 0;
    for(i = from; i <= to; ++i)
    {
        exp += 1.0/fact;
        fact *= (i + 1);
    }
    return exp;
}

int main()
{
    MPI_Init(NULL, NULL);
    int pid, np;

    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    if (pid == 0)
    {
        double tmp;
        MPI_Status status;
        float eu = 0;
        for (int i = 1; i < np; i++) {
            MPI_Recv(&tmp, 1, MPI_DOUBLE,
                     MPI_ANY_SOURCE, 0,
                     MPI_COMM_WORLD,
                     &status);
            eu += tmp;
        }
        printf("e = %.12f\n", eu);
    }
    else
    {
        int from = (pid - 1) * N;
        int to = pid * N - 1;
        double part_sum = e_foo(from, to);
        printf("%.12f\n", part_sum);
        MPI_Send(&part_sum, 1, MPI_DOUBLE,
                 0, 0, MPI_COMM_WORLD);
    }
    MPI_Finalize();
}

