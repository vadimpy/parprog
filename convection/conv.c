#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define T 9
#define L 10.0

float64_t T_0(float64_t x) {
    return sin(2 * M_PI * x / L);
}

int main(int argc, char** argv) {
    float64_t dx = 0.001;
    float64_t tau = 0.1;

    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    struct timespec start_time, end_time;
    clock_gettime(CLOCK_REALTIME, &start_time);

    int N = L / dx;
    int length = ceil((float64_t) N / world_size);
    int begin = world_rank * length;
    int end = (world_rank + 1) * length - 1;
    if (end > N) end = N;

    length = end - begin + 1;

    float64_t *current_layer_U = (float64_t*) calloc(length, sizeof(float64_t)); // for each process
    int i;
    for(i = 0; i < end - begin + 1; ++i) current_layer_U[i] = T_0(dx * (begin + i));

    int left_process = (world_size + world_rank - 1) % world_size;
    int right_process = (world_rank + 1) % world_size;

    int k;
    for(k = 0; k * tau < T; ++k) {
        float64_t *new_layer_U = (float64_t*) calloc(end - begin + 1, sizeof(float64_t)); 
        float64_t left_U, right_U;

        MPI_Request left_request, right_request;
        MPI_Isend(&current_layer_U[0], 1, MPI_DOUBLE, left_process, 0, MPI_COMM_WORLD, &left_request);
        MPI_Recv(&right_U, 1, MPI_DOUBLE, right_process, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        MPI_Isend(&current_layer_U[length - 1], 1, MPI_DOUBLE, right_process, 0, MPI_COMM_WORLD, &right_request);
        MPI_Recv(&left_U, 1, MPI_DOUBLE, left_process, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        new_layer_U[0] = current_layer_U[0] 
                        - 0.5 * (tau / dx) * (current_layer_U[1] - left_U)
                        + 0.5 * (tau / dx) * (tau / dx) * (current_layer_U[1] - 2 * current_layer_U[0] + left_U);
        for (i = 1; i < length - 1; ++i) {
            new_layer_U[i] = current_layer_U[i] 
                            - 0.5 * (tau / dx) * (current_layer_U[i + 1] - current_layer_U[i - 1])
                            + 0.5 * (tau / dx) * (tau / dx) * (current_layer_U[i + 1] - 2 * current_layer_U[i] + current_layer_U[i - 1]);
	}
        new_layer_U[length - 1] = current_layer_U[length - 1] 
                        - 0.5 * (tau / dx) * (right_U - current_layer_U[length - 2])
                        + 0.5 * (tau / dx) * (tau / dx) * (right_U - 2 * current_layer_U[length - 1] + current_layer_U[length - 2]);

        free(current_layer_U);
	current_layer_U = new_layer_U;
    }   

    float64_t *result_U = (float64_t*) calloc(N, sizeof(float64_t)); 
    MPI_Gather(current_layer_U, length, MPI_DOUBLE, result_U, length, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (world_rank == 0) {
        clock_gettime(CLOCK_REALTIME, &end_time);
        float64_t elapsed = (end_time.tv_sec - start_time.tv_sec) + (end_time.tv_nsec - start_time.tv_nsec) / 1000000000.0;
        printf("%d processes. Elapsed %f seconds.\n", world_size, elapsed);
        FILE *fpt;
        fpt = fopen("result.csv", "w+");
        fprintf(fpt,"x; U\n");
        for (i = 0; i < N; ++i) {
            fprintf(fpt,"%f; %f\n", i * dx, current_layer_U[i]);
        }
        fclose(fpt);
    }

    // Finalize the MPI environment.
    MPI_Finalize();
}
