#include <stdlib.h>
#include <mpi.h>
#include <stdio.h>

void odd_even_simple_merge(int * arr, int size, int * result) {
  int half_size = size / 2;
  int * a = arr;
  int * b = arr + half_size;

  int a_index = 0, b_index = 0, result_index = 0;
  while ((a_index < half_size) && (b_index < half_size)) {
    if (a[a_index] < b[b_index]) {
      result[result_index] = a[a_index];
      ++result_index;
      a_index += 2;
    }
    else {
      result[result_index] = b[b_index];
      ++result_index;
      b_index += 2;
    }
  }

  int * remaining;
  int remaining_index, remaining_max;
  if (a_index == half_size) {
    remaining = b;
    remaining_index = b_index;
    remaining_max = half_size;
  }
  else {
    remaining = a;
    remaining_index = a_index;
    remaining_max = half_size;
  }

  for (; remaining_index < remaining_max; ++result_index, remaining_index += 2)
    result[result_index] = remaining[remaining_index];
}

void compare_exchange(int * a, int * b) {
  if (*a > *b) {
    int tmp = *a;
    *a = *b;
    *b = tmp;
  }
}

void odd_even_merger(int * arr, int size) {
    int * left_arr = malloc((size >> 1) * sizeof(int));
    int * right_arr = malloc((size >> 1) * sizeof(int));
    odd_even_simple_merge(arr, size, left_arr);
    odd_even_simple_merge(arr + 1, size, right_arr);
    for (int i = 0; i < size >> 1; ++i) {
            if (left_arr[i] < right_arr[i]) {
            arr[2 * i] = left_arr[i];
            arr[2 * i + 1] = right_arr[i];
            }
            else
            {
                arr[2 * i] = right_arr[i];
                arr[2 * i + 1] = left_arr[i];
            }
        }

    for (int i = 0; i < size - 1; ++i)
        compare_exchange(arr + i, arr + i + 1);

    free(left_arr);
    free(right_arr);
}

void least_significant_digit_sort(int * arr, int size) {
    int bits = 8;
    int ** numbers = malloc(2 * sizeof(int*));
    numbers[0] = malloc(size * sizeof(int));
    numbers[1] = malloc(size * sizeof(int));
    int counters[] = { 0, 0 }; 
    for (int bit_num = 0; bit_num < bits; ++bit_num) {
        for (int i = 0; i < size; ++i) {
            int bit = ((arr[i] >> bit_num) & 1);
            numbers[bit][counters[bit]] = arr[i];
            ++counters[bit];
        }
        int k = 0;
        for (int i = 0; i < counters[0]; ++i, ++k)
            arr[k] = numbers[0][i];
        for (int i = 0; i < counters[1]; ++i, ++k)
            arr[k] = numbers[1][i];

        counters[0] = counters[1] = 0;
  }

  free(numbers[0]);
  free(numbers[1]);
  free(numbers);
}

void radix_batcher_sort(int * arr, int size)
{
    int length = size >> 1;

    MPI_Init(NULL, NULL);
    int pid, np;

    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    if (pid == 0)
    {
        MPI_Status status;
        MPI_Send(arr + length, length, MPI_INT,
                    1, 0, MPI_COMM_WORLD);
        least_significant_digit_sort(arr, length);
        MPI_Recv(arr + length, length, MPI_INT,
            1, 0,
            MPI_COMM_WORLD,
            &status);
    }
    if (pid == 1)
    {
        int * right_array = malloc(length * sizeof(int));
        MPI_Status status;
        MPI_Recv(right_array, length, MPI_INT,
            0, 0,
            MPI_COMM_WORLD,
            &status);
        least_significant_digit_sort(right_array, length);
        MPI_Send(right_array, length, MPI_INT,
                    0, 0, MPI_COMM_WORLD);
    }
    odd_even_merger(arr, size);
}


int main()
{
    int a[] = {5, 2, 3, 6, 8, 1, 0, 4};
    radix_batcher_sort(a, 8);
    int pid;
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    if (pid == 0)
    {
        for(int i = 0; i < 8; ++i)
            printf("%d ", a[i]);
        printf("\n");
    }
    MPI_Finalize();
}

