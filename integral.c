#define _GNU_SOURCE
#include <sched.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/sysinfo.h>
#include <pthread.h>

#define LEFT_END   0.001
#define RIGHT_END  1
#define STEP        0.0000001

typedef struct thread_args
{
  pthread_t pthread;
  double left_end;
  double right_end;
  double answer;
} thread_args_t;

double integrate(double left_end, double right_end)
{
  double sum = 0;
  for(int i = 0; i < (right_end - left_end)/STEP; i++)
  {
    sum += STEP*sin(i*STEP + left_end)/(i*STEP + left_end);
  }
  return sum;
}

void* thread_func(void* args)
{
  ((thread_args_t*)args)->answer = integrate(((thread_args_t*)args)->left_end, ((thread_args_t*)args)->right_end);
}
void* busy_wait(void* args)
{
  // I'M BORED IN THE HOUSE AND I'M IN THE HOUSE BORED
  int busywaiter = LEFT_END;
  while(1)
  busywaiter*= RIGHT_END;
}

void err_sys(const char* error)
{
	perror(error);
	exit(1);
}

int main(int argc, char* argv[])
{

  if(argc != 2)
    err_sys("INCORRECT INPUT");

  int thread_number = strtol(argv[1], NULL, 10);
  int cpu_number = get_nprocs();
  double answer = 0;

  cpu_set_t cpu_set;

  thread_args_t* thread_args_array = (thread_args_t*)calloc(thread_number, sizeof(thread_args_t));

  for(int i = thread_number; i < cpu_number; ++i)
  {
    CPU_ZERO(&cpu_set);
    CPU_SET(i, &cpu_set);

    thread_args_array[i].pthread = pthread_self();
    pthread_setaffinity_np(thread_args_array[i].pthread, sizeof(cpu_set_t), &cpu_set);
    pthread_create(&(thread_args_array[i].pthread), NULL, busy_wait, NULL);

  }
  for(int i = 0; i < thread_number; ++i)
  {
    thread_args_array[i].left_end = LEFT_END+(RIGHT_END-LEFT_END)/thread_number*i;
    thread_args_array[i].right_end = LEFT_END+(RIGHT_END-LEFT_END)/thread_number*(i+1);
    thread_args_array[i].answer = 0;
    thread_args_array[i].pthread = pthread_self();

    CPU_ZERO(&cpu_set);
    CPU_SET(i%cpu_number ,&cpu_set);

    pthread_setaffinity_np(thread_args_array[i].pthread, sizeof(cpu_set_t), &cpu_set);
    pthread_create(&(thread_args_array[i].pthread), NULL, thread_func, (void*)&thread_args_array[i]);
  }


  for (int i = 0; i < thread_number; ++i)
    pthread_join(thread_args_array[i].pthread, NULL);

  for (int i = 0; i < thread_number; ++i)
    answer += thread_args_array[i].answer;
    printf("%f\n", answer);
}
