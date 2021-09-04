#include <omp.h>
#include <utility>
#include <vector>
#include <ctime>
#include <random>
#include <iostream>
#include <cmath>
// Random Array Generation
using namespace std;

int HoarePartition(double* arr, int left_index, int right_index) {
  double pivot = arr[(right_index + left_index) / 2];
  int i = left_index - 1, j = right_index + 1;

  while (true) {
    do {
      i++;
    } while (arr[i] < pivot);

    do {
      j--;
    } while (arr[j] > pivot);

    if (i >= j)
      return j;

    swap(arr[i], arr[j]);
  }
}

// Quick Sort with Hoare separation of the array
void qHoareSort(double* arr, int left_index, int right_index) {
  if (left_index < right_index) {
    int pi = HoarePartition(arr, left_index, right_index);
    qHoareSort(arr, left_index, pi);
    qHoareSort(arr, pi + 1, right_index);
  }
}

void simple_fusion(double * arr, int n, int m)
{
  double * result = new double[n + m];
  double * a = arr;
  double * b = arr + n;

  int a_index = 0, b_index = 0, result_index = 0;
  while ((a_index < n) && (b_index < m))
  {
    if (a[a_index] < b[b_index])
    {
      result[result_index] = a[a_index];
      ++result_index;
      ++a_index;
    }
    else
    {
      result[result_index] = b[b_index];
      ++result_index;
      ++b_index;
    }
  }

  double * tmp;
  int tmp_index, tmp_max;
  if (a_index == n)
  {
    tmp = b;
    tmp_index = b_index;
    tmp_max = m;
  }
  else
  {
    tmp = a;
    tmp_index = a_index;
    tmp_max = n;
  }
  int delta = result_index - tmp_index;

  for(; tmp_index < tmp_max; ++tmp_index)
    result[tmp_index + delta] = tmp[tmp_index];

  for(int i = 0; i < n + m; ++i)
    arr[i] = result[i];
}

vector<int> get_ends(int thread_num, int n)
{
  int threads_value = omp_get_num_threads() - 1;
  int sub_length = ceil(n / float(threads_value));
  if (thread_num < threads_value)
  {
    vector<int> ends(2);
    ends[0] = (thread_num - 1) * sub_length;
    ends[1] = ends[0] + sub_length - 1;
    return ends;
  }
  else
  {
    vector<int> ends(2);
    ends[0] = (thread_num - 1) * sub_length;
    ends[1] = n - 1;
    return ends;
  }
}
 

void dump(vector<vector<int>> x)
{
  static int dump_num = 0;
  #pragma omp critical
  {
    ++dump_num;
    cout << "\n\n-----------" << '\n';
    cout << "num " << dump_num << '\n';
    cout << "thread " << omp_get_thread_num() << '\n'
         << "edge 1: " << x[0][0] << ' ' << x[0][1] << '\n'
         << "edge 2: " << x[1][0] << ' ' << x[1][1] << '\n';
    cout << "-----------" << "\n\n\n";
  }
}

void parallel_division_sort(double * arr, int n, int threads_value=4)
{
  omp_set_num_threads(threads_value + 1);
  int sub_length = ceil(n / float(threads_value));
  vector<vector<int>> splitted_ends(threads_value);
  #pragma omp parallel shared(splitted_ends, arr)
  {
    int thread_num = omp_get_thread_num();
    if (thread_num != 0)
    {
      vector<int> ends = get_ends(thread_num, n);
      splitted_ends[thread_num - 1] = ends;
      qHoareSort(arr, ends[0], ends[1]);
    }
  }
  #pragma omp barrier
  for (auto a : splitted_ends)
  {
    cout << endl;
    for (auto b : a)
      cout << ' ' << b;
  }
  cout << endl;
  int length = splitted_ends.size();

  while (length > 1)
  {
      omp_set_num_threads(length/2 + 1);
      #pragma omp parallel shared(splitted_ends, arr)
      {
        int thread_num = omp_get_thread_num();
        vector<vector<int>> sexy_variable_name(2); //!!!!!!!!!!!!!!!!!!!
        if (thread_num != 0)
        {
          int first_sub_section = (thread_num - 1) * 2 + (length % 2);
          int second_sub_section = first_sub_section + 1;
          /*cout << splitted_ends[first_sub_section][0] <<
           ' ' << splitted_ends[first_sub_section][1] <<
           ' ' << splitted_ends[second_sub_section][0] <<
           ' ' << splitted_ends[second_sub_section][1] << endl;*/
          sexy_variable_name[0] = splitted_ends[first_sub_section];
          sexy_variable_name[1] = splitted_ends[second_sub_section];
          dump(sexy_variable_name);
          int delta = sexy_variable_name[0][0],
              n = sexy_variable_name[0][1] - sexy_variable_name[0][0] + 1,
              m = sexy_variable_name[1][1] - sexy_variable_name[1][0] + 1;
          simple_fusion(arr + delta, n, m);
          sexy_variable_name[0][1] = sexy_variable_name[1][1];

          int insert_pos = thread_num - !(length % 2);
          splitted_ends[insert_pos] = sexy_variable_name[0];
        }
      }
      #pragma omp barrier
      length = (length + 1) / 2;
  }
}

int main()
{
  double a[] = {15, 10, 9, 8, 13, 11, 14, 7, 12, 6, 16, 5, 17, 4, 18, 3, 2, 1};
  parallel_division_sort(a, sizeof(a)/sizeof(double), 6);
  for(int i = 0; i < sizeof(a)/sizeof(double); ++i)
    cout << endl << a[i];
  return 0;
}
