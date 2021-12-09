cmake_minimum_required(VERSION 3.13)

enable_language(C)
find_path(omp_h "omp.h")
include_directories(omp_h)
foreach(i RANGE 1 6)
    add_executable(bugged${i} "bugged${i}.c")
    target_link_libraries(bugged${i} OpenMP::OpenMP_CXX)
endforeach()

set(CMAKE_C_FLAGS "-Xpreprocessor -fopenmp")
