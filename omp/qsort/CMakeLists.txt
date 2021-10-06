cmake_minimum_required(VERSION 3.13)

add_executable(qsort "qsort.cpp")
target_link_libraries(qsort OpenMP::OpenMP_CXX)
