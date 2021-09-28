cmake_minimum_required(VERSION 3.13)

add_executable(matmul "matmul.cpp")
target_link_libraries(matmul OpenMP::OpenMP_CXX)
