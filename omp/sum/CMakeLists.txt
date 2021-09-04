cmake_minimum_required(VERSION 3.13)

add_executable(sum "sum.cpp")
target_link_libraries(sum OpenMP::OpenMP_CXX)
