cmake_minimum_required(VERSION 3.13)

add_executable(exp "exp.cpp")
target_link_libraries(exp OpenMP::OpenMP_CXX)
