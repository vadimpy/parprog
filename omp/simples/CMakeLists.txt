cmake_minimum_required(VERSION 3.13)

add_executable(simples "simples.cpp")
target_link_libraries(simples OpenMP::OpenMP_CXX)
