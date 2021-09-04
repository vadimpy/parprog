cmake_minimum_required(VERSION 3.13)

add_executable(example "example.cpp")
target_link_libraries(example OpenMP::OpenMP_CXX)
