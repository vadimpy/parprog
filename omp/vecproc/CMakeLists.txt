cmake_minimum_required(VERSION 3.13)

add_executable(vecproc "vecproc.cpp")
target_link_libraries(vecproc OpenMP::OpenMP_CXX)
