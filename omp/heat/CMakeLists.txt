cmake_minimum_required(VERSION 3.13)

add_executable(heat "heat.cpp")
target_link_libraries(heat OpenMP::OpenMP_CXX)