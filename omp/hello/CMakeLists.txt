cmake_minimum_required(VERSION 3.13)

add_executable(hello "hello.cpp")
target_link_libraries(hello OpenMP::OpenMP_CXX)
