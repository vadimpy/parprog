cmake_minimum_required(VERSION 3.13)

project(
    parprog
    CXX
)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS})
set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS})

include_directories(${CMAKE_SOURCE_DIR}/include)

find_package(OpenMP REQUIRED)

set(
    subdirs
    example
    hello
    sum
    vecproc
    matmul
    exp
    simples
    qsort
    heat
    bugged
)

FOREACH(subdir IN LISTS subdirs)
    add_subdirectory(${subdir})
ENDFOREACH()
