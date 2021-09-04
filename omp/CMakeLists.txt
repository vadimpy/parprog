cmake_minimum_required(VERSION 3.13)

project(
    parprog
    CXX
)

set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS})

include_directories(${CMAKE_SOURCE_DIR}/include)

find_package(OpenMP REQUIRED)

set(
    subdirs
    example
    hello
    sum
)

FOREACH(subdir IN LISTS subdirs)
    add_subdirectory(${subdir})
ENDFOREACH()
