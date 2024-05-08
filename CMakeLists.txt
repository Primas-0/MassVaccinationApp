cmake_minimum_required(VERSION 3.27)
project(Project4)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(Project4
        vacdb.h
        vacdb.cpp
        mytest.cpp)