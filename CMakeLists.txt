cmake_minimum_required(VERSION 3.23)

project(numerical-solver LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(finite_difference)
add_subdirectory(laplace)
