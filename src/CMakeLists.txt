cmake_minimum_required(VERSION 3.13)

project(AbInitio)

add_subdirectory(util)
add_subdirectory(basis_sets)
#add_subdirectory(integrals)

include_directories(util)
