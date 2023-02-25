 
cmake_minimum_required(VERSION 3.13)

project(AbInitio)

add_library(basis_sets OBJECT slater_orbital.cpp gaussian_orbital.cpp)
set_target_properties(basis_sets PROPERTIES POSITION_INDEPENDEND_CODE True)

target_sources(abinit PUBLIC $<TARGET_OBJECTS:basis_sets>)

add_subdirectory(sto-ng)

install(FILES basis_set.hpp TYPE INCLUDE)
