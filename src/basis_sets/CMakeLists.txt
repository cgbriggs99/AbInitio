 
cmake_minimum_required(VERSION 3.13)

project(AbInitio)

add_library(basis_sets OBJECT slater_orbital.cpp compute.cpp electron_configs.cpp gaussian_orbital.cpp)
set_property(TARGET basis_sets PROPERTY POSITION_INDEPENDEND_CODE True)

target_sources(abinit PUBLIC $<TARGET_OBJECTS:basis_sets>)
