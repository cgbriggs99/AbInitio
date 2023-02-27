 
cmake_minimum_required(VERSION 3.13)

project(AbInitio)

add_library(input OBJECT input.cpp psi4_bs_reader.cpp)
set_target_properties(input PROPERTIES POSITION_INDEPENDENT_CODE True)

target_sources(abinit PUBLIC $<TARGET_OBJECTS:input>)

install(FILES input.hpp psi4_bs_reader.hpp TYPE INCLUDE)
