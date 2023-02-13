 
cmake_minimum_required(VERSION 3.13)

project(AbInitio)

add_library(input OBJECT input.cpp)
set_target_properties(input PROPERTIES POSITION_INDEPENDENT_CODE True)

target_sources(abinit PUBLIC $<TARGET_OBJECTS:input>)
