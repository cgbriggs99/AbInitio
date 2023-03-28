cmake_minimum_required(VERSION 3.13)

project(AbInitio)

add_library(energy OBJECT energy.cpp)
set_target_properties(energy PROPERTIES POSITION_INDEPENDEND_CODE True)

target_sources(abinit PUBLIC $<TARGET_OBJECTS:energy>)

install(FILES energy.hpp TYPE INCLUDE)
 
