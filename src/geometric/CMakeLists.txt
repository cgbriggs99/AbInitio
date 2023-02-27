 
cmake_minimum_required(VERSION 3.13)

project(AbInitio)

configure_file(rotor.pre.cpp rotor.cpp)
add_library(geometric OBJECT rotor.cpp)
set_target_properties(geometric PROPERTIES POSITION_INDEPENDENT_CODE True)

target_sources(abinit PUBLIC $<TARGET_OBJECTS:geometric>)

install(FILES rotor.hpp TYPE INCLUDE)
