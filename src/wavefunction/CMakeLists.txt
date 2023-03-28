cmake_minimum_required(VERSION 3.13)

project(AbInitio)

add_library(wavefunction OBJECT wavefunction.cpp)
set_target_properties(wavefunction PROPERTIES POSITION_INDEPENDEND_CODE True)

target_sources(abinit PUBLIC $<TARGET_OBJECTS:wavefunction>)

install(FILES wavefunction.hpp TYPE INCLUDE)
