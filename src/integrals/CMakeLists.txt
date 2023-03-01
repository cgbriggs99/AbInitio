 
cmake_minimum_required(VERSION 3.13)

project(AbInitio)

add_library(integrals OBJECT gaussian.cpp numeric.cpp)
set_target_properties(integrals PROPERTIES POSITION_INDEPENDEND_CODE True)

target_sources(abinit PUBLIC $<TARGET_OBJECTS:integrals>)

install(FILES integrals.hpp TYPE INCLUDE)
