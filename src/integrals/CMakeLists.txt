 
cmake_minimum_required(VERSION 3.13)

project(AbInitio)

add_library(integrals OBJECT gaussian_ST.cpp
  gaussian_repulsion.cpp gaussian_attraction.cpp boys.cpp)
set_target_properties(integrals PROPERTIES POSITION_INDEPENDEND_CODE True)

target_sources(abinit PUBLIC $<TARGET_OBJECTS:integrals>)

install(FILES integrals.hpp integral_factory.template.cpp TYPE INCLUDE)
