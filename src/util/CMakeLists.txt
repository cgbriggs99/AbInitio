cmake_minimum_required(VERSION 3.13)

project(AbInitio)

add_library(util OBJECT atom.cpp molecule.cpp electron_configs.cpp
  polynomial.cpp)
set_target_properties(util PROPERTIES POSITION_INDEPENDENT_CODE True)

target_sources(abinit PUBLIC $<TARGET_OBJECTS:util>)
