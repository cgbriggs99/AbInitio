cmake_minimum_required(VERSION 3.13)

project(AbInitio)

add_library(util OBJECT atom.cpp molecule.cpp)
set_property(TARGET util PROPERTY POSITION_INDEPENDENT_CODE True)

target_sources(abinit PUBLIC $<TARGET_OBJECTS:util>)
