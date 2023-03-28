 
cmake_minimum_required(VERSION 3.13)

project(AbInitio)

configure_file(rhf.pre.cpp rhf.cpp)

add_library(scf OBJECT scf.cpp rhf.cpp)
set_target_properties(scf PROPERTIES POSITION_INDEPENDEND_CODE True)

target_sources(abinit PUBLIC $<TARGET_OBJECTS:scf>)

install(FILES scf.hpp rhf.hpp TYPE INCLUDE)
