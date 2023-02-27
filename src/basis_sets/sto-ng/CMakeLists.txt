cmake_minimum_required(VERSION 3.13)

project(AbInitio)

add_library(stong OBJECT stong.cpp)
set_target_properties(stong PROPERTIES POSITION_INDEPENDENT_CODE True)
target_sources(abinit PUBLIC $<TARGET_OBJECTS:stong>)

install(FILES stong.hpp TYPE INCLUDE)
