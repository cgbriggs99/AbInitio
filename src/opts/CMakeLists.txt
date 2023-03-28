cmake_minimum_required(VERSION 3.13)

project(AbInitio)

add_library(opts OBJECT options.cpp default_options.cpp)
set_target_properties(opts PROPERTIES POSITION_INDEPENDENT_CODE True)

target_sources(abinit PUBLIC $<TARGET_OBJECTS:opts>)

install(FILES options.hpp default_options.hpp TYPE INCLUDE)
