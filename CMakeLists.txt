cmake_minimum_required(VERSION 3.13)

project(AbInitio)

include(CheckLibraryExists)
include(FindBLAS)
include(FindLAPACK)
include(CTest)

check_library_exists(m sin "" MATH)

find_library(LAPACKE_LIBRARIES lapacke)

find_path(LAPACKE_INC lapacke.h)

if(${MATH})
  set(MATH m)
endif()

add_library(abinit SHARED)
add_dependencies(abinit util)

add_subdirectory(src)
include_directories(src)

add_compile_options(-Wall -Wextra -Werror)

if(BUILD_TESTING)
  
  add_compile_options(--coverage -Og -g3 -ggdb)
  add_subdirectory(tests)

else()

  add_compile_options(-Ofast)

endif()
