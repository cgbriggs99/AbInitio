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

if(BUILD_TESTING)
  
  add_compile_options(--coverage -Og -g3 -ggdb -lgcov -pg)
  add_subdirectory(tests)

else()

  add_compile_options(-Ofast)

endif()

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
add_library(abinit SHARED)
add_dependencies(abinit util basis_sets input)

set_target_properties(abinit PROPERTIES POSITION_INDEPENDENT_CODE True)

add_subdirectory(src)
include_directories(src)

add_compile_options(-Wall -Wextra -Werror)

target_link_libraries(abinit ${LAPACKE_LIBRARIES} ${MATH})
target_include_directories(abinit PUBLIC ${LAPACKE_INC})

