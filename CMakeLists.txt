cmake_minimum_required(VERSION 3.13)

project(AbInitio)

include(CheckLibraryExists)
include(FindBLAS)
include(FindLAPACK)
include(CTest)

check_library_exists(m sin "" MATH)

message("-- Finding LAPACKE")

option(LAPACKE_LIB_DIR "directory to find liblapacke" /usr/lib/)
option(LAPACKE_INC_DIR "directory to find lapacke.h" /usr/include)

find_library(LAPACKE_LIB lapacke HINTS ${LAPACKE_LIB_DIR} REQUIRED)
message("-- Found lapacke library at " ${LAPACKE_LIB})
find_file(LAPACKE_INC lapacke.h HINTS ${LAPACKE_INC_DIR} REQUIRED)
message("-- Found lapacke.h at " ${LAPACKE_INC})

if(${MATH})
  set(MATH m)
endif()

if(BUILD_TESTING)
  message("-- Looking for valgrind in " $ENV{CONDA_PREFIX} " and " $ENV{CONDA})
  find_program(VALGRIND valgrind HINTS $ENV{CONDA_PREFIX} $ENV{CONDA})
  message("-- valgrind: " ${VALGRIND})
  add_compile_options(--coverage -Og -g3 -ggdb -lgcov -pg)
  add_subdirectory(tests)

else()

  add_compile_options(-Ofast)

endif()

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
add_library(abinit SHARED)
add_dependencies(abinit util basis_sets input)

if(BUILD_TESTING)
  target_link_libraries(abinit gcov)
endif()


set_target_properties(abinit PROPERTIES POSITION_INDEPENDENT_CODE True)

add_subdirectory(src)
include_directories(src)

add_compile_options(-Wall -Wextra -Werror)

target_link_libraries(abinit ${LAPACKE_LIB} ${MATH} ${LAPACK_LIBRARIES})

