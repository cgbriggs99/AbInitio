 
cmake_minimum_required(VERSION 3.13)

project(AbInitio)

macro(test_recipe target source)
  if(${ARGC} GREATER 2)
    set(SOURCE_FILE ${ARGV2})
    configure_file(${source} ${ARGV2})
  else()
    set(SOURCE_FILE ${source})
  endif()
  add_executable(${target} ${SOURCE_FILE})
  target_link_libraries(${target} abinit gcov ${LAPACKE_LIB} ${MATH} ${LAPACK_LIBRARIES})
  target_compile_options(${target} PUBLIC --coverage)
  add_test(NAME valgrind_${target} COMMAND valgrind --leak-check=full
    --show-leak-kinds=all --errors-for-leak-kinds=all --track-origins=yes
    $<TARGET_FILE:${target}>)
  add_test(NAME gcov_${target} COMMAND gcov ${source} -m)
  unset(SOURCE_FILE)
endmacro()

test_recipe(test_util test_util.cpp)
test_recipe(test_input test_input.pre.cpp test_input.cpp)
test_recipe(test_geometric test_geometric.pre.cpp test_geometric.cpp)
