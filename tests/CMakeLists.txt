 
cmake_minimum_required(VERSION 3.13)

project(AbInitio)

add_executable(test_util test_util.cpp)
target_link_libraries(test_util abinit)
target_link_libraries(test_util gcov)

target_compile_options(test_util PUBLIC --coverage)

add_test(NAME valgrind_test_util COMMAND valgrind --leak-check=full --show-leak-kinds=all --errors-for-leak-kinds=all --track-origins=yes $<TARGET_FILE:test_util>)
add_test(NAME util_cov COMMAND gcov test_util.cpp -m)
