 
cmake_minimum_required(VERSION 3.13)

project(AbInitio)

add_executable(test_util test_util.cpp)
target_link_libraries(test_util abinit)
target_link_libraries(test_util gcov)

target_compile_options(test_util PUBLIC --coverage)

add_test(NAME test_util COMMAND test_util)
add_test(NAME util_cov COMMAND gcov test_util.cpp -m)
