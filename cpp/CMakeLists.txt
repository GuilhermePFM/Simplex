cmake_minimum_required(VERSION 3.16)
project(simplex CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Eigen3 REQUIRED)

add_compile_options(-Wall -Wextra)

add_library(simplex_lib
    src/pivot.cpp
    src/preprocess.cpp
    src/logger.cpp
    src/phase1.cpp
    src/phase2.cpp
    src/simplex.cpp)
target_include_directories(simplex_lib PUBLIC include)
target_link_libraries(simplex_lib Eigen3::Eigen)

add_executable(run_tests
    test/main.cpp
    test/test_phase2.cpp
    test/test_phase1.cpp
    test/test_anticycling.cpp)
target_link_libraries(run_tests simplex_lib)

add_executable(bench benchmark/bench.cpp)
target_link_libraries(bench simplex_lib)
