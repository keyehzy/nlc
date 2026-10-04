cmake_minimum_required(VERSION 3.24)
project(nlc_thermo LANGUAGES CXX)

if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

option(NLC_NATIVE "Tune for the build machine's CPU" ON)

find_package(Boost 1.81 REQUIRED CONFIG)
find_package(Threads REQUIRED)
if(APPLE)
  find_library(ACCELERATE Accelerate REQUIRED)
  set(NLC_LAPACK ${ACCELERATE})
else()
  find_package(LAPACK REQUIRED)
  set(NLC_LAPACK LAPACK::LAPACK)
endif()

add_library(nlc STATIC
  src/graph.cpp
  src/lattice.cpp
  src/expansion.cpp
  src/ed.cpp
  src/thermo.cpp
  src/resum.cpp
  src/nlc.cpp
)
target_include_directories(nlc PUBLIC src)
target_link_libraries(nlc PUBLIC Boost::headers Threads::Threads ${NLC_LAPACK})
target_compile_options(nlc PUBLIC -Wall -Wextra -Wno-sign-compare)

if(NLC_NATIVE)
  include(CheckCXXCompilerFlag)
  check_cxx_compiler_flag("-mcpu=native" NLC_HAS_MCPU_NATIVE)
  if(NLC_HAS_MCPU_NATIVE)
    target_compile_options(nlc PUBLIC -mcpu=native)
  else()
    check_cxx_compiler_flag("-march=native" NLC_HAS_MARCH_NATIVE)
    if(NLC_HAS_MARCH_NATIVE)
      target_compile_options(nlc PUBLIC -march=native)
    endif()
  endif()
endif()

add_executable(nlc_run src/main.cpp)
target_link_libraries(nlc_run PRIVATE nlc)
set_target_properties(nlc_run PROPERTIES OUTPUT_NAME nlc)

enable_testing()
add_executable(nlc_tests tests/unit_tests.cpp)
target_link_libraries(nlc_tests PRIVATE nlc)
add_test(NAME unit COMMAND nlc_tests)
