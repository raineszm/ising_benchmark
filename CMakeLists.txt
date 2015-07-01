cmake_minimum_required(VERSION 2.6)
project(ising)

find_package(GSL REQUIRED)
set(ISING_SRCS
    src/ising.cpp
    src/Lattice.cpp)

set(INCLUDE_DIRS src/ ${GSL_INCLUDE_DIRS})
list(APPEND LIBS ${GSL_LIBRARIES})

add_executable(ising ${ISING_SRCS})
target_link_libraries(ising ${LIBS})
set_property(TARGET ising PROPERTY CXX_STANDARD 11)
