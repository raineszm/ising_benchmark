cmake_minimum_required(VERSION 3.8)

project(ising_cpp)
find_package(Threads)
add_executable(ising src/ising.cpp)
target_link_libraries(ising Threads::Threads)
target_compile_features(ising PRIVATE cxx_std_14)
