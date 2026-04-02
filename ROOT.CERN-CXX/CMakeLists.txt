cmake_minimum_required(VERSION 3.16)
project(root_benchmark)

# Point this to your ROOT build directory
set(ROOT_DIR "c:/Users/proek/Documents/GitHub/root/out/build/x64-Debug")

# Find ROOT using its CMake config
find_package(ROOT REQUIRED COMPONENTS MathCore)

add_executable(root_benchmark root_benchmark.cpp)
target_link_libraries(root_benchmark ROOT::MathCore)