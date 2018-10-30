cmake_minimum_required(VERSION 3.8.2)

set(CMAKE_C_COMPILER gcc)
project(heisenberg)

set(CMAKE_C_FLAGS "-O5")
set(CMAKE_C_FLAGS_DEBUG "-O0 -ggdb")
set(CMAKE_C_FLAGS_RELEASE "-O5 -ggdb")

# Sources directories
add_executable(heisenberg heisenberg.c)
target_link_libraries(heisenberg m)
