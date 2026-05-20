cmake_minimum_required(VERSION 3.15)
project(gauss_solver LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(MSVC)
    add_compile_options(/W4 /O2)
else()
    add_compile_options(-Wall -Wextra -pedantic -O2)
endif()

include_directories(include)

add_library(solver_lib
    src/gauss_solver.cpp
    src/lu_solver.cpp
    src/matrix_utils.cpp
    src/timer.cpp
)

add_executable(gauss_solver_app src/main.cpp)
target_link_libraries(gauss_solver_app solver_lib)