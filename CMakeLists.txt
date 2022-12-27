cmake_minimum_required(VERSION 2.8)
project(term_project)

set(CMAKE_BUILD_TYPE "Debug")

add_library(tool SHARED tool.cpp)

add_executable(LU LU.cpp)

target_link_libraries(LU tool)

add_executable(Jacobi Jacobi.cpp)

target_link_libraries(Jacobi tool)

add_executable(GS GS.cpp)

target_link_libraries(GS tool)

add_executable(data data.cpp)

target_link_libraries(data tool)

add_executable(sor sor.cpp)

target_link_libraries(sor tool)

add_executable(gmres gmres.cpp)

target_link_libraries(gmres tool)

add_executable(exp1 exp1.cpp)

target_link_libraries(exp1 tool)