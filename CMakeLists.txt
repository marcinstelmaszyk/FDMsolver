cmake_minimum_required(VERSION 3.16)
project(FDMsolver)

set(CMAKE_CXX_STANDARD 17)

add_executable(FDMsolver "src/main.cpp"
        "src/cGrid.cpp"
        "src/cGrid.h"
        "src/cPoint.cpp"
        "src/cPoint.h"
        "src/cSolution.cpp"
        "src/cSolution.h"
        "src/cSolver.cpp"
        "src/cSolver.h"
        "src/functions.h"
        )
