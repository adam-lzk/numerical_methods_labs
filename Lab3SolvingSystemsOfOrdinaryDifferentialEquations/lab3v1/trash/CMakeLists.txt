cmake_minimum_required(VERSION 3.23)
project(Eiler)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(Eiler
        Eiler.h
        eilerExplicit.cpp
        eilerImplicit.cpp
        Gauss.cpp
        Gauss.h
        Main.cpp
        Newton.cpp
        Newton.h)
