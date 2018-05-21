cmake_minimum_required(VERSION 3.9)
project(oop_nodes_genetic)

set(CMAKE_CXX_STANDARD 11)

add_executable(oop_nodes_genetic
        Functions.h
        Functions.cpp
        Node.h
        Node.cpp
        oop_nodes_genetic.cpp ModelParameters.hpp)
