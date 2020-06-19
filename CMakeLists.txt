cmake_minimum_required(VERSION 2.8)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project( PlaneDetector )

find_package(Eigen3 CONFIG REQUIRED)

add_executable( PlaneDetector 
planeFinder.cpp
PointCloud.cpp
PointCloud.h )

target_link_libraries( PlaneDetector PRIVATE Eigen3::Eigen)