cmake_minimum_required(VERSION 2.8)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project( PlaneDetector )

find_package(Eigen3 CONFIG REQUIRED)
find_package(tinyply CONFIG REQUIRED)

find_package(OpenMP)
if (OPENMP_FOUND)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

add_executable( PlaneDetector 
planeFinder.cpp
PointCloud.cpp
PointCloud.h   "PointCloud.cpp" "VectorPC.h")

target_link_libraries( PlaneDetector PRIVATE Eigen3::Eigen)
target_link_libraries( PlaneDetector PRIVATE tinyply)