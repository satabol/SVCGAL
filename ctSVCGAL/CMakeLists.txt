# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12)
project( ctSVCGAL )

SET(CGAL_CMAKE_EXACT_NT_BACKEND BOOST_BACKEND CACHE STRING "")
SET(CGAL_DISABLE_GMP ON CACHE BOOL "")
SET(CMAKE_DISABLE_FIND_PACKAGE_GMP ON CACHE BOOL "")

set(Boost_USE_STATIC_LIBS   ON)
find_package(Boost REQUIRED COMPONENTS thread)
find_package(CGAL REQUIRED)
add_library(ctSVCGAL SHARED ctSVCGAL.cpp)
#target_link_directories(ctSVCGAL PRIVATE "E:/github.com/boost_1_86_0/stage/lib") # ex. for Windows
target_link_directories(ctSVCGAL PRIVATE "/opt/github.com/boost_1_86_0/stage/lib") # ex. for Linux
target_link_libraries(ctSVCGAL CGAL::CGAL Boost::thread)