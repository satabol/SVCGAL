## Description

Binary background for nodes of Sverchok: https://github.com/nortikin/sverchok

- SvStraightSkeleton2DOffset
- SvStraightSkeleton2DExtrude

## Dependency

- CGAL (https://www.cgal.org/)
- Boost (https://www.boost.org/)

## Build

1. Download and build CGAL

git clone https://github.com/CGAL/cgal.git
git branch -m 6.0.x-branch

2. Download and build Boost

Build boost

Windows:

bootstrap.bat
b2 link=static

Linux:

./bootstrap.sh
./b2 cxxflags=-fPIC cflags=-fPIC install --prefix=/opt/github.com/boost_1_86_0/stage

3. Build ctypes_SVCGAL
3.1. Build Windows
3.2. Build Linux

4. Publish pypi