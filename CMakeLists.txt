CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

FILE(TO_CMAKE_PATH "$ENV{IMAGINEPP_ROOT}" d)

IF(NOT EXISTS "${d}")

    MESSAGE(FATAL_ERROR "Error: environment variable IMAGINEPP_ROOT=" "${d}")

ENDIF(NOT EXISTS "${d}")

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${d}/CMake")
FIND_PACKAGE(Imagine)



ImagineAddExecutable(Console test.cpp options.cpp options.h)
