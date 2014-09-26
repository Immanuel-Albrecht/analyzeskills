cmake_minimum_required (VERSION 2.6)

project( ganter )

SET(CMAKE_C_FLAGS_DEBUG "-O2 -Wall -std=c9x -g -D_REENTRANT")
SET(CMAKE_C_FLAGS_RELEASE "-O3 -Wall -std=c9x -g -D_REENTRANT")

add_subdirectory(src)

