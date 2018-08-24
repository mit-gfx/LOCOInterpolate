cmake_minimum_required (VERSION 3.0)
add_library(LCUtils STATIC LCError.h LCError.cpp LCTimer.h LCMathHelper.h LCMathHelper.cpp)

include_directories(${CMAKE_SOURCE_DIR}/libs/eigen)
add_definitions(-D_MBCS)
