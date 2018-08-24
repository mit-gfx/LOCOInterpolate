include_directories(${CMAKE_SOURCE_DIR}/LCFunction)
include_directories(${CMAKE_SOURCE_DIR}/LCRealFunction)
include_directories(${CMAKE_SOURCE_DIR}/libs/protobuf/include)
include_directories(${CMAKE_SOURCE_DIR}/libs/eigen)
include_directories(${CMAKE_SOURCE_DIR}/LCUtils)

add_library(LCRealFunction STATIC LCRealFunction.cpp LCRealFunction.h LCRealFunctionDeriv.h LCRealFunctionValue.h LCRealFunctionValueMap.h LCRealFunctionDeriv.cpp LCRealFunctionValue.cpp LCRealFunctionValueMap.cpp)