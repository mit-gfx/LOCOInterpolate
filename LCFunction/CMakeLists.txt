include_directories(${CMAKE_SOURCE_DIR}/LCFunction)
include_directories(${CMAKE_SOURCE_DIR}/LCRealFunction)
include_directories(${CMAKE_SOURCE_DIR}/libs/protobuf/include)
include_directories(${CMAKE_SOURCE_DIR}/libs/eigen)
include_directories(${CMAKE_SOURCE_DIR}/LCUtils)

add_library(LCFunction STATIC LCFunctionDeriv.h LCFunctionValue.h LCFunctionValueMap.h LCOneDimLinearFunction.h LCFunction.h LCShapeInfoProtoConverter.h PrecomputedParametricShape.pb.h LCFunctionDeriv.cpp LCFunctionValue.cpp LCFunctionValueMap.cpp LCOneDimLinearFunction.cpp LCFunction.cpp LCShapeInfoProtoConverter.cpp PrecomputedParametricShape.pb.cc)