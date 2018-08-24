add_library(LCAdaptiveSampling STATIC LCAdaptiveGrid.h LCAdaptiveGridCell.h LCBasisFunction.h LCSampledFunction.h LCSample.h LCProtoConverter.h LCAdaptiveSamplingParams.h LCAdaptiveGrid.cpp LCAdaptiveGridCell.cpp LCBasisFunction.cpp LCSampledFunction.cpp LCSample.cpp LCProtoConverter.cpp)

include_directories(${CMAKE_SOURCE_DIR}/LCFunction)
include_directories(${CMAKE_SOURCE_DIR}/LCRealFunction)
include_directories(${CMAKE_SOURCE_DIR}/libs/protobuf/include)
include_directories(${CMAKE_SOURCE_DIR}/libs/eigen)
include_directories(${CMAKE_SOURCE_DIR}/LCUtils)