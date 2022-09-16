cmake_minimum_required(VERSION 3.4)
project(TMCs)

# C++
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -Wno-return-type -fPIC ")
set(CMAKE_CXX_COMPILER g++)

include_directories(
    ./TMCsrc
    ./cubature/
)
file(GLOB_RECURSE LIB_CUBATURE_SRC "cubature/*.c")
file(GLOB_RECURSE LIB_TMC_SRC "TMCsrc/*.cpp")

#compile cubature libs
add_library(cubature STATIC ${LIB_CUBATURE_SRC})
set_target_properties(cubature PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ./lib)

#compile TMC lib
add_library(TMC STATIC ${LIB_TMC_SRC})
target_link_libraries(TMC cubature)
set_target_properties(TMC PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ./lib)

#compile MainTMC.cpp
add_executable(TMCdummy example/MainTMC.cpp)
target_link_libraries(TMCdummy TMC)
set_target_properties(TMCdummy PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/example")