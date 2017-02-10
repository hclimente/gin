cmake_minimum_required(VERSION 3.6)
project(scones2)

# set(CMAKE_CXX_STANDARD 11)
# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/cmake-build-debug/)

# add extra include directories
include_directories(/bioinfo/users/hcliment/libs/boost_1_63/include ${PROJECT_SOURCE_DIR}/src
                    ${PROJECT_SOURCE_DIR}/External ${PROJECT_SOURCE_DIR}/External/Eigen
                    ${PROJECT_SOURCE_DIR}/External/maxflow ${PROJECT_SOURCE_DIR}/External/Cephes)

add_subdirectory (${PROJECT_SOURCE_DIR}/External/Cephes)
add_subdirectory (${PROJECT_SOURCE_DIR}/External/maxflow)
add_subdirectory (${PROJECT_SOURCE_DIR}/tests)

file(GLOB_RECURSE SOURCES ${PROJECT_SOURCE_DIR}/src/CEasyGWAS/*.cc ${PROJECT_SOURCE_DIR}/src/*.h
     ${PROJECT_SOURCE_DIR}/External/Eigen/*.cc ${PROJECT_SOURCE_DIR}/External/Eigen/*.h)

# add extra lib directories
link_directories(/bioinfo/users/hcliment/libs/boost_1_63/lib)

add_executable(scones2 ${PROJECT_SOURCE_DIR}/src/tools/scones2.cc ${SOURCES})
target_link_libraries (scones2 boost_program_options)
target_link_libraries (scones2 cephes)
target_link_libraries (scones2 maxflow)