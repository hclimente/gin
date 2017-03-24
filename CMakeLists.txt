cmake_minimum_required(VERSION 3.6)

set(CMAKE_CXX_STANDARD 11)
# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/cmake-build-debug/)
file(GLOB_RECURSE SOURCES ${PROJECT_SOURCE_DIR}/src/CEasyGWAS/*.cc ${PROJECT_SOURCE_DIR}/src/*.h
        ${PROJECT_SOURCE_DIR}/src/CEasyGWAS/*/*.cc ${PROJECT_SOURCE_DIR}/src/*/*.h
        ${PROJECT_SOURCE_DIR}/src/External/Eigen/*.cc ${PROJECT_SOURCE_DIR}/src/External/Eigen/*.h)
find_package(Boost)
include_directories(${Boost_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR}/src/
        ${PROJECT_SOURCE_DIR}/External ${PROJECT_SOURCE_DIR}/External/Eigen
        ${PROJECT_SOURCE_DIR}/External/maxflow ${PROJECT_SOURCE_DIR}/External/Cephes)

add_library(CEasyGWAS SHARED ${SOURCES})
target_link_libraries (CEasyGWAS cephes maxflow)

add_subdirectory (${PROJECT_SOURCE_DIR}/rscones2)
add_subdirectory (${PROJECT_SOURCE_DIR}/External/Cephes)
add_subdirectory (${PROJECT_SOURCE_DIR}/External/maxflow)
add_subdirectory (${PROJECT_SOURCE_DIR}/tests)
add_subdirectory (${PROJECT_SOURCE_DIR}/src)

install(TARGETS CEasyGWAS
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib)