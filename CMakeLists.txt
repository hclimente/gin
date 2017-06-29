cmake_minimum_required(VERSION 3.6)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_STANDARD 11)

if (CMAKE_BUILD_TYPE EQUAL "DEBUG")
	SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -g" )
endif (CMAKE_BUILD_TYPE EQUAL "DEBUG")

include_directories(${PROJECT_SOURCE_DIR}/include
		${PROJECT_SOURCE_DIR}/lib ${PROJECT_SOURCE_DIR}/lib/eigen
		${PROJECT_SOURCE_DIR}/lib/maxflow ${PROJECT_SOURCE_DIR}/lib/Cephes)

add_subdirectory (lib)
add_subdirectory (test)
add_subdirectory (src)