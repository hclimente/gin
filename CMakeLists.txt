cmake_minimum_required(VERSION 3.6)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_STANDARD 11)

if (CMAKE_BUILD_TYPE EQUAL "DEBUG")
	SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -g" )
endif (CMAKE_BUILD_TYPE EQUAL "DEBUG")

find_package(Boost)
include_directories(${Boost_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR}/include/
		${PROJECT_SOURCE_DIR}/lib ${PROJECT_SOURCE_DIR}/lib/eigen
		${PROJECT_SOURCE_DIR}/lib/maxflow ${PROJECT_SOURCE_DIR}/lib/Cephes)

# create gin library
if (UNIX)
	if (APPLE)
		SET( ginlib  "gin-macosx-x64" )
	else (APPLE)
		SET( ginlib  "gin-linux-x64" )
	endif (APPLE)
endif (UNIX)

file(GLOB_RECURSE SOURCES ${PROJECT_SOURCE_DIR}/src/*.cc)
add_library(${ginlib} STATIC ${SOURCES} $<TARGET_OBJECTS:maxflow> $<TARGET_OBJECTS:cephes>)

add_subdirectory (${PROJECT_SOURCE_DIR}/lib/Cephes)
add_subdirectory (${PROJECT_SOURCE_DIR}/lib/maxflow)
add_subdirectory (${PROJECT_SOURCE_DIR}/lib/googletest)


install(TARGETS ${ginlib}
		RUNTIME DESTINATION bin
		ARCHIVE DESTINATION lib)

# get tools and tests
add_subdirectory (${PROJECT_SOURCE_DIR}/test)
add_subdirectory (${PROJECT_SOURCE_DIR}/src/tools)
