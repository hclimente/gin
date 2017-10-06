########################################################################
#
# Project-wide settings

# Name of the project.
#
# CMake files in this project can refer to the root source directory
# as ${LIBGIN_SOURCE_DIR} and to the root binary directory as
# ${LIBGIN_BINARY_DIR}.
project (LIBGIN)

########################################################################
#
# GIN library libgin

# Different names for release and debug
if (CMAKE_BUILD_TYPE MATCHES Debug)
	set ( libgin  "gin-devel" )
else (CMAKE_BUILD_TYPE MATCHES Debug)
	set ( libgin  "gin" )
endif (CMAKE_BUILD_TYPE MATCHES Debug)

# Where gin's .h files can be found.
include_directories (${GIN_SOURCE_DIR}/include ${GIN_SOURCE_DIR}/lib/eigen ${GIN_SOURCE_DIR}/lib)

# *.cc files for gin, removing the CLI tools
file (GLOB_RECURSE SOURCES "${GIN_SOURCE_DIR}/src/*.cc")
list (REMOVE_ITEM SOURCES
		${GIN_SOURCE_DIR}/src/tools/emmax.cc
		${GIN_SOURCE_DIR}/src/tools/filter_snps.cc
		${GIN_SOURCE_DIR}/src/tools/random_effect_meta_analysis.cc
		${GIN_SOURCE_DIR}/src/tools/shake.cc
		${GIN_SOURCE_DIR}/src/tools/fastlmm.cc
		${GIN_SOURCE_DIR}/src/tools/linreg.cc
		${GIN_SOURCE_DIR}/src/tools/scones.cc )

# locate boost
find_package (Boost)
include_directories (${Boost_INCLUDE_DIRS})

# create library
add_library (${libgin} SHARED ${SOURCES} $<TARGET_OBJECTS:maxflow> $<TARGET_OBJECTS:cephes>)
target_link_libraries (${libgin} boost_program_options)
target_include_directories (${libgin} INTERFACE ${GIN_SOURCE_DIR}/include ${GIN_SOURCE_DIR}/lib)

# installation configuration
install (TARGETS ${libgin}
		 LIBRARY DESTINATION lib INCLUDES DESTINATION include)
install (DIRECTORY ${GIN_SOURCE_DIR}/include/ DESTINATION include
		 FILES_MATCHING PATTERN "*.h")
install (DIRECTORY ${GIN_SOURCE_DIR}/lib/Cephes DESTINATION include
		 FILES_MATCHING PATTERN "*.h")
install (DIRECTORY ${GIN_SOURCE_DIR}/lib/eigen/Eigen DESTINATION include)

# define alias for linking purposes (in TESTS)
add_library (gin ALIAS ${libgin})

########################################################################
#
# Compile the CLI tools

# compile tools
add_executable (evo tools/evo.cc)
target_link_libraries (evo LINK_PUBLIC boost_program_options gin)

add_executable (shake tools/shake.cc)
target_link_libraries (shake LINK_PUBLIC boost_program_options gin)

add_executable (scones tools/scones.cc)
target_link_libraries (scones LINK_PUBLIC gin)

# installation configuration
install (TARGETS shake scones
		RUNTIME DESTINATION bin)