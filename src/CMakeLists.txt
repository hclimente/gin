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
	set ( libgin  "gin-db" )
else (CMAKE_BUILD_TYPE MATCHES Debug)
	set ( libgin  "gin" )
endif (CMAKE_BUILD_TYPE MATCHES Debug)

# Where gin's .h files can be found.
include_directories (${GIN_SOURCE_DIR}/include ${GIN_SOURCE_DIR}/lib/eigen ${GIN_SOURCE_DIR}/lib)

# *.cc files for gin, removing the CLI tools
file (GLOB_RECURSE SOURCES "${GIN_SOURCE_DIR}/src/*.cc")
list (FILTER SOURCES EXCLUDE REGEX "${GIN_SOURCE_DIR}/src/tools/")

# create library
add_library (${libgin} STATIC ${SOURCES} $<TARGET_OBJECTS:maxflow> $<TARGET_OBJECTS:cephes>)
target_include_directories (${libgin} INTERFACE ${GIN_SOURCE_DIR}/include ${GIN_SOURCE_DIR}/lib)

# installation configuration
install (TARGETS ${libgin}
		 ARCHIVE DESTINATION lib INCLUDES DESTINATION include)
install (DIRECTORY ${GIN_SOURCE_DIR}/include/ DESTINATION include
		 FILES_MATCHING PATTERN "*.h")
install (DIRECTORY ${GIN_SOURCE_DIR}/lib/Cephes DESTINATION include
		 FILES_MATCHING PATTERN "*.h")

# define alias for linking purposes (in TESTS)
add_library (gin ALIAS ${libgin})

########################################################################
#
# Compile the CLI tools

# fetch boost
find_package (Boost)
include_directories (${Boost_INCLUDE_DIRS})

# compile tools
add_executable (shake tools/shake.cc)
target_link_libraries (shake LINK_PUBLIC boost_program_options gin)

add_executable (scones tools/scones.cc)
target_link_libraries (scones LINK_PUBLIC gin)

add_executable (emmax tools/emmax.cc)
target_link_libraries (emmax LINK_PUBLIC gin)

add_executable (fastlmm tools/fastlmm.cc)
target_link_libraries (fastlmm LINK_PUBLIC gin)

add_executable (filter_snps tools/filter_snps.cc)
target_link_libraries (filter_snps LINK_PUBLIC gin)

add_executable (linreg tools/linreg.cc)
target_link_libraries (linreg LINK_PUBLIC gin)

add_executable (random_effect_meta_analysis tools/random_effect_meta_analysis.cc)
target_link_libraries (random_effect_meta_analysis LINK_PUBLIC gin)

# installation configuration
install (TARGETS shake scones
		RUNTIME DESTINATION bin)