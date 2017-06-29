find_package(Boost)
include_directories(${Boost_INCLUDE_DIRS})

# create gin library
file(GLOB_RECURSE SOURCES ${PROJECT_SOURCE_DIR}/src/*.cc)
add_library(libgin STATIC ${SOURCES} $<TARGET_OBJECTS:maxflow> $<TARGET_OBJECTS:cephes>)

install(TARGETS libgin
		RUNTIME DESTINATION bin
		ARCHIVE DESTINATION lib)

# get tools
add_executable(gin tools/gin.cc)
target_link_libraries (gin boost_program_options libgin)

add_executable(scones tools/scones.cc)
target_link_libraries (scones libgin)

add_executable(emmax tools/emmax.cc)
target_link_libraries (emmax libgin)

add_executable(fastlmm tools/fastlmm.cc)
target_link_libraries (fastlmm libgin)

add_executable(filter_snps tools/filter_snps.cc)
target_link_libraries (filter_snps libgin)

add_executable(linreg tools/linreg.cc)
target_link_libraries (linreg libgin)

add_executable(random_effect_meta_analysis tools/random_effect_meta_analysis.cc)
target_link_libraries (random_effect_meta_analysis libgin)
