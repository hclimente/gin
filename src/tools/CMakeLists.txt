project(tools)

# add extra lib directories
add_executable(gin gin.cc)
target_link_libraries (gin boost_program_options ${ginlib})

add_executable(scones scones.cc)
target_link_libraries (scones ${ginlib})

add_executable(emmax emmax.cc)
target_link_libraries (emmax ${ginlib})

add_executable(fastlmm fastlmm.cc)
target_link_libraries (fastlmm ${ginlib})

add_executable(filter_snps filter_snps.cc)
target_link_libraries (filter_snps ${ginlib})

add_executable(linreg linreg.cc)
target_link_libraries (linreg ${ginlib})

add_executable(random_effect_meta_analysis random_effect_meta_analysis.cc)
target_link_libraries (random_effect_meta_analysis ${ginlib})
