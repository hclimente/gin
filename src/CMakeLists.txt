project(apps)

# add extra lib directories
add_executable(gin ${PROJECT_SOURCE_DIR}/tools/gin.cc ${SOURCES})
target_link_libraries (gin boost_program_options cephes maxflow CEasyGWAS)

add_executable(scones ${PROJECT_SOURCE_DIR}/tools/scones.cc ${SOURCES})
target_link_libraries (scones cephes maxflow CEasyGWAS)

add_executable(emmax ${PROJECT_SOURCE_DIR}/tools/emmax.cc ${SOURCES})
target_link_libraries (emmax cephes maxflow CEasyGWAS)

add_executable(fastlmm ${PROJECT_SOURCE_DIR}/tools/fastlmm.cc ${SOURCES})
target_link_libraries (fastlmm cephes maxflow CEasyGWAS)

add_executable(filter_snps ${PROJECT_SOURCE_DIR}/tools/fastlmm.cc ${SOURCES})
target_link_libraries (filter_snps cephes maxflow CEasyGWAS)

add_executable(linreg ${PROJECT_SOURCE_DIR}/tools/fastlmm.cc ${SOURCES})
target_link_libraries (linreg cephes maxflow CEasyGWAS)

add_executable(random_effect_meta_analysis ${PROJECT_SOURCE_DIR}/tools/fastlmm.cc ${SOURCES})
target_link_libraries (random_effect_meta_analysis cephes maxflow CEasyGWAS)
