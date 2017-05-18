project(gin)

# add extra lib directories
add_executable(gin ${PROJECT_SOURCE_DIR}/tools/gin.cc ${SOURCES})

target_link_libraries (gin boost_program_options cephes maxflow CEasyGWAS)
