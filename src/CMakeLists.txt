project(scones2)

# add extra lib directories
link_directories(/bioinfo/users/hcliment/libs/boost_1_63/lib)

add_executable(scones2 ${PROJECT_SOURCE_DIR}/tools/scones2.cc ${SOURCES})

target_link_libraries (scones2 boost_program_options cephes maxflow CEasyGWAS)