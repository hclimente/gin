include_directories(${GIN_SOURCE_DIR}/lib/eigen)

add_executable (runGinTests ginTests.cc)
target_link_libraries (runGinTests gtest gtest_main gin)