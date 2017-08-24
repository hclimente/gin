project(TESTS)

include_directories(${GIN_SOURCE_DIR}/lib/eigen ${TESTS_SOURCE_DIR})

# unit tests
file (GLOB_RECURSE UNIT_TEST_SOURCES unit/*/*.cc)

add_executable (unit_tests ${UNIT_TEST_SOURCES})
target_link_libraries (unit_tests LINK_PUBLIC gtest gtest_main gmock gin)

# integration tests
file (GLOB_RECURSE INTEGRATION_TEST_SOURCES integration/*.cc)

add_executable (integration_tests ${INTEGRATION_TEST_SOURCES})
target_link_libraries (integration_tests LINK_PUBLIC gtest gtest_main gin)