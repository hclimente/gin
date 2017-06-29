include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

add_executable(runGinTests ginTests.cc)

target_link_libraries(runGinTests gtest gtest_main libgin)