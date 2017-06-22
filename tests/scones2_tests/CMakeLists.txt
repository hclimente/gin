include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

add_executable(runGinTests ginTests.cc)
add_executable(runCChi2Tests CChi2Tests.cc)

target_link_libraries(runGinTests gtest gtest_main CEasyGWAS)
target_link_libraries(runCChi2Tests gtest gtest_main CEasyGWAS)