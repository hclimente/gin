include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

add_executable(runCChi2Tests CChi2Tests.cc)

target_link_libraries(runCChi2Tests gtest gtest_main ${ginlib})