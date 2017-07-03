include_directories(${GIN_SOURCE_DIR}/lib/eigen)

add_executable (runCChi2Tests CChi2Tests.cc)
target_link_libraries (runCChi2Tests gtest gtest_main gin)