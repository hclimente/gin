include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

add_executable(runScones2Tests SCones2Test.cc)

target_link_libraries(runScones2Tests gtest gtest_main CEasyGWAS)