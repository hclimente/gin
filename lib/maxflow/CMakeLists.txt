project(MAXFLOW CXX)

add_library(maxflow OBJECT maxflow.cpp maxflow.h)
target_include_directories (maxflow INTERFACE ${MAXFLOW_SOURCE_DIR})