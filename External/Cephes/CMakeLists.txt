file(GLOB_RECURSE CEPHES_SOURCES
		${PROJECT_SOURCE_DIR}/External/Cephes/cephes.h ${PROJECT_SOURCE_DIR}/External/Cephes/double/*.c )
list(REMOVE_ITEM CEPHES_SOURCES "${PROJECT_SOURCE_DIR}/External/Cephes/double/setpmsvc.c")

add_library(cephes STATIC ${CEPHES_SOURCES})