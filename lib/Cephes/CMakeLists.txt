file(GLOB_RECURSE CEPHES_SOURCES cephes.h double/*.c )
list(REMOVE_ITEM CEPHES_SOURCES "${PROJECT_SOURCE_DIR}/lib/Cephes/double/setpmsvc.c"
		"${PROJECT_SOURCE_DIR}/lib/Cephes/double/lmdif.c")

add_library(cephes OBJECT ${CEPHES_SOURCES})