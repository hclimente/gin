project (CEPHES C)

file (GLOB_RECURSE CEPHES_SOURCES double/*.c )
list (REMOVE_ITEM CEPHES_SOURCES
		${CEPHES_SOURCE_DIR}/double/setpmsvc.c
		${CEPHES_SOURCE_DIR}/double/lmdif.c)

add_library (cephes OBJECT ${CEPHES_SOURCES} cephes.h)
target_include_directories (cephes PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})