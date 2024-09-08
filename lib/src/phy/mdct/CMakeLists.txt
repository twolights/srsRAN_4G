add_subdirectory(test)
file(GLOB SOURCES "*.c")
add_library(srsran_mdct OBJECT ${SOURCES})
