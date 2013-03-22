file(GLOB Moirai_SOURCES "*.cpp")
include_directories(${Moirai_INCLUDE_DIRECTORIES})


add_library(Moirai SHARED ${Moirai_SOURCES})
TARGET_LINK_LIBRARIES(Moirai "/scratch/robinsonm/pkg/lib/libtrng4.so" ${VTK_LIBRARIES} ${MPI_CXX_LIBRARIES} ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})
