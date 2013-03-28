file(GLOB Moirai_SOURCES "*.cpp")
file(GLOB Moirai_HEADERS "*.h")
include_directories(${Moirai_INCLUDE_DIRECTORIES})


add_library(Moirai SHARED ${Moirai_SOURCES})
TARGET_LINK_LIBRARIES(Moirai ${VTK_LIBRARIES} ${MPI_CXX_LIBRARIES} ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})

install (TARGETS Moirai DESTINATION lib)
install (FILES ${Moirai_HEADERS} DESTINATION include)