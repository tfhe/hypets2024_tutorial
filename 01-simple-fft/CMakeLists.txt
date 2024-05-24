add_executable(simple-fft simple-fft.cpp)
target_link_libraries(simple-fft libspqlios)
target_include_directories(simple-fft PRIVATE ${SPQLIOS_INCLUDE_DIR})
