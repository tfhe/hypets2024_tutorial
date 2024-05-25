cmake_minimum_required(VERSION 3.0)
project(spqlios-hypets-samples)

add_subdirectory(lib/spqlios-arithmetic)
set(SPQLIOS_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/lib/spqlios-arithmetic)

add_subdirectory(01-simple-fft)
add_subdirectory(04-automorphisms)
add_subdirectory(05-onionpir)

