set(CMAKE_CXX_STANDARD 17)

set(SAMPLES_COMMON
        onionpir.h
        onionpir.cpp
        tiny_fhe.cpp
        tiny_fhe.h
)

# the regular test mimics what a real user would write: stick to the public api!
add_executable(samples_test samples_test.cpp ${SAMPLES_COMMON})
target_include_directories(samples_test PRIVATE ${SPQLIOS_INCLUDE_DIR})
target_link_libraries(samples_test libspqlios gtest gtest_main)
add_test(samples_test samples_test)

# the internal test may use the testlib
add_executable(samples_test_internal samples_test_internal.cpp ${SAMPLES_COMMON})
target_include_directories(samples_test_internal PRIVATE ${SPQLIOS_INCLUDE_DIR})
target_link_libraries(samples_test_internal libspqlios spqlios-testlib gtest gtest_main pthread m)
add_test(samples_test_internal samples_test_internal)

add_executable(samples_test_internal_bigN samples_test_internal.cpp ${SAMPLES_COMMON})
target_include_directories(samples_test_internal_bigN PRIVATE ${SPQLIOS_INCLUDE_DIR})
target_compile_definitions(samples_test_internal_bigN PRIVATE BIG_N_MODE)
target_link_libraries(samples_test_internal_bigN libspqlios spqlios-testlib gtest gtest_main pthread m)
# big N is not part of the default test suite of course!!
