set(CMAKE_CXX_STANDARD 17)

# the regular test mimics what a real user would write: stick to the public api!
add_executable(automorphism automorphism.cpp)
target_include_directories(automorphism PRIVATE ${SPQLIOS_INCLUDE_DIR})
target_link_libraries(automorphism libspqlios gtest gtest_main)
