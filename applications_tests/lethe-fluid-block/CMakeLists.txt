set(TEST_TARGET lethe-fluid-block)

file(COPY 2d_channel.msh DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")
file(COPY cylinder_structured.msh DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")
file(COPY square.msh DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")

deal_ii_pickup_tests()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set_tests_properties(lethe-fluid-block/mms3d_gd.mpirun=2.debug PROPERTIES TIMEOUT 2000)
endif()

