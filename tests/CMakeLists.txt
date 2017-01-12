add_executable(test_gsdplugin test_gsdplugin.cc catch.cc)
target_link_libraries(test_gsdplugin gsdplugin dl)
add_test(gsdplugin test_gsdplugin)

add_custom_command(TARGET test_gsdplugin POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy
                           ${CMAKE_CURRENT_SOURCE_DIR}/test.gsd
                           ${CMAKE_CURRENT_BINARY_DIR}/test.gsd)
