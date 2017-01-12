add_executable(test_gsdplugin test_gsdplugin.cc catch.cc)
target_link_libraries(test_gsdplugin gsdplugin)
add_test(gsdplugin test_gsdplugin)
