

enable_testing()
set_directory_properties(PROPERTIES CORRADE_USE_PEDANTIC_FLAGS ON)

#corrade_add_test(modica_mortola_gradient_test modica_mortola_gradient_test.cpp)
#target_link_libraries(modica_mortola_gradient_test PRIVATE phasefields::cost Magnum::Primitives MagnumIntegration::Eigen)

corrade_add_test(fw_test fw_test.cpp)
target_link_libraries(fw_test PRIVATE phasefields::cost)

corrade_add_test(dijkstra_test dijkstra_test.cpp)
target_link_libraries(dijkstra_test PRIVATE phasefields::cost)

#corrade_add_test(connectedness_constraint_test connectedness_constraint_test.cpp)
#target_link_libraries(connectedness_constraint_test PRIVATE phasefields::cost Magnum::Primitives MagnumIntegration::Eigen phasefields::io)
