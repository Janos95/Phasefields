find_package(GSL REQUIRED)
find_package(CGAL CONFIG REQUIRED COMPONENTS Core)

add_library(polygonizer STATIC
        polygonize_wireframe.cpp
        polygonize_wireframe.hpp
        isotropic_remeshing.cpp
        isotropic_remeshing.hpp
        )
#target_precompile_headers(polygonizer PRIVATE
#        exprtk.hpp
#        )
target_include_directories(polygonizer PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_link_libraries(polygonizer PRIVATE GSL::gsl CGAL::CGAL CGAL::CGAL_Core)
target_link_libraries(polygonizer PUBLIC Phasefield::visualization)
set_property(TARGET polygonizer PROPERTY POSITION_INDEPENDENT_CODE ON)
add_library(Phasefield::polygonizer ALIAS polygonizer)
