find_package(GSL REQUIRED)
find_package(CGAL CONFIG REQUIRED COMPONENTS Core)

add_library(polygonizer STATIC
        polygonize_expression.hpp
        polygonize_expression.cpp
        )
#target_precompile_headers(polygonizer PRIVATE
#        exprtk.hpp
#        )
target_include_directories(polygonizer PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_link_libraries(polygonizer PRIVATE GSL::gsl CGAL::CGAL CGAL::CGAL_Core)
set_property(TARGET polygonizer PROPERTY POSITION_INDEPENDENT_CODE ON)
add_library(phasefields::polygonizer ALIAS polygonizer)
