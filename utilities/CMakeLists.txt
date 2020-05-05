find_package(boost_mp11 CONFIG REQUIRED)
add_library(utilities INTERFACE)
target_sources(utilities INTERFACE
        robin_map.h
        robin_growth_policy.h
        robin_hash.h
        unique_function.h
        detail/pointer_union.h
        detail/pointer_int_pair.h
        detail/pointer_like_type_traits.h
        variant.hpp
        shared_ressource.hpp
        function_ref.hpp
        )

target_link_libraries(utilities INTERFACE Boost::mp11)
target_include_directories(utilities INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
add_library(phasefields::utilities ALIAS utilities)

#add_executable(map_test test.cpp)