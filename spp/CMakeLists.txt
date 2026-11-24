cmake_minimum_required(VERSION 3.14)
project(spp VERSION 1.1.1 DESCRIPTION "Segmented primorial prefilter for Miller-Rabin primality testing" LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()

option(SPP_BUILD_TOOLS "Build command-line tools" ON)
option(SPP_BUILD_TESTS "Build unit tests" ON)

# ---- GMP
find_path(GMP_INCLUDE_DIR gmp.h)
find_library(GMP_LIBRARY NAMES gmp libgmp gmp-10 libgmp-10)
if(NOT GMP_INCLUDE_DIR OR NOT GMP_LIBRARY)
  message(FATAL_ERROR "GMP not found. Install it (apt install libgmp-dev / pacman -S mingw-w64-ucrt-x86_64-gmp / "
                      "vcpkg install gmp) or set GMP_INCLUDE_DIR and GMP_LIBRARY.")
endif()
add_library(GMP::gmp UNKNOWN IMPORTED)
set_target_properties(GMP::gmp PROPERTIES IMPORTED_LOCATION "${GMP_LIBRARY}"
                                          INTERFACE_INCLUDE_DIRECTORIES "${GMP_INCLUDE_DIR}")

# ---- library
add_library(spp src/spp.cpp)
add_library(spp::spp ALIAS spp)
target_include_directories(spp PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                                      $<INSTALL_INTERFACE:include>)
target_link_libraries(spp PUBLIC GMP::gmp)
if(WIN32)
  target_link_libraries(spp PRIVATE bcrypt)
endif()
if(MSVC)
  target_compile_options(spp PRIVATE /W4)
else()
  target_compile_options(spp PRIVATE -Wall -Wextra)
endif()

# ---- tools
if(SPP_BUILD_TOOLS)
  foreach(t spp_check spp_gen spp_bench)
    add_executable(${t} tools/${t}.cpp)
    target_link_libraries(${t} PRIVATE spp)
  endforeach()
  add_executable(spp_example examples/example.cpp)
  target_link_libraries(spp_example PRIVATE spp)
endif()

# ---- tests
if(SPP_BUILD_TESTS)
  enable_testing()
  add_executable(spp_tests tests/spp_tests.cpp)
  target_link_libraries(spp_tests PRIVATE spp)
  add_test(NAME spp_tests COMMAND spp_tests)
endif()

# ---- install
include(GNUInstallDirs)
install(TARGETS spp EXPORT sppTargets ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
if(SPP_BUILD_TOOLS)
  install(TARGETS spp_check spp_gen spp_bench RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
install(EXPORT sppTargets NAMESPACE spp:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/spp)
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/sppConfig.cmake
"include(CMakeFindDependencyMacro)
find_path(GMP_INCLUDE_DIR gmp.h)
find_library(GMP_LIBRARY NAMES gmp libgmp gmp-10 libgmp-10)
if(NOT TARGET GMP::gmp)
  add_library(GMP::gmp UNKNOWN IMPORTED)
  set_target_properties(GMP::gmp PROPERTIES IMPORTED_LOCATION \"\${GMP_LIBRARY}\" INTERFACE_INCLUDE_DIRECTORIES \"\${GMP_INCLUDE_DIR}\")
endif()
include(\"\${CMAKE_CURRENT_LIST_DIR}/sppTargets.cmake\")
")
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/sppConfig.cmake DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/spp)
