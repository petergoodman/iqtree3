# Unit tests of the pure module model/fixedeqchart.cpp (decision 012; docs/agent/CODE_PLAN.md 3.2).
# A standalone project: it changes no IQ-TREE build file and compiles only the module and the tests.
#
#   cmake -S unittest -B BUILD -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_BUILD_TYPE=Release
#   cmake --build BUILD -j 2 && BUILD/fixedeqchart_tests
#
# test_scripts/fixedeq/unittests/unittests.py runs this in the Release and sanitizer variants.
cmake_minimum_required(VERSION 3.14)
if (POLICY CMP0135)
    cmake_policy(SET CMP0135 NEW)
endif()
project(fixedeqchart_tests CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# IQ-TREE's clang Release flags (CMakeLists.txt:420), without CMake's default -DNDEBUG, which
# IQ-TREE never defines, so Eigen's assertions stay on as they are in IQ-TREE; and the SSE flags
# IQ-TREE gives the model target on x86-64 (CMakeLists.txt:711-718 and 984)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ffunction-sections -fdata-sections")
if (CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64")
    add_compile_options(-D__SSE3 -msse3)
endif()

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# googletest at the commit cmaple pins (cmaple/CMakeLists.txt:287-293), version 1.12.1
set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
include(FetchContent)
FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
    URL_HASH SHA256=edd885a1ab32b6999515a880f669efadb80b3f880215f315985fa3f6eca7c4d3
)
FetchContent_MakeAvailable(googletest)

add_executable(fixedeqchart_tests test_fixedeqchart.cpp ${CMAKE_CURRENT_SOURCE_DIR}/../model/fixedeqchart.cpp)
target_include_directories(fixedeqchart_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../model)
target_compile_options(fixedeqchart_tests PRIVATE -Wall -Wextra)
target_compile_definitions(fixedeqchart_tests PRIVATE FIXTURE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/fixtures")
target_link_libraries(fixedeqchart_tests PRIVATE Eigen3::Eigen GTest::gtest_main)

enable_testing()
add_test(NAME fixedeqchart_tests COMMAND fixedeqchart_tests)
