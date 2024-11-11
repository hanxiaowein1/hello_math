cmake_minimum_required(VERSION 3.20.0)

set(CMAKE_CXX_STANDARD 20)

set(CMAKE_TOOLCHAIN_FILE "D:\\Library\\vcpkg\\scripts\\buildsystems\\vcpkg.cmake")

project(hello_symengine)

add_executable(${PROJECT_NAME} hello_world.cpp)

#  gmp
find_package(PkgConfig REQUIRED)
pkg_check_modules(gmp REQUIRED IMPORTED_TARGET gmp)
target_link_libraries(${PROJECT_NAME} PkgConfig::gmp)

find_package(SymEngine CONFIG REQUIRED)
target_link_libraries(${PROJECT_NAME} symengine)

# find_package(flint CONFIG REQUIRED)
# target_link_libraries(${PROJECT_NAME} flint)

# link_directories("D:\\Library\\vcpkg\\packages\\flint_x64-windows\\lib")
target_link_libraries(${PROJECT_NAME} "D:\\Library\\vcpkg\\packages\\flint_x64-windows\\lib\\flint.lib")


enable_testing()
find_package(GTest CONFIG REQUIRED)
add_executable(unit_test unit_test.cpp multivariable_extream.cpp)
target_link_libraries(unit_test GTest::gtest GTest::gmock PkgConfig::gmp symengine)
target_link_libraries(unit_test "D:\\Library\\vcpkg\\packages\\flint_x64-windows\\lib\\flint.lib")
add_test(AllTestsInMain unit_test)

# there are alwarys some symbol problem in manually build library, some was symbol not found, some was stdxxx duplicate symbols in lib, so use symengine installed by vcpkg, and it works fine.
# include_directories("D:\\Library\\symengine-0.13.0")
# include_directories("D:\\Library\\symengine-0.13.0\\build")
# include_directories("D:\\Library\\symengine-0.13.0\\build\\symengine\\utilities\\teuchos")

# # link_directories("D:\\Library\\symengine-0.13.0\\build\\symengine\\Debug")
# # link_directories("D:\\Library\\symengine-0.13.0\\build\\symengine\\utilities\\catch\\Debug")
# # link_directories("D:\\Library\\symengine-0.13.0\\build\\benchmarks\\Debug")
# target_link_libraries(${PROJECT_NAME} "D:\\Library\\symengine-0.13.0\\build\\symengine\\Debug\\symengine.lib")
# target_link_libraries(${PROJECT_NAME} "D:\\Library\\symengine-0.13.0\\build\\symengine\\utilities\\teuchos\\Debug\\teuchos.lib")

