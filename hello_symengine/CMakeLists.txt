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


# include_directories("D:\\Library\\symengine-0.13.0")
# include_directories("D:\\Library\\symengine-0.13.0\\build")
# include_directories("D:\\Library\\symengine-0.13.0\\build\\symengine\\utilities\\teuchos")

# # link_directories("D:\\Library\\symengine-0.13.0\\build\\symengine\\Debug")
# # link_directories("D:\\Library\\symengine-0.13.0\\build\\symengine\\utilities\\catch\\Debug")
# # link_directories("D:\\Library\\symengine-0.13.0\\build\\benchmarks\\Debug")
# target_link_libraries(${PROJECT_NAME} "D:\\Library\\symengine-0.13.0\\build\\symengine\\Debug\\symengine.lib")
# target_link_libraries(${PROJECT_NAME} "D:\\Library\\symengine-0.13.0\\build\\symengine\\utilities\\teuchos\\Debug\\teuchos.lib")

