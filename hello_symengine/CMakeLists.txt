cmake_minimum_required(VERSION 3.10.0)

set(CMAKE_CXX_STANDARD 20)

set(CMAKE_TOOLCHAIN_FILE "D:\\Library\\vcpkg\\scripts\\buildsystems\\vcpkg.cmake")

set(CMAKE_INSTALL_PREFIX "D:/OwnLibrary")

project(hello_symengine)

add_executable(hello_world hello_world.cpp)

#  gmp
find_package(PkgConfig REQUIRED)
pkg_check_modules(gmp REQUIRED IMPORTED_TARGET gmp)
target_link_libraries(hello_world PkgConfig::gmp)

find_package(SymEngine CONFIG REQUIRED)
target_link_libraries(hello_world symengine)

# find_package(flint CONFIG REQUIRED)
# target_link_libraries(${PROJECT_NAME} flint)

# link_directories("D:\\Library\\vcpkg\\packages\\flint_x64-windows\\lib")
target_link_libraries(hello_world "D:\\Library\\vcpkg\\packages\\flint_x64-windows\\lib\\flint.lib")

enable_testing()
find_package(GTest CONFIG REQUIRED)
add_executable(unit_test unit_test.cpp multivariable_extream.cpp charles_symengine_common.cpp)
target_link_libraries(unit_test GTest::gtest GTest::gmock PkgConfig::gmp symengine)
target_link_libraries(unit_test "D:\\Library\\vcpkg\\packages\\flint_x64-windows\\lib\\flint.lib")
add_test(AllTestsInMain unit_test)


# set(EXPORT_NAME multivariable_extream)
add_library(${PROJECT_NAME} multivariable_extream.cpp charles_symengine_common.cpp)
target_link_libraries(${PROJECT_NAME} PkgConfig::gmp symengine)
target_link_libraries(${PROJECT_NAME} "D:\\Library\\vcpkg\\packages\\flint_x64-windows\\lib\\flint.lib")
target_link_libraries(${PROJECT_NAME} GTest::gtest GTest::gmock)

# 引入要用到的 CMake 模块
include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

# 基本安装及 Targets 文件的生成
install(TARGETS ${PROJECT_NAME}
  EXPORT ${PROJECT_NAME}-targets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Targets 文件的安装
install(EXPORT ${PROJECT_NAME}-targets
  FILE ${PROJECT_NAME}-targets.cmake
  NAMESPACE ${PROJECT_NAME}::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
)

# Config 文件的导出（与安装？）
configure_package_config_file(
  ${PROJECT_SOURCE_DIR}/cmake/config.cmake.in
  ${PROJECT_BINARY_DIR}/${PROJECT_NAME}-config.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
)

# ConfigVersion 文件的导出
write_basic_package_version_file(
  ${PROJECT_NAME}-config-version.cmake
  VERSION 0.1
  COMPATIBILITY AnyNewerVersion
)

# Config 和 ConfigVersion 文件的安装
install(FILES
  ${PROJECT_BINARY_DIR}/${PROJECT_NAME}-config.cmake
  ${PROJECT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
)
# install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/${PROJECT_NAME} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(FILES "${PROJECT_SOURCE_DIR}/multivariable_extream.h" DESTINATION include)
install(FILES "${PROJECT_SOURCE_DIR}/charles_symengine_common.h" DESTINATION include)


# there are alwarys some symbol problem in manually build library, some was symbol not found, some was stdxxx duplicate symbols in lib, so use symengine installed by vcpkg, and it works fine.
# include_directories("D:\\Library\\symengine-0.13.0")
# include_directories("D:\\Library\\symengine-0.13.0\\build")
# include_directories("D:\\Library\\symengine-0.13.0\\build\\symengine\\utilities\\teuchos")

# # link_directories("D:\\Library\\symengine-0.13.0\\build\\symengine\\Debug")
# # link_directories("D:\\Library\\symengine-0.13.0\\build\\symengine\\utilities\\catch\\Debug")
# # link_directories("D:\\Library\\symengine-0.13.0\\build\\benchmarks\\Debug")
# target_link_libraries(${PROJECT_NAME} "D:\\Library\\symengine-0.13.0\\build\\symengine\\Debug\\symengine.lib")
# target_link_libraries(${PROJECT_NAME} "D:\\Library\\symengine-0.13.0\\build\\symengine\\utilities\\teuchos\\Debug\\teuchos.lib")

