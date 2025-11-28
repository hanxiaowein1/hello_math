# Welcome!
This is a project for learning some math, and building some functions based on other basic math library.
# Hello Symengine
This project is corrently main for solving quadric multivariable extream with some boundary condition. For example, given a quadric multivariable function f(x, y, z), and limit its domain in a tetrahedron area. The method is quite primitive, here is the detailed explanation:
```text
1. compute partial derivative for x, y, z, we get three linear equations
2. solve equations, we get a point, if it's in the limited area, return
3. if we cannot find extream values in the inner volume area, then we should find extream on each face(if cannot find extream value in face, then find extream value in the edge of face, if we still cannot find the extream value in edge, we can get two points from the edge, get the max/min one).
4. compare all values from all faces, then we get the extream value.
```
However, after watching Justin Solomon's Shape Analysis course, I realized this problem can be solved more easily, so probably this project will not be continuously developed.
## Build
Currently this project is only for self use, however, if you want to compile it locally, I can provide some ideas.
```text
1. download vcpkg, and install symengine/PkgConfig/Eigen/GTest/flint
2. modify the CMAKE_TOOLCHAIN_FILE to your own vcpkg toolchain location
3. modify the target_link_libaries(xxx xxx/flint.lib) to your flint location
4. if you want to install, modify the CMAKE_INSTALL_PREFIX variable to your wanted location
```
I'm pretty sure the above steps will be enough to build locally, however if you have some problems, feel free to contact me at [email](hanwein2@gmail.com).