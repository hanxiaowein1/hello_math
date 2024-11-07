#include <iostream>
#include <vector>
#include <symengine/expression.h>
#include "symengine/solve.h"
#include <symengine/basic.h>
#include <symengine/add.h>
#include <symengine/matrix.h>
// using SymEngine::Expression;

void basic_test()
{
    SymEngine::Expression x("x");
    auto ex = pow(x+sqrt(SymEngine::Expression(2)), 6);
    SymEngine::expand(ex);
    std::cout << ex << std::endl;
}

void derivative_test()
{
    // total 3 points
    std::vector<std::vector<double>> points{
        {1, 0, 0},
        {0, 2, 0},
        {0, 0, 3}
    };
    std::vector<std::vector<double>> normals{
        {1, 2, 3},
        {4, 5, 7},
        {9, 10, 3}
    };

    SymEngine::Expression x("x");
    SymEngine::Expression y("y");
    SymEngine::Expression z("z");
    SymEngine::Expression ex;
    // for(const auto& point: points)
    for(auto point_iterator = points.begin(); point_iterator != points.end(); point_iterator++)
    {
        auto point = *point_iterator;
        auto index = std::distance(points.begin(), point_iterator);
        auto normal = normals[index];
        ex = ex + pow(normal[0] * (x - point[0]) + normal[1] * (y - point[1]) + normal[2] * (z - point[2]), 2);
        SymEngine::expand(ex);
    }
    SymEngine::expand(ex);
    std::cout << ex << std::endl;
    SymEngine::Expression derivative_x = ex.diff(x);
    // std::cout << "partial derivative of x: " << derivative_x << std::endl;
    SymEngine::Expression derivative_y = ex.diff(y);
    SymEngine::Expression derivative_z = ex.diff(z);
}

void solve_poly_test()
{     
    // total 3 points
    std::vector<std::vector<double>> points{
        {1, 0, 0},
        {0, 2, 0},
        {0, 0, 3}
    };
    std::vector<std::vector<double>> normals{
        {1, 2, 3},
        {4, 5, 7},
        {9, 10, 3}
    };

    auto x = SymEngine::symbol("x");
    auto y = SymEngine::symbol("y");
    auto z = SymEngine::symbol("z");
    SymEngine::Expression x_(x);
    SymEngine::Expression y_(y);
    SymEngine::Expression z_(z);

    SymEngine::Expression ex;
    for (auto point_iterator = points.begin(); point_iterator != points.end(); point_iterator++)
    {
        auto point = *point_iterator;
        auto index = std::distance(points.begin(), point_iterator);
        auto normal = normals[index];
        ex = ex + pow(normal[0] * (x_ - point[0]) + normal[1] * (y_ - point[1]) + normal[2] * (z_ - point[2]), 2);
        SymEngine::expand(ex);
    }
    SymEngine::expand(ex);
    std::cout << ex << std::endl;
    SymEngine::Expression derivative_x = ex.diff(x);
    std::cout << "partial derivative of x: " << derivative_x << std::endl;
    SymEngine::Expression derivative_y = ex.diff(y);
    std::cout << "partial derivative of y: " << derivative_y << std::endl;
    SymEngine::Expression derivative_z = ex.diff(z);
    std::cout << "partial derivative of z: " << derivative_z << std::endl;
    auto solutions = SymEngine::solve_poly(derivative_x.get_basic(), x);
    std::cout << *(solutions.get()) << std::endl;
}

void linear_solver_test()
{
    // total 3 points
    std::vector<std::vector<double>> points{
        {1, 0, 0},
        {0, 2, 0},
        {0, 0, 3}
    };
    std::vector<std::vector<double>> normals{
        {1, 2, 3},
        {4, 5, 7},
        {9, 10, 3}
    };

    auto x = SymEngine::symbol("x");
    auto y = SymEngine::symbol("y");
    auto z = SymEngine::symbol("z");
    SymEngine::Expression x_(x);
    SymEngine::Expression y_(y);
    SymEngine::Expression z_(z);

    SymEngine::Expression ex;
    for (auto point_iterator = points.begin(); point_iterator != points.end(); point_iterator++)
    {
        auto point = *point_iterator;
        auto index = std::distance(points.begin(), point_iterator);
        auto normal = normals[index];
        ex = ex + pow(normal[0] * (x_ - point[0]) + normal[1] * (y_ - point[1]) + normal[2] * (z_ - point[2]), 2);
        SymEngine::expand(ex);
    }
    SymEngine::expand(ex);
    std::cout << ex << std::endl;
    SymEngine::Expression derivative_x = ex.diff(x);
    std::cout << "partial derivative of x: " << derivative_x << std::endl;
    SymEngine::Expression derivative_y = ex.diff(y);
    std::cout << "partial derivative of y: " << derivative_y << std::endl;
    SymEngine::Expression derivative_z = ex.diff(z);
    std::cout << "partial derivative of z: " << derivative_z << std::endl;
    auto solutions = SymEngine::linsolve({ derivative_x.get_basic(), derivative_y.get_basic(), derivative_z.get_basic() }, { x, y, z });
    //auto solutions = SymEngine::solve_poly_quadratic({ derivative_x.get_basic(), derivative_y.get_basic(), derivative_z.get_basic() });
    //std::cout << *(solutions.get()) << std::endl;
    for (const auto& solution : solutions)
    {
        std::cout << *(solution.get()) << std::endl;
    }

    SymEngine::Expression substituter;
    substituter = 1 - x_ - y_;
    auto boundary = ex.subs({ {z, substituter} });
    std::cout << boundary << std::endl;
}

int main()
{
    return 0;
}