#include <iostream>
#include <vector>
#include <symengine/expression.h>
#include "symengine/solve.h"
#include <symengine/basic.h>
#include <symengine/add.h>
#include <symengine/matrix.h>
#include <symengine/integer.h>

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
        auto value = solution.get();
        // auto value1 = value->loads;
        // std::cout << value1 << std::endl;
        //std::cout << typeid(value).name() << std::endl;
        //auto solution_type = typeid(*value).name();
        if (typeid(*value) == typeid(SymEngine::RealDouble))
        {
            std::cout << "is real double" << std::endl;
        }
        auto number = dynamic_cast<const SymEngine::RealDouble*>(value);
        if(number == nullptr)
        {
            std::cout << "cast failed" << std::endl;
        }
        else
        {
            std::cout << number->as_double() << std::endl;
        }
        std::cout << *(solution.get()) << std::endl;
    }

    SymEngine::Expression substituter;
    substituter = 1 - x_ - y_;
    auto boundary = ex.subs({ {z, substituter} });
    std::cout << boundary << std::endl;
}

void linear_solver_test2()
{
    auto x = SymEngine::symbol("x");
    auto y = SymEngine::symbol("y");
    auto z = SymEngine::symbol("z");

    SymEngine::Expression x_(x);
    SymEngine::Expression y_(y);
    SymEngine::Expression z_(z);
    SymEngine::Expression ex;
    ex = ex + 1 * x_ + 0 * y_ + 0 * z_ - 1;
    SymEngine::expand(ex);
    // auto solutions = SymEngine::solve_poly(ex.get_basic(), x);
    auto solutions = SymEngine::solve(ex.get_basic(), x);
    auto integers_set = dynamic_cast<const SymEngine::Integers*>(solutions.get());
    if(integers_set)
    {
        std::cout << "cast succeed" << std::endl;
    }
    std::cout << typeid(*(solutions.get())).name() << std::endl;
    std::cout << *(solutions.get()) << std::endl;
    SymEngine::expand(solutions);
    auto solutions_args = solutions->get_args();
    std::cout << solutions_args << std::endl;

    //NOTE: not implement
    //auto solutions_ex = solutions->expand_as_exp();
    //std::cout << *(solutions_ex.get()) << std::endl;

    SymEngine::Expression ex2;
    ex2 = ex2 + 1 * x_ + 1 * y_ + 1 * z_ - 1;
    SymEngine::expand(ex2);
    auto equation1 = ex2.subs({{x, solutions}});
    SymEngine::expand(equation1);
    std::cout << equation1 << std::endl;

    SymEngine::Expression ex3;
    ex3 = ex3 + x_ - y_ + z_;
    SymEngine::expand(ex3);
    auto equation2 = ex3.subs({ { x, solutions } });
    
    auto solutions2 = SymEngine::linsolve({ equation1.get_basic(), equation2.get_basic() }, { y, z });
    for (const auto& solution : solutions2)
    {
        SymEngine::expand(solution);
        std::cout << typeid(*(solution.get())).name() << std::endl;
        std::cout << *(solution.get()) << std::endl;
    }

}

void linear_solver_test3()
{
    auto x = SymEngine::symbol("x");
    auto y = SymEngine::symbol("y");
    auto z = SymEngine::symbol("z");

    SymEngine::Expression x_(x);
    SymEngine::Expression y_(y);
    SymEngine::Expression z_(z);
    SymEngine::Expression ex;
    ex = ex + 1 * x_ + 0 * y_ + 0 * z_ - 1;
    SymEngine::expand(ex);
    // auto solutions = SymEngine::solve_poly(ex.get_basic(), x);
    auto solutions = SymEngine::solve(ex.get_basic(), x);
    auto solutions_args = solutions->get_args();
    std::cout << solutions_args << std::endl;
    for(const auto& solutions_arg: solutions_args)
    {
        SymEngine::Expression ex2;
        ex2 = ex2 + 1 * x_ + 1 * y_ + 1 * z_ - 1;
        SymEngine::expand(ex2);
        auto equation1 = ex2.subs({{x, solutions_arg}});
        SymEngine::expand(equation1);
        std::cout << equation1 << std::endl;

        SymEngine::Expression ex3;
        ex3 = ex3 + x_ - y_ + z_;
        SymEngine::expand(ex3);
        auto equation2 = ex3.subs({ { x, solutions_arg } });
        
        auto solutions2 = SymEngine::linsolve({ equation1.get_basic(), equation2.get_basic() }, { y, z });
        for (const auto& solution : solutions2)
        {
            SymEngine::expand(solution);
            std::cout << typeid(*(solution.get())).name() << std::endl;
            std::cout << *(solution.get()) << std::endl;
        }
    }
}

void substitute_with_normal_value()
{
    auto x = SymEngine::symbol("x");
    auto y = SymEngine::symbol("y");
    auto z = SymEngine::symbol("z");

    SymEngine::Expression x_(x);
    SymEngine::Expression y_(y);
    SymEngine::Expression z_(z);
    SymEngine::Expression ex;
    ex = ex + 1 * x_ + 1 * y_ + 1 * z_ - 1;
    std::map<SymEngine::RCP<const SymEngine::Basic>, SymEngine::RCP<const SymEngine::Basic>, SymEngine::RCPBasicKeyLess> substituter;
    substituter[x] = SymEngine::integer(3);
    substituter[y] = SymEngine::integer(5);
    substituter[z] = SymEngine::integer(5);
    auto subs_ex = ex.subs(substituter);
    std::cout << subs_ex << std::endl;
    auto subs_double = dynamic_cast<const SymEngine::Integer*>(subs_ex.get_basic().get());
    if(subs_double)
    {
        std::cout << subs_double->as_int() << std::endl;
    }
}

int main()
{
    //linear_solver_test2();
    auto x = SymEngine::symbol("x");
    auto y = SymEngine::symbol("y");
    auto z = SymEngine::symbol("z");

    SymEngine::Expression x_(x);
    SymEngine::Expression y_(y);
    SymEngine::Expression z_(z);

    SymEngine::Expression ex = x_ + y_ + z_ + 1;
    auto solution = SymEngine::solve(ex.get_basic(), z);
    std::cout << *(solution.get()) << std::endl;

    return 0;
}