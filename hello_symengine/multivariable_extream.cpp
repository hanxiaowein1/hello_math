#include <iostream>
#include <vector>
#include <stdexcept>
#include <random>
#include <limits>

#include "Eigen/Dense"
#include <symengine/expression.h>
#include "symengine/solve.h"
#include <symengine/basic.h>
#include <symengine/add.h>
#include <symengine/matrix.h>
#include <symengine/integer.h>

#include "unit_test.h"

class Polyhedron3D
{
public:
    std::vector<Eigen::Vector3d> vertices;
    std::vector<Eigen::Vector3i> triangles;
};

class Triangle2D
{
public:
    Eigen::Vector2d point_a;
    Eigen::Vector2d point_b;
    Eigen::Vector2d point_c;
};

class Line1D
{
public:

};

class ExtreamLocation
{
public:
    double value;
    std::map<SymEngine::RCP<const SymEngine::Symbol>, double, SymEngine::RCPBasicKeyLess> location;
    // Overloading the < operator for sorting
    bool operator<(const ExtreamLocation& other) const {
        return this->value < other.value; // Sort by value
    }
    void add_location(const SymEngine::RCP<const SymEngine::Symbol> &key, const double &value)
    {
        this->location.emplace(key, value);
    }
};

std::ostream& operator<<(std::ostream& os, const ExtreamLocation& m)
{
    for (const auto& [key, value] : m.location)
    {
        os << key->get_name() << ": " << value << ", ";
    }
    return os;
}


std::ostream& operator<<(std::ostream& os, const std::map<SymEngine::RCP<const SymEngine::Symbol>, double, SymEngine::RCPBasicKeyLess>& m)
{
    for (const auto& [key, value] : m)
    {
        os << key->get_name() << ": " << value << ", ";
    }
    return os;
}


double substitute_with_number(
    const SymEngine::Expression& ex,
    const std::map<SymEngine::RCP<const SymEngine::Basic>, SymEngine::RCP<const SymEngine::Basic>, SymEngine::RCPBasicKeyLess>& double_substituter
)
{
    auto substitute_result = ex.subs(double_substituter);
    auto substitute_result_double = dynamic_cast<const SymEngine::RealDouble*>(substitute_result.get_basic().get());
    if(substitute_result_double)
    {
        auto value = substitute_result_double->as_double();
        return value;
    }
    else
    {
        auto substitute_result_int = dynamic_cast<const SymEngine::Integer*>(substitute_result.get_basic().get());
        if (substitute_result_int)
        {
            return double(substitute_result_int->as_int());
        }
        else
        {
            throw std::exception("cannot subtitute expression with double");
        }
    }
}

std::tuple<bool, bool> update_min_max(
    const SymEngine::Expression& ex,
    const std::map<SymEngine::RCP<const SymEngine::Basic>, SymEngine::RCP<const SymEngine::Basic>, SymEngine::RCPBasicKeyLess>& double_substituter,
    double& desired_min, double& desired_max
)
{
    double value = substitute_with_number(ex, double_substituter);
    bool min_modified = false, max_modified = false;
    if(value > desired_max)
    {
        desired_max = value;
        max_modified = true;
    }
    if(value < desired_min)
    {
        desired_min = value;
        min_modified = true;
    }
    return {min_modified, max_modified};
}

double get_double_from_solution(const SymEngine::RCP<const SymEngine::Basic>& solution)
{
    auto number = dynamic_cast<const SymEngine::RealDouble*>(solution.get());
    if(number == nullptr)
    {
        // try to cast into integer
        auto int_number = dynamic_cast<const SymEngine::Integer*>(solution.get());
        if (int_number == nullptr)
        {
            throw std::exception("solution 3d cannot cast to double and integer");
        }
        else
        {
            return double(int_number->as_int());
        }
    }
    else
    {
        return number->as_double();
    }
}

template <typename VectorXd>
bool point_on_line(const VectorXd& point, const VectorXd& v0, const VectorXd& v1)
{
    auto distance0 = (v0 - point).norm();
    auto distance1 = (v1 - point).norm();
    auto total_distance = (v1 - v0).norm();
    if(distance0 + distance1 == total_distance)
    {
        return true;
    }
    else
    {
        return false;
    }
}

template <typename RandomType>
RandomType get_random_value(const RandomType& lower, const RandomType& upper)
{
    //using namespace std;
    std::uniform_real_distribution<RandomType> unif(lower, upper);

    std::random_device rd;
    std::default_random_engine generator(rd());
    // Getting a random double value
    RandomType random_double = unif(generator);
    return random_double;
}

bool ray_triangle_intersect(
    const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& orig,
    const Eigen::Vector3d& dir, double& tnear, double& u, double& v)
{
    // TODO: Implement this function that tests whether the triangle
    // that's specified bt v0, v1 and v2 intersects with the ray (whose
    // origin is *orig* and direction is *dir*)
    // Also don't forget to update tnear, u and v.
    auto p0 = v0;
    auto p1 = v1;
    auto p2 = v2;
    auto e1 = p1 - p0;
    auto e2 = p2 - p0;
    auto s = orig - p0;
    // auto s1 = crossProduct(dir, e2);
    auto s1 = dir.cross(e2);
    // auto s2 = crossProduct(s, e1);
    auto s2 = s.cross(e1);
    // auto s1_e1_dot = dotProduct(s1, e1);
    auto s1_e1_dot = s1.dot(e1);
    // tnear = (1.0f / s1_e1_dot) * dotProduct(s2, e2);
    tnear = (1.0f / s1_e1_dot) * s2.dot(e2);
    // u = (1.0f / s1_e1_dot) * dotProduct(s1, s);
    u = (1.0f / s1_e1_dot) * s1.dot(s);
    // v = (1.0f / s1_e1_dot) * dotProduct(s2, dir);
    v = (1.0f / s1_e1_dot) * s2.dot(dir);
    // printf("u=%.6f, v=%.6f, t=%.6f\n", u, v, tnear);
    if(u >= 0.0f && v >=0.0f && (1.0f - u - v) >=0.0f && tnear >= 0.0f){
        return true;
    }
    return false;
}

/**
 * @brief 
 * compute if point inside mesh(mesh is manifold)
 * algorithm: cast a ray from point, and count numbers of intersection, if odd, then inside mesh, from: https://stackoverflow.com/questions/17125082/3d-point-inside-object
 * @param vertices 
 * @param triangles 
 * @param orig 
 * @return true 
 * @return false 
 */
bool point_inside_mesh(
    const std::vector<Eigen::Vector3d>& vertices,
    const std::vector<Eigen::Vector3i>& triangles,
    const Eigen::Vector3d& orig
)
{
    // check if orig on vertex or edge, if on, then inside and return true
    for(const auto& vertex: vertices)
    {
        if(vertex == orig)
        {
            return true;
        }
    }

    // TODO: the following use simplest method, just iterate all triangle and test if they has intersection(if intersection lie on boundary, then drop this case, and cast a new ray)
    int intersect_count = 0;
    bool intersect_with_boundary = true;
    // if ray intersect with mesh on boundary or vertex, then try again
    while(intersect_with_boundary)
    {
        intersect_with_boundary = false;
        Eigen::Vector3d dir;
        dir[0] = get_random_value<double>(0.0f, 10.0f);
        dir[1] = get_random_value<double>(0.0f, 10.0f);
        dir[2] = get_random_value<double>(0.0f, 10.0f);
        dir = dir.normalized();
        for(const auto& triangle: triangles)
        {
            // if orig on edge, it's inside mesh, return
            std::vector<std::vector<Eigen::Vector3d>> edges = {
                {vertices[triangle.x()], vertices[triangle.y()]},
                {vertices[triangle.x()], vertices[triangle.z()]},
                {vertices[triangle.y()], vertices[triangle.z()]}
            };
            for(const auto& edge: edges)
            {
                if(point_on_line(orig, edge[0], edge[1]))
                {
                    return true;
                }
            }
            // initialize a random direction
            double tnear = 0.0f, u = 0.0f, v = 0.0f;
            if(ray_triangle_intersect(vertices[triangle.x()], vertices[triangle.y()], vertices[triangle.z()], orig, dir, tnear, u, v))
            {
                // check if point on boundary
                auto intersect_point = orig + tnear * dir;
                if(point_on_line<Eigen::Vector3d>(intersect_point, vertices[triangle.x()], vertices[triangle.y()]) || point_on_line<Eigen::Vector3d>(intersect_point, vertices[triangle.x()], vertices[triangle.z()]) || point_on_line<Eigen::Vector3d>(intersect_point, vertices[triangle.y()], vertices[triangle.z()]))
                {
                    // if on boundary, the clean intersect count
                    intersect_with_boundary = true;
                    intersect_count = 0;
                    break; // jump to try again, and get new dir
                }
                else
                {
                    intersect_count++;
                }
            }
        }
    }
    if(intersect_count % 2 == 1)
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool ray_line_intersect(const Eigen::Vector2d& v0, const Eigen::Vector2d& v1, const Eigen::Vector2d& orig, const Eigen::Vector2d& dir, double &tnear, double &b)
{
    auto symbol_t = SymEngine::symbol("t");
    auto symbol_b = SymEngine::symbol("b");
    SymEngine::Expression ex_t(symbol_t);
    SymEngine::Expression ex_b(symbol_b);
    // construct equation: O + t * D = (1 - b) V0 + b * V1;
    SymEngine::Expression equation1 = orig[0] + ex_t * dir[0] - (1 - ex_b) * v0[0] - ex_b * v1[0];
    SymEngine::Expression equation2 = orig[1] + ex_t * dir[1] - (1 - ex_b) * v0[1] - ex_b * v1[1];
    auto solutions = SymEngine::linsolve({equation1.get_basic(), equation2.get_basic()}, {symbol_t, symbol_b});
    tnear = get_double_from_solution(solutions[0]);
    b = get_double_from_solution(solutions[1]);
    if(tnear >= 0.0f && b >= 0.0f && b <= 1.0f)
    {
        return true;
    }
    return false;
}

bool point_inside_triangle(const Eigen::Vector2d& v0, const Eigen::Vector2d& v1, const Eigen::Vector2d& v2, const Eigen::Vector2d& orig)
{
    // if orig is vertex, return true
    if(orig == v0 || orig == v1 || orig == v2)
    {
        return true;
    }
    int intersect_count = 0;
    bool intersect_with_boundary = true;
    while(intersect_with_boundary)
    {
        intersect_with_boundary = false;
        Eigen::Vector2d dir;
        dir[0] = get_random_value<double>(0.0f, 10.0f);
        dir[1] = get_random_value<double>(0.0f, 10.0f);
        dir = dir.normalized();
        std::vector<std::vector<Eigen::Vector2d>> lines{
            {v0, v1}, {v0, v2}, {v1, v2}
        };
        for(const auto& line: lines)
        {
            double tnear = 0.0f, b = 0.0f;
            if(ray_line_intersect(line[0], line[1], orig, dir, tnear, b))
            {
                auto intersect_point = orig + tnear * dir;
                if(intersect_point == line[0] || intersect_point == line[1])
                {
                    intersect_with_boundary = true;
                    intersect_count = 0;
                    break;
                }
                else
                {
                    intersect_count++;
                }
            }
        }
    }
    if(intersect_count % 2 == 1)
    {
        return true;
    }
    else
    {
        return false;
    }
}

std::tuple<double, double, std::map<SymEngine::RCP<const SymEngine::Symbol>, double, SymEngine::RCPBasicKeyLess>, std::map<SymEngine::RCP<const SymEngine::Symbol>, double, SymEngine::RCPBasicKeyLess>> min_max_in_line(
    const Eigen::Vector2d& v0, const Eigen::Vector2d& v1,
    const std::vector<SymEngine::RCP<const SymEngine::Symbol>> symbol_axes,
    const SymEngine::Expression& ex/** one axis has already been substitued, composed by symbol_axes */
)
{
    double desired_max = std::numeric_limits<double>::min(), desired_min = std::numeric_limits<double>::max();

    std::map<SymEngine::RCP<SymEngine::Symbol>, double> axis_minimum_location;
    std::map<SymEngine::RCP<SymEngine::Symbol>, double> axis_maximum_location;

    std::vector<ExtreamLocation> multi_extream_location;

    // calculate line equation
    // (y1-y2) * x + (x2-x1) * y + (x1-x2)*y1 + (y2-y1)*x1 = 0
    double a = v0[1] - v1[1];
    double b = v1[0] - v0[0];
    double c = (v0[0] - v1[0]) * v0[1] + (v1[1] - v0[1]) * v0[0];
    SymEngine::Expression sub_ex0(symbol_axes[0]);
    SymEngine::Expression sub_ex1(symbol_axes[1]);
    SymEngine::Expression line_ex = a * sub_ex0 + b * sub_ex1 + c;

    std::cout << "line equation: " << line_ex << std::endl;

    auto calculate_extream_value_in_1d_line = [&](int cutted_axis, int retained_axis){
        auto cutted_symbol = symbol_axes[cutted_axis];
        auto retained_symbol = symbol_axes[retained_axis];

        auto reduce_dim_solutions = SymEngine::solve(line_ex.get_basic(), symbol_axes[cutted_axis]);
        for(const auto& reduce_dim_solution: reduce_dim_solutions->get_args())
        {
            std::cout << "reduce_dim_solution: " << *(reduce_dim_solution.get()) << std::endl;
            // one dimention: line
            auto i_d_ex = ex.subs({{cutted_symbol, reduce_dim_solution}});
            std::cout << "i_d_ex: " << i_d_ex << std::endl;
            auto derivative = i_d_ex.diff(retained_symbol);
            std::cout << "derivative: " << derivative << std::endl;
            auto i_d_solutions = SymEngine::linsolve({derivative.get_basic()}, {retained_symbol});
            std::cout << "i_d_solution: " << *(i_d_solutions[0].get()) << std::endl;
            std::cout << "type i_d_solution: " << typeid(*(i_d_solutions[0].get())).name() << std::endl;
            double i_critical_point = get_double_from_solution(i_d_solutions[0]);
            // critical point inside interval
            if(std::abs(i_critical_point - v0[retained_axis] + i_critical_point - v1[retained_axis]) < std::abs(v0[retained_axis] - v1[retained_axis]))
            {
                std::map<SymEngine::RCP<const SymEngine::Basic>, SymEngine::RCP<const SymEngine::Basic>, SymEngine::RCPBasicKeyLess> substituter{
                    {retained_symbol, SymEngine::real_double(i_critical_point)},
                };
                ExtreamLocation extream_location;
                extream_location.value = substitute_with_number(i_d_ex, substituter);
                //extream_location.location[retained_symbol] = i_critical_point;
                extream_location.location.emplace(retained_symbol, i_critical_point);
                //extream_location.location[cutted_symbol] = substitute_with_number(reduce_dim_solution, substituter);
                extream_location.location.emplace(std::make_pair(cutted_symbol, substitute_with_number(reduce_dim_solution, substituter)));
                multi_extream_location.emplace_back(std::move(extream_location));

                // auto [min_modified, max_modified] = update_min_max(i_d_ex, substituter, desired_min, desired_max);
                // if(min_modified)
                // {
                //     axis_minimum_location[retained_symbol] = i_critical_point;
                //     axis_minimum_location[cutted_symbol] = substitute_with_number(reduce_dim_solution, substituter);
                // }
                // if(max_modified)
                // {
                //     axis_maximum_location[retained_symbol] = i_critical_point;
                //     axis_maximum_location[cutted_symbol] = substitute_with_number(reduce_dim_solution, substituter);
                // }
            }
            // calculate value in endpoint 0
            {
                std::map<SymEngine::RCP<const SymEngine::Basic>, SymEngine::RCP<const SymEngine::Basic>, SymEngine::RCPBasicKeyLess> substituter1{
                    {retained_symbol, SymEngine::real_double(v0[retained_axis])},
                };
                ExtreamLocation extream_location;
                extream_location.value = substitute_with_number(i_d_ex, substituter1);
                //extream_location.location[retained_symbol] = i_critical_point;
                extream_location.location.emplace(std::make_pair(retained_symbol, v0[retained_axis]));
                //extream_location.location[cutted_symbol] = substitute_with_number(reduce_dim_solution, substituter1);
                extream_location.location.emplace(std::make_pair(cutted_symbol, substitute_with_number(reduce_dim_solution, substituter1)));
                multi_extream_location.emplace_back(std::move(extream_location));
            }

            // auto [min_modified, max_modified] = update_min_max(i_d_ex, substituter1, desired_min, desired_max);
            // if(min_modified)
            // {
            //     axis_minimum_location[retained_symbol] = i_critical_point;
            //     axis_minimum_location[cutted_symbol] = substitute_with_number(reduce_dim_solution, substituter1);
            // }
            // if(max_modified)
            // {
            //     axis_maximum_location[retained_symbol] = i_critical_point;
            //     axis_maximum_location[cutted_symbol] = substitute_with_number(reduce_dim_solution, substituter1);
            // }
            // calculate value in endpoint 1
            {
                std::map<SymEngine::RCP<const SymEngine::Basic>, SymEngine::RCP<const SymEngine::Basic>, SymEngine::RCPBasicKeyLess> substituter2{
                    {retained_symbol, SymEngine::real_double(v1[retained_axis])},
                };
                ExtreamLocation extream_location;
                extream_location.value = substitute_with_number(i_d_ex, substituter2);
                //extream_location.location[retained_symbol] = i_critical_point;
                extream_location.location.emplace(std::make_pair(retained_symbol, v1[retained_axis]));
                //extream_location.location[cutted_symbol] = substitute_with_number(reduce_dim_solution, substituter2);
                extream_location.location.emplace(std::make_pair(cutted_symbol, substitute_with_number(reduce_dim_solution, substituter2)));
                multi_extream_location.emplace_back(std::move(extream_location));
            }

            // std::tie(min_modified, max_modified) = update_min_max(i_d_ex, substituter2, desired_min, desired_max);
            // if(min_modified)
            // {
            //     axis_minimum_location[retained_symbol] = i_critical_point;
            //     axis_minimum_location[cutted_symbol] = substitute_with_number(reduce_dim_solution, substituter2);
            // }
            // if(max_modified)
            // {
            //     axis_maximum_location[retained_symbol] = i_critical_point;
            //     axis_maximum_location[cutted_symbol] = substitute_with_number(reduce_dim_solution, substituter2);
            // }
        }
    };
    if(b != 0.0f)
    {
        calculate_extream_value_in_1d_line(1, 0);
    }
    else if(a != 0.0f)
    {
        calculate_extream_value_in_1d_line(0, 1);
    }
    else
    {
        throw std::exception("both line a and b is zero, there must be problems");
    }
    // sort all possible value
    std::sort(multi_extream_location.begin(), multi_extream_location.end());
    return {multi_extream_location[0].value, multi_extream_location[multi_extream_location.size() - 1].value, multi_extream_location[0].location, multi_extream_location[multi_extream_location.size() - 1].location};
    // return {desired_min, desired_max, axis_minimum_location, axis_maximum_location};
}

std::tuple<double, double, std::map<SymEngine::RCP<const SymEngine::Symbol>, double, SymEngine::RCPBasicKeyLess>, std::map<SymEngine::RCP<const SymEngine::Symbol>, double, SymEngine::RCPBasicKeyLess>> min_max_in_triangle(
    const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2,
    // const SymEngine::RCP<const SymEngine::Symbol> &x, const SymEngine::RCP<const SymEngine::Symbol> &y, const SymEngine::RCP<const SymEngine::Symbol> &z,
    const std::vector<SymEngine::RCP<const SymEngine::Symbol>> symbol_axes,
    const SymEngine::Expression& ex
)
{
    double desired_max = std::numeric_limits<double>::min(), desired_min = std::numeric_limits<double>::max();
    std::map<SymEngine::RCP<SymEngine::Symbol>, double> axis_minimum_location;
    std::map<SymEngine::RCP<SymEngine::Symbol>, double> axis_maximum_location;
    std::vector<ExtreamLocation> multi_extream_location;

    // calculate plain equation: ax + by + cz + d = 0, (a, b, c) is normal
    // calculate the normal of triangle
    Eigen::Vector3d vec01 = v1 - v0;
    Eigen::Vector3d vec02 = v2 - v0;
    Eigen::Vector3d plane_normal = vec01.cross(vec02);
    double a = plane_normal.x(), b = plane_normal.y(), c = plane_normal.z();
    // bring in one point to calculate d
    double d = 0 - a * v0.x() - b * v0.y() - c * v0.z();

    SymEngine::Expression x_(symbol_axes[0]);
    SymEngine::Expression y_(symbol_axes[1]);
    SymEngine::Expression z_(symbol_axes[2]);
    SymEngine::Expression plane_ex = a * x_ + b * y_ + c * z_ + d;

    std::cout << "plane equation: " << plane_ex << std::endl;

    auto calculate_extream_value_in_2d_triangle = [&](int cutted_axis, int retained_axis0, int retained_axis1){
        auto cutted_symbol = symbol_axes[cutted_axis];
        auto retained_symbol0 = symbol_axes[retained_axis0];
        auto retained_symbol1 = symbol_axes[retained_axis1];


        auto reduce_dim_solutions = SymEngine::solve(plane_ex.get_basic(), symbol_axes[cutted_axis]);
        for(const auto& reduce_dim_solution : reduce_dim_solutions->get_args())
        {
            std::cout << *(reduce_dim_solution.get()) << std::endl;
            // two dimention
            auto ii_d_ex = ex.subs({{cutted_symbol, reduce_dim_solution}});
            // calculate derivative of retained axes
            auto derivative0 = ii_d_ex.diff(retained_symbol0);
            auto derivative1 = ii_d_ex.diff(retained_symbol1);
            // solve linear equation in 2d
            auto ii_d_solutions = SymEngine::linsolve({derivative0.get_basic(), derivative1.get_basic()}, {retained_symbol0, retained_symbol1});
            double ii_coor0 = get_double_from_solution(ii_d_solutions[0]), ii_coor1 = get_double_from_solution(ii_d_solutions[1]);
            Eigen::Vector2d ii_critical_point{ii_coor0, ii_coor1};
            Eigen::Vector2d ii_d_v0{v0[retained_axis0], v0[retained_axis1]}, ii_d_v1{v1[retained_axis0], v1[retained_axis1]}, ii_d_v2{v2[retained_axis0], v2[retained_axis1]};
            if(point_inside_triangle(ii_d_v0, ii_d_v1, ii_d_v2, ii_critical_point))
            {
                std::map<SymEngine::RCP<const SymEngine::Basic>, SymEngine::RCP<const SymEngine::Basic>, SymEngine::RCPBasicKeyLess> substituter{
                    {retained_symbol0, SymEngine::real_double(ii_critical_point.x())},
                    {retained_symbol1, SymEngine::real_double(ii_critical_point.y())},
                };
                // auto substitute_result = ii_d_ex.subs(substituter);
                // auto substitute_result_double = dynamic_cast<const SymEngine::RealDouble*>(substitute_result.get_basic().get());
                ExtreamLocation extream_location;
                extream_location.value = substitute_with_number(ii_d_ex, substituter);
                //extream_location.location[retained_symbol0] = ii_critical_point.x();
                extream_location.location.emplace(std::make_pair(retained_symbol0, ii_critical_point.x()));
                //extream_location.location[retained_symbol1] = ii_critical_point.y();
                extream_location.location.emplace(std::make_pair(retained_symbol1, ii_critical_point.y()));
                //extream_location.location[cutted_symbol] = substitute_with_number(reduce_dim_solution, substituter);
                extream_location.location.emplace(std::make_pair(cutted_symbol, substitute_with_number(reduce_dim_solution, substituter)));
                multi_extream_location.emplace_back(std::move(extream_location));

                // auto [min_modified, max_modified] = update_min_max(ii_d_ex, substituter, desired_min, desired_max);
                // if(min_modified)
                // {
                //     axis_minimum_location[retained_symbol0] = ii_critical_point.x();
                //     axis_minimum_location[retained_symbol1] = ii_critical_point.y();
                //     axis_minimum_location[cutted_symbol] = substitute_with_number(reduce_dim_solution, substituter);
                // }
                // if(max_modified)
                // {
                //     axis_maximum_location[retained_symbol0] = ii_critical_point.x();
                //     axis_maximum_location[retained_symbol1] = ii_critical_point.y();
                //     axis_maximum_location[cutted_symbol] = substitute_with_number(reduce_dim_solution, substituter);
                // }
            }
            // calculate min max in line
            // std::vector<double> min_array, max_array;
            std::vector<std::vector<Eigen::Vector2d>> lines = {
                {ii_d_v0, ii_d_v1},
                {ii_d_v0, ii_d_v2},
                {ii_d_v1, ii_d_v2},
            };
            for(const auto& line: lines)
            {
                auto [min_value, max_value, ii_d_minimum_location, ii_d_maximum_location] = min_max_in_line(line[0], line[1], {retained_symbol0, retained_symbol1}, ii_d_ex);
                auto generate_new_extream_location = [&](double &value, std::map<SymEngine::RCP<const SymEngine::Symbol>, double, SymEngine::RCPBasicKeyLess>& location){
                    ExtreamLocation extream_location;
                    extream_location.value = value;
                    extream_location.location = location;
                    std::map<SymEngine::RCP<const SymEngine::Basic>, SymEngine::RCP<const SymEngine::Basic>, SymEngine::RCPBasicKeyLess> substituter;
                    for(const auto& [key, value]: extream_location.location)
                    {
                        substituter[key] = SymEngine::real_double(value);
                    }
                    //extream_location.location[cutted_symbol] = substitute_with_number(reduce_dim_solution, substituter);
                    extream_location.location.emplace(std::make_pair(cutted_symbol, substitute_with_number(reduce_dim_solution, substituter)));
                    multi_extream_location.emplace_back(extream_location);
                };
                generate_new_extream_location(min_value, ii_d_minimum_location);
                generate_new_extream_location(max_value, ii_d_maximum_location);
                // min_array.emplace_back(min_value);
                // max_array.emplace_back(max_value);
            }
            // desired_min = std::min(desired_min, *std::min_element(min_array.begin(), min_array.end()));
            // desired_max = std::max(desired_max, *std::max_element(max_array.begin(), max_array.end()));
        }
    };

    if(c != 0.0f)
    {
        // substitute z
        calculate_extream_value_in_2d_triangle(2, 0, 1);
    }
    else if(b != 0.0f)
    {
        calculate_extream_value_in_2d_triangle(1, 0, 2);
    }
    else if(a != 0.0f)
    {
        calculate_extream_value_in_2d_triangle(0, 1, 2);
    }
    else
    {
        throw std::exception("all axis in normal is zero, there are must problems!");
    }
    return {multi_extream_location[0].value, multi_extream_location[multi_extream_location.size() - 1].value, multi_extream_location[0].location, multi_extream_location[multi_extream_location.size() - 1].location};
    // return {desired_min, desired_max};
}

// write code from beginning to end
std::tuple<double, double, std::map<SymEngine::RCP<const SymEngine::Symbol>, double, SymEngine::RCPBasicKeyLess>, std::map<SymEngine::RCP<const SymEngine::Symbol>, double, SymEngine::RCPBasicKeyLess>> min_max_3_variable_quadratic_polynomial(
    const SymEngine::Expression& equation,
    const std::vector<SymEngine::RCP<const SymEngine::Symbol>>& symbol_axes,
    const std::vector<Eigen::Vector3d>& vertices,
    const std::vector<Eigen::Vector3i>& triangles
)
{
    // calculate deriative of x, y, z
    SymEngine::Expression derivative_x = equation.diff(symbol_axes[0]);
    std::cout << "partial derivative of x: " << derivative_x << std::endl;
    SymEngine::Expression derivative_y = equation.diff(symbol_axes[1]);
    std::cout << "partial derivative of y: " << derivative_y << std::endl;
    SymEngine::Expression derivative_z = equation.diff(symbol_axes[2]);
    std::cout << "partial derivative of z: " << derivative_z << std::endl;

    // calculate 3d extream point
    auto solutions_3d = SymEngine::linsolve({ derivative_x.get_basic(), derivative_y.get_basic(), derivative_z.get_basic() }, { symbol_axes[0], symbol_axes[1], symbol_axes[2] });
    Eigen::Vector3d extream_3d_point;
    extream_3d_point[0] = get_double_from_solution(solutions_3d[0]);
    extream_3d_point[1] = get_double_from_solution(solutions_3d[1]);
    extream_3d_point[2] = get_double_from_solution(solutions_3d[2]);

    //auto get_3d_coor_func = [&](int index) {
    //    auto solution_3d_coor = solutions_3d[index];
    //    auto number = dynamic_cast<const SymEngine::RealDouble*>(solution_3d_coor.get());
    //    if(number == nullptr)
    //    {
    //        throw std::exception("solution 3d cannot cast to double");
    //    }
    //    else
    //    {
    //        extream_3d_point[index] = number->as_double();
    //    }
    //};
    //get_3d_coor_func(0);
    //get_3d_coor_func(1);
    //get_3d_coor_func(2);

    // TODO: this is fake 3d mesh, should be replace in the future
    // check if 3d extream point is in 3d range
    // std::vector<Eigen::Vector3d> vertices;
    // std::vector<Eigen::Vector3i> triangles;
    double desired_max = std::numeric_limits<double>::min(), desired_min = std::numeric_limits<double>::max();
    std::vector<ExtreamLocation> multi_extream_location;
    if(point_inside_mesh(vertices, triangles, extream_3d_point))
    {
        // get extream value and substitute max and min
        std::map<SymEngine::RCP<const SymEngine::Basic>, SymEngine::RCP<const SymEngine::Basic>, SymEngine::RCPBasicKeyLess> substituter{
            {symbol_axes[0], SymEngine::real_double(extream_3d_point.x())},
            {symbol_axes[1], SymEngine::real_double(extream_3d_point.y())},
            {symbol_axes[2], SymEngine::real_double(extream_3d_point.z())},
        };
        // auto substitute_result = equation.subs(substituter);
        // auto substitute_result_double = dynamic_cast<const SymEngine::RealDouble*>(substitute_result.get_basic().get());
        ExtreamLocation extream_location;
        extream_location.value = substitute_with_number(equation, substituter);
        //extream_location.location[symbol_axes[0]] = extream_3d_point.x();
        extream_location.location.emplace(std::make_pair(symbol_axes[0], extream_3d_point.x()));
        //extream_location.location[symbol_axes[1]] = extream_3d_point.y();
        extream_location.location.emplace(std::make_pair(symbol_axes[1], extream_3d_point.y()));
        //extream_location.location[symbol_axes[2]] = extream_3d_point.z();
        extream_location.location.emplace(std::make_pair(symbol_axes[2], extream_3d_point.z()));
        multi_extream_location.emplace_back(std::move(extream_location));
        // if(substitute_result_double)
        // {
        //     auto value = substitute_result_double->as_double();
        //     if(value > desired_max)
        //     {
        //         desired_max = value;
        //     }
        //     if(value < desired_min)
        //     {
        //         desired_min = value;
        //     }
        // }
        // else
        // {
        //     throw std::exception("cannot subtitute expression with double");
        // }
    }

    // get max/min value from 3d boundary
    // for triangle in triangles
    // compute min/max value from each triangle
    // std::vector<double> min_vec, max_vec;
    for(const auto& triangle: triangles)
    {
        auto [min_value, max_value, min_location, max_location] = min_max_in_triangle(vertices[triangle[0]], vertices[triangle[1]], vertices[triangle[2]], {symbol_axes[0], symbol_axes[1], symbol_axes[2]}, equation);
        auto generate_new_extream_location = [&](double &value, std::map<SymEngine::RCP<const SymEngine::Symbol>, double, SymEngine::RCPBasicKeyLess>& location){
            ExtreamLocation extream_location;
            extream_location.value = value;
            extream_location.location = location;
            multi_extream_location.emplace_back(extream_location);
        };
        generate_new_extream_location(min_value, min_location);
        generate_new_extream_location(max_value, max_location);
        // min_vec.emplace_back(min_value);
        // max_vec.emplace_back(max_value);
    }
    // desired_min = std::min(desired_min, *std::min_element(min_vec.begin(), min_vec.end()));
    // desired_max = std::max(desired_max, *std::max_element(max_vec.begin(), max_vec.end()));
    // return { desired_min, desired_max };
    return {multi_extream_location[0].value, multi_extream_location[multi_extream_location.size() - 1].value, multi_extream_location[0].location, multi_extream_location[multi_extream_location.size() - 1].location};
}

TEST(GlobalTest, update_min_max)
{
    double desired_min = 5, desired_max = 1;
    auto x = SymEngine::symbol("x");
    auto y = SymEngine::symbol("y");
    auto z = SymEngine::symbol("z");

    SymEngine::Expression x_(x);
    SymEngine::Expression y_(y);
    SymEngine::Expression z_(z);

    SymEngine::Expression ex;
    ex = x_ + y_ + z_;

    std::map<SymEngine::RCP<const SymEngine::Basic>, SymEngine::RCP<const SymEngine::Basic>, SymEngine::RCPBasicKeyLess> substituter{
        {x, SymEngine::real_double(1)},
        {y, SymEngine::real_double(1)},
        {z, SymEngine::real_double(1)},
    };

    update_min_max(ex, substituter, desired_min, desired_max);
    ASSERT_EQ(3.0f, desired_min);
    ASSERT_EQ(3.0f, desired_max);
}

TEST(GlobalTest, get_double_from_solution)
{
    auto x = SymEngine::symbol("x");
    SymEngine::Expression x_(x);
    SymEngine::Expression ex;
    ex = ex + x_ - 1.0f;
    auto solutions = SymEngine::linsolve({ ex.get_basic() }, { x });
    for (const auto& solution : solutions)
    {
        double result = get_double_from_solution(solution);
        ASSERT_EQ(result, 1.0f);
    }
}

TEST(GlobalTest, point_on_line)
{
    // test point on line
    {
        Eigen::Vector3d point0{ 0, 0, 0 };
        Eigen::Vector3d point1{ 1, 1, 1 };
        Eigen::Vector3d point2{ 2, 2, 2 };
        auto result = point_on_line<Eigen::Vector3d>(point1, point0, point2);
        ASSERT_EQ(true, result);
    }
    // test point not on line
    {
        Eigen::Vector3d point0{ 0, 0, 0 };
        Eigen::Vector3d point1{ 1, 1, 2 };
        Eigen::Vector3d point2{ 2, 2, 2 };
        auto result = point_on_line<Eigen::Vector3d>(point1, point0, point2);
        ASSERT_EQ(false, result);
    }
}

TEST(GlobalTest, get_random_value)
{
    double result = get_random_value<double>(1.0f, 10.0f);
    EXPECT_TRUE((result >= 1.0f) && (result <= 10.0f));
}

TEST(GlobalTest, ray_triangle_intersect)
{
    // test intersect
    {
        Eigen::Vector3d v0{ 1, 0, 0 };
        Eigen::Vector3d v1{ 0, 1, 0 };
        Eigen::Vector3d v2{ 0, 0, 1 };

        Eigen::Vector3d orig{ 0,0,0 };
        Eigen::Vector3d dir{ 1,1,1 };

        double tnear = 0.0f, u = 0.0f, v = 0.0f;
        bool is_intersect = ray_triangle_intersect(v0, v1, v2, orig, dir, tnear, u, v);
        ASSERT_EQ(true, is_intersect);
    }

    // test not intersect(reverse direction)
    {
        Eigen::Vector3d v0{ 1, 0, 0 };
        Eigen::Vector3d v1{ 0, 1, 0 };
        Eigen::Vector3d v2{ 0, 0, 1 };

        Eigen::Vector3d orig{ 0,0,0 };
        Eigen::Vector3d dir{ -1,-1,-1 };

        double tnear = 0.0f, u = 0.0f, v = 0.0f;
        bool is_intersect = ray_triangle_intersect(v0, v1, v2, orig, dir, tnear, u, v);
        ASSERT_EQ(false, is_intersect);
    }

    // test not intersect
    {
        Eigen::Vector3d v0{ 1, 0, 0 };
        Eigen::Vector3d v1{ 0, 1, 0 };
        Eigen::Vector3d v2{ 0, 0, 1 };

        Eigen::Vector3d orig{ 2,2,2 };
        Eigen::Vector3d dir{ 1,1,1 };

        double tnear = 0.0f, u = 0.0f, v = 0.0f;
        bool is_intersect = ray_triangle_intersect(v0, v1, v2, orig, dir, tnear, u, v);
        ASSERT_EQ(false, is_intersect);
    }
}

TEST(GlobalTest, point_inside_mesh)
{
    Eigen::Vector3d v0{ 0, 0, 0 };
    Eigen::Vector3d v1{ 1, 0, 0 };
    Eigen::Vector3d v2{ 0, 1, 0 };
    Eigen::Vector3d v3{ 0, 0, 1 };
    std::vector<Eigen::Vector3d> vertices = {
        v0, v1, v2, v3
    };
    std::vector<Eigen::Vector3i> triangles = {
        Eigen::Vector3i{0, 1, 2},
        Eigen::Vector3i{0, 1, 3},
        Eigen::Vector3i{0, 2, 3},
        Eigen::Vector3i{1, 2, 3},
    };

    // inside
    {
        Eigen::Vector3d orig = { 0.1, 0.1, 0.1 };
        bool result = point_inside_mesh(vertices, triangles, orig);
        ASSERT_EQ(true, result);
    }

    // outside
    {
        Eigen::Vector3d orig = { 2, 2, 2 };
        bool result = point_inside_mesh(vertices, triangles, orig);
        ASSERT_EQ(false, result);
    }
}

TEST(GlobalTest, ray_line_intersect)
{
    Eigen::Vector2d v0{ 1, 0 };
    Eigen::Vector2d v1{ 0, 1 };

    double tnear = 0.0f, b = 0.0f;

    // intersect
    {
        Eigen::Vector2d orig{ 0, 0 };
        Eigen::Vector2d dir{ 1, 1 };
        dir = dir.normalized();

        bool result = ray_line_intersect(v0, v1, orig, dir, tnear, b);
        ASSERT_EQ(true, result);
    }

    // not intersect(reverse direction)
    {
        Eigen::Vector2d orig{ 0, 0 };
        Eigen::Vector2d dir{ -1, -1 };
        dir = dir.normalized();

        bool result = ray_line_intersect(v0, v1, orig, dir, tnear, b);
        ASSERT_EQ(false, result);
    }

    // not intersect
    {
        Eigen::Vector2d orig{ 2, 2 };
        Eigen::Vector2d dir{ 1, 1 };
        dir = dir.normalized();

        bool result = ray_line_intersect(v0, v1, orig, dir, tnear, b);
        ASSERT_EQ(false, result);
    }
}

TEST(GlobalTest, point_inside_triangle)
{
    Eigen::Vector2d v0{ 0, 0 };
    Eigen::Vector2d v1{ 1, 0 };
    Eigen::Vector2d v2{ 0, 1 };

    // inside
    {
        Eigen::Vector2d orig{ 0.4f, 0.4f };
        bool result = point_inside_triangle(v0, v1, v2, orig);
        ASSERT_EQ(true, result);
    }

    // outside
    {
        Eigen::Vector2d orig{ 0.6f, 0.7f };
        bool result = point_inside_triangle(v0, v1, v2, orig);
        ASSERT_EQ(false, result);
    }
}

TEST(GlobalTest, min_max_in_line)
{
    Eigen::Vector2d v0{ 1,0 };
    Eigen::Vector2d v1{ 0,1 };
    
    auto x = SymEngine::symbol("x");
    auto y = SymEngine::symbol("y");

    SymEngine::Expression x_(x);
    SymEngine::Expression y_(y);

    SymEngine::Expression ex = SymEngine::pow(x_, 2) + SymEngine::pow(y_, 2);
    auto [min_value, max_value, min_location, max_location] = min_max_in_line(v0, v1, { x, y }, ex);
    ASSERT_EQ(0.5f, min_value);
    ASSERT_EQ(min_location.at(x), 0.5f);
    ASSERT_EQ(1.0f, max_value);
}

TEST(GlobalTest, min_max_in_triangle)
{
    double double_deviation = 0.00000000000000010;
    Eigen::Vector3d v0{ 1, 0, 0 };
    Eigen::Vector3d v1{ 0, 1, 0 };
    Eigen::Vector3d v2{ 0, 0, 1 };

    auto x = SymEngine::symbol("x");
    auto y = SymEngine::symbol("y");
    auto z = SymEngine::symbol("z");

    SymEngine::Expression x_(x);
    SymEngine::Expression y_(y);
    SymEngine::Expression z_(z);

    SymEngine::Expression ex = SymEngine::pow(x_, 2) + SymEngine::pow(y_, 2) + SymEngine::pow(z_, 2);
    auto [min_value, max_value, min_location, max_location] = min_max_in_triangle(v0, v1, v2, { x, y, z }, ex);
    //ASSERT_EQ(double(1) / double(3), min_value);
    EXPECT_TRUE(std::abs(double(1) / double(3) - min_value) < double_deviation);
    ASSERT_EQ(1, max_value);
    std::cout << min_location << std::endl;
    std::cout << max_location << std::endl;
}

TEST(GlobalTest, min_max_3_variable_quadratic_polynomial)
{
    std::cout << "----------------------------test min_max_3_variable_quadratic_polynomial--------------------------" << std::endl;
    // compute MIN/MAX SUM((N*(x - A))^2)
    // interpolated points and its normals
    std::vector<Eigen::Vector3d> interpolated_points = {
        {0.5f, 0, 0},
        {0, 0.5f, 0},
        {0, 0, 0.5f},
    };
    std::vector<Eigen::Vector3d> interpolated_normals = {
        {0, 1, 2},
        {2, 0, 1},
        {1, 2, 0},
    };

    // get expression
    auto x = SymEngine::symbol("x");
    auto y = SymEngine::symbol("y");
    auto z = SymEngine::symbol("z");
    SymEngine::Expression x_(x);
    SymEngine::Expression y_(y);
    SymEngine::Expression z_(z);
    SymEngine::Expression ex;
    for(auto iter = interpolated_points.begin(); iter != interpolated_points.end(); iter++)
    {
        auto index = std::distance(interpolated_points.begin(), iter);
        auto interpolated_point = *iter;
        auto interpolated_normal = interpolated_normals[index];
        ex = ex + SymEngine::pow(interpolated_normal.x() * (x_ - interpolated_point.x()) + interpolated_normal.y() * (y_ - interpolated_point.y()) + interpolated_normal.z() * (z_ - interpolated_point.z()), 2);
    }
    std::cout << ex << std::endl;

    std::vector<Eigen::Vector3d> vertices = {
        {0, 0, 0},
        {1, 0, 0},
        {0, 1, 0},
        {0, 0, 1},
    };
    std::vector<Eigen::Vector3i> triangles = {
        {0, 1, 2},
        {0, 1, 3},
        {0, 2, 3},
        {1, 2, 3},
    };

    auto [min_value, max_value, min_location, max_location] = min_max_3_variable_quadratic_polynomial(ex, {x, y, z}, vertices, triangles);
    std::cout << "min_value" << ": " << min_value << std::endl;
    std::cout << min_location << std::endl;
    // std::cout << "min location" << std::endl;
    // std::cout << "x: " << min_location.at(x);
    // std::cout << "y: " << min_location.at(y);
    // std::cout << "z: " << min_location.at(z);

    std::cout << "max_value" << ": " << max_value << std::endl;
    std::cout << "max location" << std::endl;
    std::cout << max_location << std::endl;
    // std::cout << "x: " << max_location.at(x);
    // std::cout << "y: " << max_location.at(y);
    // std::cout << "z: " << max_location.at(z);
}