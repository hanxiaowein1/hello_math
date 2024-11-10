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

void update_min_max(const SymEngine::Expression& ex, const std::map<SymEngine::RCP<const SymEngine::Basic>, SymEngine::RCP<const SymEngine::Basic>, SymEngine::RCPBasicKeyLess>& double_substituter, double& desired_min, double& desired_max)
{
    auto substitute_result = ex.subs(double_substituter);
    auto substitute_result_double = dynamic_cast<const SymEngine::RealDouble*>(substitute_result.get_basic().get());
    if(substitute_result_double)
    {
        auto value = substitute_result_double->as_double();
        if(value > desired_max)
        {
            desired_max = value;
        }
        if(value < desired_min)
        {
            desired_min = value;
        }
    }
    else
    {
        throw std::exception("cannot subtitute expression with double");
    }
}

double get_double_from_solution(const SymEngine::RCP<const SymEngine::Basic>& solution)
{
    auto number = dynamic_cast<const SymEngine::RealDouble*>(solution.get());
    if(number == nullptr)
    {
        throw std::exception("solution 3d cannot cast to double");
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
    if(distance0 == distance1)
    {
        return true;
    }
    else
    {
        return false;
    }
}

template <typename RandomType>
RandomType get_random_value(const RandomType& upper, const RandomType& lower)
{
    using namespace std;
    uniform_real_distribution<RandomType> unif(lower_bound, upper_bound);

    default_random_engine re;
    // Getting a random double value
    RandomType random_double = unif(re);
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
    // TODO: the following use simplest method, just iterate all triangle and test if they has intersection(if intersection lie on boundary, then drop this case, and cast a new ray)

    int intersect_count = 0;
    // if ray intersect with mesh on boundary or vertex, then try again
    while(true)
    {
        Eigen::Vector3d dir;
        dir[0] = get_random_value<double>(0.0f, 10.0f);
        dir[1] = get_random_value<double>(0.0f, 10.0f);
        dir[2] = get_random_value<double>(0.0f, 10.0f);
        dir = dir.normalized();
        for(const auto& triangle: triangles)
        {
            // initialize a random direction
            double tnear = 0.0f, u = 0.0f, v = 0.0f;
            if(ray_triangle_intersect(vertices[triangle.x()], vertices[triangle.y()], vertices[triangle.z()], orig, dir, tnear, u, v))
            {
                // check if point on boundary
                auto intersect_point = orig + tnear * dir;
                if(point_on_line<Eigen::Vector3d>(intersect_point, vertices[triangle.x()], vertices[triangle.y()]) || point_on_line<Eigen::Vector3d>(intersect_point, vertices[triangle.x()], vertices[triangle.z()]) || point_on_line<Eigen::Vector3d>(intersect_point, vertices[triangle.y()], vertices[triangle.z()]))
                {
                    // if on boundary, the clean intersect count
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
    // construct equation: O + t * D = (1 - b) V0 + V1;
    SymEngine::Expression equation1 = orig[0] + ex_t * dir[0] - (1 - ex_b) * v0[0] + v1[0];
    SymEngine::Expression equation2 = orig[1] + ex_t * dir[1] - (1 - ex_b) * v0[1] + v1[1];
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
    int intersect_count = 0;
    while(true)
    {
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

std::pair<double, double> min_max_in_line(
    const Eigen::Vector2d& v0, const Eigen::Vector2d& v1,
    const std::vector<SymEngine::RCP<const SymEngine::Symbol>> symbol_axes,
    const SymEngine::Expression& ex/** one axis has already been substitued, composed by symbol_axes */
)
{
    double desired_max = std::numeric_limits<double>::min(), desired_min = std::numeric_limits<double>::max();

    // calculate line equation
    // (y1-y2) * x + (x2-x1) * y + (x1-x2)*y1 + (y2-y1)*x1 = 0
    double a = v0[1] - v1[1];
    double b = v1[0] - v0[0];
    double c = (v0[0] - v1[0]) * v0[1] + (v1[1] - v0[1]) * v0[0];
    SymEngine::Expression sub_ex0(symbol_axes[0]);
    SymEngine::Expression sub_ex1(symbol_axes[1]);
    SymEngine::Expression line_ex = a * sub_ex0 + b * sub_ex1 + c;
    auto calculate_extream_value_in_1d_line = [&](int cutted_axis, int retained_axis){
        auto cutted_symbol = symbol_axes[cutted_axis];
        auto retained_symbol = symbol_axes[retained_axis];

        auto reduce_dim_solutions = SymEngine::solve(line_ex.get_basic(), symbol_axes[cutted_axis]);
        for(const auto& reduce_dim_solution: reduce_dim_solutions->get_args())
        {
            // one dimention: line
            auto i_d_ex = ex.subs({{cutted_symbol, reduce_dim_solution}});
            auto derivative = i_d_ex.diff(retained_symbol);
            auto i_d_solutions = SymEngine::linsolve({derivative}, {retained_symbol});
            double i_critical_point = get_double_from_solution(i_d_solutions[0]);
            // critical point inside interval
            if(std::abs(i_critical_point - v0[retained_axis] + i_critical_point - v1[retained_axis]) < std::abs(v0[retained_axis] - v1[retained_axis]))
            {
                std::map<SymEngine::RCP<const SymEngine::Basic>, SymEngine::RCP<const SymEngine::Basic>, SymEngine::RCPBasicKeyLess> substituter{
                    {retained_symbol, SymEngine::real_double(i_critical_point)},
                };
                update_min_max(i_d_ex, substituter, desired_min, desired_max);
            }
            // calculate value in endpoint 0
            std::map<SymEngine::RCP<const SymEngine::Basic>, SymEngine::RCP<const SymEngine::Basic>, SymEngine::RCPBasicKeyLess> substituter{
                {retained_symbol, SymEngine::real_double(v0[retained_axis])},
            };
            update_min_max(i_d_ex, substituter, desired_min, desired_max);
            // calculate value in endpoint 1
            std::map<SymEngine::RCP<const SymEngine::Basic>, SymEngine::RCP<const SymEngine::Basic>, SymEngine::RCPBasicKeyLess> substituter{
                {retained_symbol, SymEngine::real_double(v1[retained_axis])},
            };
            update_min_max(i_d_ex, substituter, desired_min, desired_max);
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
    return std::make_pair(desired_min, desired_max);
}

std::pair<double, double> min_max_in_triangle(
    const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2,
    // const SymEngine::RCP<const SymEngine::Symbol> &x, const SymEngine::RCP<const SymEngine::Symbol> &y, const SymEngine::RCP<const SymEngine::Symbol> &z,
    const std::vector<SymEngine::RCP<const SymEngine::Symbol>> symbol_axes,
    const SymEngine::Expression& ex
)
{
    double desired_max = std::numeric_limits<double>::min(), desired_min = std::numeric_limits<double>::max();

    // calculate plain equation: ax + by + cz + d = 0, (a, b, c) is normal
    // calculate the normal of triangle
    Eigen::Vector3d vec01 = v1 - v0;
    Eigen::Vector3d vec02 = v2 - v0;
    Eigen::Vector3d plane_normal = vec01.cross(vec02);
    double a = plane_normal.x(), b = plane_normal.y(), c = plane_normal.z();
    // bring in one point to calculate d
    double d = 0 - a * v0.x() - b * v1.y() - c * v1.z();

    SymEngine::Expression x_(symbol_axes[0]);
    SymEngine::Expression y_(symbol_axes[1]);
    SymEngine::Expression z_(symbol_axes[2]);
    SymEngine::Expression plane_ex = a * x_ + b * y_ + c * z_ + d;

    auto calculate_extream_value_in_2d_triangle = [&](int cutted_axis, int retained_axis0, int retained_axis1){
        auto cutted_symbol = symbol_axes[cutted_axis];
        auto retained_symbol0 = symbol_axes[retained_axis0];
        auto retained_symbol1 = symbol_axes[retained_axis1];

        auto reduce_dim_solutions = SymEngine::solve(plane_ex.get_basic(), symbol_axes[cutted_axis]);
        for(const auto& reduce_dim_solution : reduce_dim_solutions->get_args())
        {
            // two dimention
            auto ii_d_ex = ex.subs({{cutted_symbol, reduce_dim_solution}});
            // calculate derivative of retained axes
            auto derivative0 = ii_d_ex.diff(retained_symbol0);
            auto derivative1 = ii_d_ex.diff(retained_symbol1);
            // solve linear equation in 2d
            auto ii_d_solutions = SymEngine::linsolve({derivative0.get_basic(), derivative1.get_basic()}, {retained_symbol0, retained_symbol0});
            double ii_coor0 = get_double_from_solution(ii_d_solutions[0]), ii_coor1 = get_double_from_solution(ii_d_solutions[1]);
            Eigen::Vector2d ii_critical_point{ii_coor0, ii_coor1};
            Eigen::Vector2d ii_d_v0{v0[retained_axis0], v0[retained_axis1]}, ii_d_v1{v1[retained_axis0], v1[retained_axis1]}, ii_d_v2{v2[retained_axis0], v2[retained_axis1]};
            if(point_inside_triangle(ii_d_v0, ii_d_v1, ii_d_v2, ii_critical_point))
            {
                std::map<SymEngine::RCP<const SymEngine::Basic>, SymEngine::RCP<const SymEngine::Basic>, SymEngine::RCPBasicKeyLess> substituter{
                    {retained_symbol0, SymEngine::real_double(ii_critical_point.x())},
                    {retained_symbol0, SymEngine::real_double(ii_critical_point.y())},
                };
                auto substitute_result = ii_d_ex.subs(substituter);
                auto substitute_result_double = dynamic_cast<const SymEngine::RealDouble*>(substitute_result.get_basic().get());
                update_min_max(ii_d_ex, substituter, desired_min, desired_max);
            }
            // calculate min max in line
            std::vector<double> min_array, max_array;
            std::vector<std::vector<Eigen::Vector2d>> lines = {
                {ii_d_v0, ii_d_v1},
                {ii_d_v0, ii_d_v2},
                {ii_d_v1, ii_d_v2},
            };
            for(const auto& line: lines)
            {
                auto min_max_value = min_max_in_line(line[0], line[1], {retained_symbol0, retained_symbol1}, ii_d_ex);
                min_array.emplace_back(min_max_value.first);
                max_array.emplace_back(min_max_value.second);
            }
            desired_min = std::min(desired_min, *std::min_element(min_array.begin(), min_array.end()));
            desired_max = std::min(desired_max, *std::max_element(max_array.begin(), max_array.end()));
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
}

// write code from beginning to end
void min_max_multivariable_polynomial()
{
    // TODO: this is fake points and normals, should be replace
    // interpolated points and its normals
    std::vector<Eigen::Vector3d> interpolated_points;
    std::vector<Eigen::Vector3d> interpolated_normals;

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
        ex = ex + pow(interpolated_normal.x() * (x_ - interpolated_point.x()) + interpolated_normal.y() * (y_ - interpolated_point.y()) + interpolated_normal.z() * (z_ - interpolated_point.z()), 2);
    }
    std::cout << ex << std::endl;

    // calculate deriative of x, y, z
    SymEngine::Expression derivative_x = ex.diff(x);
    std::cout << "partial derivative of x: " << derivative_x << std::endl;
    SymEngine::Expression derivative_y = ex.diff(y);
    std::cout << "partial derivative of y: " << derivative_y << std::endl;
    SymEngine::Expression derivative_z = ex.diff(z);
    std::cout << "partial derivative of z: " << derivative_z << std::endl;

    // calculate 3d extream point
    auto solutions_3d = SymEngine::linsolve({ derivative_x.get_basic(), derivative_y.get_basic(), derivative_z.get_basic() }, { x, y, z });
    Eigen::Vector3d extream_3d_point;
    auto get_3d_coor_func = [&](int index) {
        auto solution_3d_coor = solutions_3d[index];
        auto number = dynamic_cast<const SymEngine::RealDouble*>(solution_3d_coor.get());
        if(number == nullptr)
        {
            throw std::exception("solution 3d cannot cast to double");
        }
        else
        {
            extream_3d_point[index] = number->as_double();
        }
    };
    get_3d_coor_func(0);
    get_3d_coor_func(1);
    get_3d_coor_func(2);

    // TODO: this is fake 3d mesh, should be replace in the future
    // check if 3d extream point is in 3d range
    std::vector<Eigen::Vector3d> vertices;
    std::vector<Eigen::Vector3i> triangles;
    double desired_max = std::numeric_limits<double>::min(), desired_min = std::numeric_limits<double>::max();
    if(point_inside_mesh(vertices, triangles, extream_3d_point))
    {
        // get extream value and substitute max and min
        std::map<SymEngine::RCP<const SymEngine::Basic>, SymEngine::RCP<const SymEngine::Basic>, SymEngine::RCPBasicKeyLess> substituter{
            {x, SymEngine::real_double(extream_3d_point.x())},
            {y, SymEngine::real_double(extream_3d_point.y())},
            {z, SymEngine::real_double(extream_3d_point.z())},
        };
        auto substitute_result = ex.subs(substituter);
        auto substitute_result_double = dynamic_cast<const SymEngine::RealDouble*>(substitute_result.get_basic().get());
        if(substitute_result_double)
        {
            auto value = substitute_result_double->as_double();
            if(value > desired_max)
            {
                desired_max = value;
            }
            if(value < desired_min)
            {
                desired_min = value;
            }
        }
        else
        {
            throw std::exception("cannot subtitute expression with double");
        }
    }

    // get max/min value from 3d boundary
    // for triangle in triangles
    // compute min/max value from each triangle
    std::vector<double> min_vec, max_vec;
    for(const auto& triangle: triangles)
    {
        auto min_max_value = min_max_in_triangle(vertices[triangle[0]], vertices[triangle[1]], vertices[triangle[2]], {x, y, z}, ex);
        min_vec.emplace_back(min_max_value.first);
        max_vec.emplace_back(min_max_value.second);
    }
    desired_min = std::min(desired_min, *std::min_element(min_vec.begin(), min_vec.end()));
    desired_max = std::min(desired_max, *std::max_element(min_vec.begin(), min_vec.end()));
}