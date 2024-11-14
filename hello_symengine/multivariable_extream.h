#ifndef __CHARLES_MULTIVARIABLE_EXTREAM_H__
#define __CHARLES_MULTIVARIABLE_EXTREAM_H__

#include <tuple>
#include <map>
#include <symengine/expression.h>
#include <symengine/solve.h>
#include <symengine/basic.h>
#include <symengine/add.h>
#include <symengine/matrix.h>
#include <symengine/integer.h>
#include <Eigen/Dense>

std::tuple<double, double, std::map<SymEngine::RCP<const SymEngine::Symbol>, double, SymEngine::RCPBasicKeyLess>, std::map<SymEngine::RCP<const SymEngine::Symbol>, double, SymEngine::RCPBasicKeyLess>> min_max_3_variable_quadratic_polynomial(
    const SymEngine::Expression& equation,
    const std::vector<SymEngine::RCP<const SymEngine::Symbol>>& symbol_axes,
    const std::vector<Eigen::Vector3d>& vertices,
    const std::vector<Eigen::Vector3i>& triangles
);

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

std::ostream& operator<<(std::ostream& os, const ExtreamLocation& m);

std::ostream& operator<<(std::ostream& os, const std::map<SymEngine::RCP<const SymEngine::Symbol>, double, SymEngine::RCPBasicKeyLess>& m);

#endif