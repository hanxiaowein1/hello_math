/**
 * @file charles_symengine_common.h
 * @author your name (you@domain.com)
 * @brief extend some function of symengine for self use
 * @version 0.1
 * @date 2024-12-10
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef __CHARLES_SYMENGINE_COMMON_H__
#define __CHARLES_SYMENGINE_COMMON_H__

#include <exception>

#include <symengine/expression.h>
#include <symengine/solve.h>
#include <symengine/basic.h>
#include <symengine/add.h>
#include <symengine/matrix.h>
#include <symengine/integer.h>
#include <Eigen/Dense>

namespace hello_symengine
{
extern double deviation;
};

double get_double_from_solution(const SymEngine::RCP<const SymEngine::Basic>& solution);

double substitute_with_number(
    const SymEngine::Expression& ex,
    const std::map<SymEngine::RCP<const SymEngine::Basic>, SymEngine::RCP<const SymEngine::Basic>, SymEngine::RCPBasicKeyLess>& double_substituter
);

double get_coefficient_of_symbol(const SymEngine::Expression& expression, const SymEngine::RCP<const SymEngine::Symbol>& symbol);

class NoUniqueSolutionException : public std::exception
{
public:
    const char* what() const noexcept override {
        return "no unique solution";
    }
};

bool in_deviation(const double& v1, const double& v2);
bool in_deviation(const Eigen::Vector2d& vec1, const Eigen::Vector2d& vec2);
bool in_deviation(const Eigen::Vector3d& vec1, const Eigen::Vector3d& vec2);

#endif