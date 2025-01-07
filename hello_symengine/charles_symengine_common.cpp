#include "charles_symengine_common.h"
#include <symengine/visitor.h>

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

double get_coefficient_of_symbol(const SymEngine::Expression& expression, const SymEngine::RCP<const SymEngine::Symbol>& symbol)
{
    using namespace SymEngine;
    auto ret = coeff(expression, *symbol, *integer(1));
    try
    {
        double coeff = get_double_from_solution(ret);
        return coeff;
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        return 0.0f;
    }
}
