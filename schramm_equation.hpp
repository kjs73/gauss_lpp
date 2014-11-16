#ifndef SCHRAMM_EQUATION_HPP
#define SCHRAMM_EQUATION_HPP

#include <algorithm>

#include "gauss_transformations.hpp"

class schramm_equation{
public:
    typedef long double float_t;
private:
    const float_t kappa;
    const float_t EPS;
    const bool space_filling;
    const float_t constant;
    gauss_transformations gh;
public:
    schramm_equation(const float_t kappa_)
        : kappa(kappa_),
          EPS(1e-13),
          space_filling(fabs(8 - kappa_) < EPS || kappa_ > 8),
          constant(space_filling ? 0 : tgamma(4 / kappa_) / std::sqrt(M_PI) / tgamma((8 - kappa_) / (2 * kappa_))),
          gh(gauss_transformations(float_t(1) / float_t(2), 4 / kappa_, float_t(3) / float_t(2)))
    {}
    float_t operator()(const float_t phi) const
    {
        if (phi < 0 || phi > M_PI) {
            throw std::runtime_error("schramm_equation::operator(): illegal input");
        }
        if (phi < EPS) {
            return 1;
        }
        if (space_filling) {
            return float_t(0.5);
        }
        const float_t tmp = cot(phi);
        const float_t result = float_t(1) / float_t(2) + constant * tmp * gh(-tmp * tmp);
        if (result < -EPS || result > 1 + EPS) {
            throw std::runtime_error("schramm_equation::operator() evaluation failed");
        }
        return std::min<float_t>(std::max<float_t>(result, 0), 1);
    }
    float_t cot(const float_t phi) const
    {
        return float_t(1) / std::tan(phi);
    }
};

#endif//SCHRAMM_EQUATION_HPP
