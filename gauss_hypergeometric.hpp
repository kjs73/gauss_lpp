#ifndef GAUSS_HYPERGEOMETRIC_HPP
#define GAUSS_HYPERGEOMETRIC_HPP

/*2F1(a,b,c,x)*/

#include <cmath>

class gauss_hypergeometric{
public:
	typedef long double float_t;
	typedef int int_t;
private:
	const float_t a;
	const float_t b;
	const float_t c;
	const float_t precision;
	const float_t inv_beta_val;
public:
	gauss_hypergeometric(const float_t a_, const float_t b_, const float_t c_)
	    : a(a_),
          b(b_),
          c(c_),
          precision(std::pow(float_t(10), -4 * static_cast<int>(std::numeric_limits<float_t>::digits10))),
          inv_beta_val(std::tgamma(c_) / std::tgamma(b_) / std::tgamma(c_ - b_))
	{
		if (c_ == 0) throw std::runtime_error("2F1: illegal input for parameter c");
		if (b_ == 0) throw std::runtime_error("2F1: illegal input for parameter c");
		if (c_ == b_) throw std::runtime_error("2F1: illegal parameter combination input b==c");
		if (a_ == c_) throw std::runtime_error("2F1: illegal parameter combination input a==c");
	}
	float_t operator()(const float_t x) const
    {
		return inv_beta_val * c_factor(x);
	}
	float_t c_factor(const float_t x) const
    {
		float_t t(0);
		float_t t_(0);
		float_t sum(integrand(t, x));
		float_t h(float_t(1) / float_t(128));
		int_t n(1);
		float_t term;
		float_t term_;
		do {
			t += h;
			t_-= h;
			term = integrand(t, x);
			term_= integrand(t_, x);
			sum += term;
			sum += term_;
			++n;
		} while(std::abs(term) / std::abs(sum) > precision || std::abs(term_) / std::abs(sum) > precision );
		return h * sum;
	}
	float_t integrand(const float_t t, const float_t x) const
    {
		const float_t tmp = alpha(t);
		return std::abs(alpha_(t)) * std::pow(tmp, b-1) * std::pow(1 - tmp, c - b - 1) * std::pow(1 - x * tmp, -a);
	}
	float_t alpha(const float_t t) const
    {
		return float_t(0.5) * (float_t(1) + std::tanh(std::sinh(t)));
	}
	float_t alpha_(const float_t t) const
    {
		const float_t tmp = float_t(1) / std::cosh(std::sinh(t));
		return float_t(0.5) * std::cosh(t) * tmp * tmp;
	}
};

#endif //GAUSS_HYPERGEOMETRIC_HPP
