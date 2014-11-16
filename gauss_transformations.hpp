#ifndef GAUSS_TRANSFORMATIONS_HPP
#define GAUSS_TRANSFORMATIONS_HPP

/*2F1(a,b,c,x) with some transformations to make the integral more well-behaved*/

#include "gauss_hypergeometric.hpp"

class gauss_transformations{
public:
	typedef gauss_hypergeometric::float_t float_t;
private:
	const float_t a;
	const float_t b;
	const float_t c;
	const float_t EPS;
public:
	gauss_transformations(const float_t a_, const float_t b_, const float_t c_)
        : a(a_),
          b(b_),
          c(c_),
          EPS(1e-13)
    {}
	bool close_enough(const float_t x, const float_t y) const
    {
	    return std::fabs(x - y) < EPS;
	}
	float_t operator()(const float_t x) const
    {
		if (close_enough(a, c)) {
			return std::pow(1 - x, -b);
		}
		else if (close_enough(b, c)) {
			return std::pow(1 - x, -a);
		}
		else if (a >= 0) {
			const float_t top1 = a - 1 - c;
			const float_t top2 = c - 2 * (a - 1) + (a - 1 - b) * x;
			const float_t bottom = (a - 1) * (x - 1);
			if (close_enough(bottom, 0)) {
                throw std::runtime_error("2F1: would divide by zero: (a-1)*(x-1)==0");
            }
			return (top1 / bottom) * gauss_transformations(a - 2, b, c)(x) + (top2 / bottom) * gauss_transformations(a - 1, b, c)(x);
		}
		else if (b < 1) {
			const float_t top1 = (b + 1) * (x - 1);
			const float_t top2 = c - 2 * (b + 1) + (b + 1 - a) * x;
			const float_t bottom = b + 1 - c;
			if (close_enough(bottom, 0)) {
                throw std::runtime_error("2F1: would divide by zero: b+1-c==0");
            }
			return (top1 / bottom) * gauss_transformations(a, b + 2, c)(x) - (top2 / bottom) * gauss_transformations(a, b + 1, c)(x);
		}
		else if(c <= (b + 1)) {
			const float_t top1 = (a - c - 1) * (b - c - 1) * x;
			const float_t top2 = -c + (2 * c + 1 - a - b) * x;
			const float_t bottom1 = c * (c + 1) * (1 - x);
			const float_t bottom2 = c * (1 - x);
			if (close_enough(bottom1, 0)) {
                throw std::runtime_error("2F1: would divide by zero: c*(c+1)*(1-x)==0");
            }
			if (close_enough(bottom2, 0)) {
                throw std::runtime_error("2F1: would divide by zero: c*(1-x)==0");
            }
			return (top1 / bottom1) * gauss_transformations(a, b, c + 2)(x) - (top2 / bottom2) * gauss_transformations(a, b, c + 1)(x);
		}
		else {
			return gauss_hypergeometric(a, b, c)(x);
		}
	}		
};

#endif //GAUSS_TRANSFORMATIONS_HPP
