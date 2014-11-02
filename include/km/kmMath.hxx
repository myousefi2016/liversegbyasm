#ifndef KM_MATH_HXX
#define KM_MATH_HXX

#include "kmMath.h"

namespace km
{
	double cdf(double x)
	{
		// constants
		double a1 =  0.254829592;
		double a2 = -0.284496736;
		double a3 =  1.421413741;
		double a4 = -1.453152027;
		double a5 =  1.061405429;
		double p  =  0.3275911;

		// Save the sign of x
		int sign = 1;
		if (x < 0)
			sign = -1;
		x = fabs(x)/sqrt(2.0);

		// A&S formula 7.1.26
		double t = 1.0/(1.0 + p*x);
		double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

		return 0.5*(1.0 + sign*y);
	}

	double cdf_inside(double x)
	{
		if (x == 0)
		{
			return 0;
		}

		return std::abs(cdf(x) - cdf(-x)); 
	}

	double cdf_outside(double x)
	{
		return 1.0 - cdf_inside(x);
	}
}

#endif