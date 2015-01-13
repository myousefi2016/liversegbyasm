#ifndef __KM_UTILITY_H
#define __KM_UTILITY_H

#include <sstream>
#include <string>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

namespace km
{
	const double pi = 3.14159265;

	class NormalDistribution
	{
	public:
		double mu;
		double sigma;

		NormalDistribution(const double _mu, const double _sigma) : mu(_mu), sigma(_sigma) {}
		NormalDistribution() : mu(0), sigma(1) {}

		void set_mu(const double _mu){this->mu = _mu;}
		void set_sigma(const double _sigma){this->sigma = _sigma;}

		double pdf(const double x)
		{
			return exp( -1 * (x - mu) * (x - mu) / (2 * sigma * sigma)) / (sigma * sqrt(2 * pi));
		}

		double erf(const double x)
		{
			double y = 1.0 / ( 1.0 + 0.3275911 * x);   
			return 1 - (((((
				+ 1.061405429  * y
				- 1.453152027) * y
				+ 1.421413741) * y
				- 0.284496736) * y 
				+ 0.254829592) * y) 
				* exp (-x * x);
		}

		static double cdf(double x)
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

		//double cdf(const double x)
		//{
		//	//// Integral from a to x;
		//	//const double ninf = mu - 10 * sigma; // Dependent on the type of distribution, tune as appropriate
		//	//double sum = 0;
		//	//double n = 1e2; // tune for speed/accuracy
		//	//double c = (x - ninf) / n;

		//	//for (double k = 1.; k < n-1; k++)
		//	//	sum += pdf( ninf + k*c);

		//	//return c * ((pdf(x) + pdf(ninf))/2 + sum);

		//	return 0.5 * (1 + erf((x - mu) / (sigma * sqrt(2.))));
		//}

		double cdf_inside(const double x)
		{
			return std::abs(cdf(x)-cdf(x+2.0*(mu-x)));
		}

		double cdf_outside(const double x)
		{
			return std::abs(1.0 - cdf_inside(x));
			//return std::min(0.0, 1.0 - cdf_inside(x));
		}

		void print()
		{
			std::cout<<"[Normal Distribution] mu: "<<mu<<", sigma: "<<sigma<<std::endl;
		}
	};

	class Math
	{
	public:
		static double setToBetween(double val, double lowerBound, double upperBound)
		{
			if(val<lowerBound)
				return lowerBound;
			else if(val>upperBound)
				return upperBound;
			else
				return val;
		}
		
		static void calculateMeanAndSigma(const std::vector<double> & list, double & mean, double & sigma)
		{
			mean = sigma = 0.0;
			for(int i=0;i<list.size();i++)
			{
				mean += list[i];
			}
			mean /= list.size();
			for(int i=0;i<list.size();i++)
			{
				sigma += std::pow(list[i]-mean, 2);
			}
			sigma = std::sqrt( sigma/list.size() );
		}
	};
	
	/************************************************************************/
	/* Get Data List                                                        */
	/************************************************************************/
	int getDataList( std::string dataListFileName, std::vector<std::string> &files );
}

#undef kmStaticImageMacro

#endif
