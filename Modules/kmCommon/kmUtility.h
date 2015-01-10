#ifndef __KM_UTILITY_H
#define __KM_UTILITY_H

#include <sstream>
#include <string>
#include <fstream>

using namespace std;

namespace km
{
	class Math
	{
	public:
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

		static double cdf_inside(double x)
		{
			if (x == 0)
			{
				return 0;
			}

			return std::abs(cdf(x) - cdf(-x)); 
		}

		static double cdf_outside(double x)
		{
			return 1.0 - cdf_inside(x);
		}
		
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
