#include <iostream>
#include <string>
#include <map>

#include "vtkTimerLog.h"
#include "vtkNew.h"

#include "kmUtility.h"

int main(int argc, char* argv[])
{
	std::map<int, km::NormalDistribution> TempMap;

	km::NormalDistribution shapeDistribution(0,1);
	std::cout<<shapeDistribution.cdf(0)<<std::endl;
	std::cout<<shapeDistribution.cdf_inside(2.5)<<std::endl;
	//std::cout<<shapeDistribution.cdf_inside(0)<<std::endl;
	//std::cout<<shapeDistribution.cdf_outside(0)<<std::endl;
	//std::cout<<shapeDistribution.cdf_outside(3.0)<<std::endl;

	//std::cout<<"---------------------Insert-----------------------"<<std::endl;
	//int N = 10000000;
	//for (int i=0;i<N;i++)
	//{
	//	TempMap[i] = Gauss;
	//}

	//std::cout<<"---------------------Calculate-----------------------"<<std::endl;
	//vtkNew<vtkTimerLog> timer;
	//timer->StartTimer();
	//for (int i=0;i<N;i++)
	//{
	//	double pdfval = TempMap[i].cdf_outside(i);
	//	//std::cout<<pdfval<<std::endl;
	//}
	//timer->StopTimer();
	//std::cout<<"Calculate pdf done. "<<timer->GetElapsedTime() << " s." <<endl;


	//timer->StartTimer();
	//for (int i=0;i<10000;i++)
	//{
	//	double cdfoutsideval = TempMap[0].cdf_outside(30);
	//	//std::cout<<cdfoutsideval<<std::endl;
	//}
	//timer->StopTimer();
	//std::cout<<"Calculate cdf_outside done. "<<timer->GetElapsedTime() << " s." <<endl;

	system("pause");

	return 0;
}