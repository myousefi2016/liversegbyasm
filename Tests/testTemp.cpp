#include <iostream>
#include <string>

#include "itkVector.h"
#include "itkImage.h"

int main(int argc, char* argv[])
{

	itk::Vector<double, 4> vec;
	vec[0] = 1;
	vec[1] = 20;
	vec[2] = 200;
	vec[3] = 300;

	vec[0] = 1.0/(vec[0]+1.0);
	vec[1] = 1.0/(vec[1]+1.0);
	vec[2] = 1.0/(vec[2]+1.0);
	vec[3] = 1.0/(vec[3]+1.0);

	double sum = vec[0]+vec[1]+vec[2]+vec[3];

	vec /= sum;

	std::cout<<vec<<std::endl;

	const char* tmp = NULL;
	char str[1024];
	sprintf(str, "ABCD%sEFG", tmp);

	std::cout<<str<<std::endl;

	typedef itk::Image<int, 3> ImageType;
	ImageType::Pointer img;
	std::cout<<img.IsNull()<<std::endl;

	img = NULL;
	std::cout<<img.IsNull()<<std::endl;

	system("pause");

	return 0;
}