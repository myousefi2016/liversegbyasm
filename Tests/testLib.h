#ifndef __testLib_h
#define __testLib_h

#include <iostream>
#include "itkImage.h"

class Container
{
	public:
		void print()
		{
			std::cout<<"This is a container."<<std::endl;
		}
};

typedef itk::Image<float,3> ImageType;

inline void func1()
{
	ImageType::Pointer img = ImageType::New();
	img->Print(std::cout);
	std::cout<<"I'm function. I'm in testLib.h"<<std::endl;
}

template<class T>
void
	func3(typename T data)
{
	std::cout<<"I'm function 3. I'm in testLib.h"<<std::endl;
	std::cout<<"Data: "<<data<<std::endl;
}

void func2();

#endif