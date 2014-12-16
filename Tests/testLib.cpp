#ifndef __testLib_cpp
#define __testLib_cpp

#include "testLib.h"

void func2()
{
	ImageType::Pointer img = ImageType::New();
	img->Print(std::cout);
	std::cout<<"I'm func2. I'm in testLib.cpp"<<std::endl;
}

#endif