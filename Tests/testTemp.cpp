#include <iostream>
#include <string>

#include "testLib.h"

#define TestMacro(X) std::cout<<X<<std::endl;

void test(const std::string inputstr)
{
	//std::cout<<inputstr.<<std::endl;
}

int main(int argc, char* argv[])
{
	func1();
	func2();

	Container c;
	c.print();

	func3<double>(3.0);

	TestMacro(3.14);

	test(NULL);

	system("pause");

	return 0;
}