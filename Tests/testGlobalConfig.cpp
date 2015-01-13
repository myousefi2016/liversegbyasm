#include <iostream>
#include <string>
#include <map>

#include "kmCommon.h"

int main(int argc, char* argv[])
{
	const char* configFile = "D:\\Workspace\\LiverSegByASM\\liversegbyasm-v2\\Data\\config.txt";
	
	KM_DEBUG_INFO("Load configuration file...");
	km::Config::loadConfig(configFile);
	
	

	system("pause");

	return 0;
}