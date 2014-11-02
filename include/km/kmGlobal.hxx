#ifndef KM_GLOBAL_HXX
#define KM_GLOBAL_HXX

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include "kmGlobal.h"

namespace km
{
	static void loadConfig(const char* filename)
	{
		string line;
		ifstream myfile (filename);

		if (myfile.is_open())
		{
			try
			{
				while ( getline (myfile,line) )
				{
					if (line == "#shape_penalty")
					{
						getline (myfile,line);
						g_shape_penalty = atof( line.c_str() );
					}
				}

				myfile.close();
			}
			catch (...)
			{
				std::cerr<<"Exception thrown when read from configuration file"<<std::endl;
			}
		}
		else
		{
			std::cerr<<"Unable to open configuration file: "<<filename<<std::endl;
		}

		std::cout<<"g_shape_penalty: "<<g_shape_penalty<<std::endl;
	}

}

#endif