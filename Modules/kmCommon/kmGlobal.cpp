#ifndef KM_GLOBAL_CPP
#define KM_GLOBAL_CPP

#include "kmGlobal.h"

namespace km
{
	void loadConfig(const char* filename)
	{
		std::string line;
		std::ifstream myfile (filename);

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
					else if (line == "#fitting_error_threshold")
					{
						getline (myfile,line);
						g_fitting_error_threshold = atof( line.c_str() );
					}
					else if (line == "#disable_abnormal")
					{
						getline (myfile,line);
						double val = atof(line.c_str());
						g_disable_abnormal = val>0?true:false;
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

		std::cout<<"shape_penalty: "<<g_shape_penalty<<std::endl;
		std::cout<<"fitting_error_threshold: "<<g_fitting_error_threshold<<std::endl;
		std::cout<<"diable_abnormal: "<<g_disable_abnormal<<std::endl;
	}
}

#endif