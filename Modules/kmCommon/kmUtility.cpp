#ifndef __KM_UTILITY_CPP
#define __KM_UTILITY_CPP

#include <sstream>
#include <string>
#include <fstream>
#include <vector>

namespace km
{
	/************************************************************************/
	/* Get Data List                                                        */
	/************************************************************************/
	int getDataList( std::string dataListFileName, std::vector<std::string> &files )
	{
		std::ifstream  fin(dataListFileName.c_str(), std::ios::in);  
		char  filename[1024]={0};
		int num = 0;
		while(fin.getline(filename, sizeof(filename)))  {
			std::string str_filename( filename );
			if(str_filename.size()==0){
				continue;
			}
			if( str_filename.find('#') != std::string::npos ){
				continue;
			}
			if ( str_filename.find( "END" ) != std::string::npos ){
				break;
			}
			files.push_back(str_filename);
			num ++;
		}  
		fin.clear();  
		fin.close();
		return num;
	}
}

#endif
