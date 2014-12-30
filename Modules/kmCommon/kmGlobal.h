#ifndef KM_GLOBAL_H
#define KM_GLOBAL_H

#define KM_DEBUG_ERROR(X) std::cout<<"[DEBUG ERROR]"<<(X)<<std::endl;
#define KM_DEBUG_INFO(X) std::cout<<"[DEBUG INFO] "<<(X)<<std::endl;
#define KM_DEBUG_PRINT(X,Y) std::cout<<"[DEBUG PRINT] "<<X<<":"<<(Y)<<std::endl;
#define KM_DEBUG_PRINT_VALUE(X) std::cout<<"[DEBUG PRINT] "<<"X:"<<(X)<<std::endl;
#define KM_PRINT_EXCEPTION(E) std::cout<<"[EXCEPTION] "<<(E)<<std::endl;
#define KM_ASSERT(X) if(!(X)) std::cout<<"[ASSERT FAIL] "<<std::endl; 

#include "itkVector.h"
#include "itkPoint.h"
#include <vector>
#include <map>
#include <fstream>
#include <string>

static double SIGMA = 1.0; //�����ݶ�profileʱ�������ݶȳ������Ĳ���sigma
static int PROFILE_DIM = 9;
static double PROFILE_SPACING = 1.5; //���� PROFILE_LENTH/(PROFILE_DIM-1)
static double SHIFT_INSIDE = 3.0; //����������profileʱ�������ӱ������������ƫ�Ƶľ���
static double SHIFT_OUTSIDE = 3.0; //����������profileʱ�������ӱ������������ƫ�Ƶľ���
static double SHIFT_BOUNDARY = 1.0;
static int NUMBER_OF_INSIDE_PER_POINT = 3; // ��ÿ��MESH�㸽������inside profile��������
static int NUMBER_OF_BOUNDARY_PER_POINT = 6; // ��ÿ��MESH�㸽������boundary profile��������
static int NUMBER_OF_OUTSIDE_PER_POINT = 9; // ��ÿ��MESH�㸽������outside profile��������

#define RESAMPLE_SPACING 2.0

#ifndef B_CLASS_TYPE
typedef int BClassType;
#define IPClass 0
#define BPClass 1
#define OPClass 2
#endif

enum PROCESS_PHASE
{
	INITIALIZATION = 0,
	ROI_LOCATION,
	ADABOOST_REGION_SEGMENTATION,
	MODEL_FITTING_BY_DISTMAP,
	TUMOR_DETECTION,
	MODEL_FITTING_BY_LIVER_PROFILE,
	DEFORMATION_BY_LIVER_PROFILE,
	DEFORMATION_BY_BOUNDARY_PROFILE,
	TRAINING
};

enum PROFILE_CATEGORY
{
	BOUNDARY = 0,
	PLAIN,
	LIVER,
	COORDINATE,
	DEFAULT
};

namespace km
{
	class ProfileCategoryUtils
	{
	public:
		static PROFILE_CATEGORY string2category(const std::string str)
		{
			if (str.compare("BOUNDARY") == 0){
				return BOUNDARY;
			}else if (str.compare("PLAIN") == 0){
				return PLAIN;
			}else if (str.compare("LIVER") == 0){
				return LIVER;
			}else if (str.compare("COORDINATE") == 0){
				return COORDINATE;
			}else{
				return DEFAULT;
			}
		}

		static std::string category2string(PROFILE_CATEGORY category)
		{
			if (category == BOUNDARY){
				return "BOUNDARY";
			}else if (category == PLAIN){
				return "PLAIN";
			}else if (category == LIVER){
				return "LIVER";
			}else if (category == COORDINATE){
				return "COORDINATE";
			}else{
				return "DEFAUTL";
			}
		}
	};

	static PROCESS_PHASE g_phase;
	static std::vector<std::pair<double, double>> g_liverThresholds;
	static itk::Point<double> g_liverCentroid;
	static std::map<int, double> g_varianceMap;
	static std::map<int, double> g_fittingErrorMap;
	static double g_shape_penalty = 0.0;
	static double g_fitting_error_threshold = 0.6;
	static bool   g_disable_abnormal = true;
	static char   g_output_dir[1024];
	static int    g_number_clusters = 5;
	static double g_cluster_min_dist = 5.0;

	class Config
	{
	public:
		static void loadConfig(const char* filename)
		{
			std::string line;
			std::ifstream myfile (filename);

			if (myfile.is_open())
			{
				try{
					while ( getline (myfile,line) ){
						if (line == "#shape_penalty"){
							getline (myfile,line);
							g_shape_penalty = atof( line.c_str() );
						}else if (line == "#fitting_error_threshold"){
							getline (myfile,line);
							g_fitting_error_threshold = atof( line.c_str() );
						}else if (line == "#disable_abnormal"){
							getline (myfile,line);
							double val = atof(line.c_str());
							g_disable_abnormal = val>0?true:false;
						}else if (line == "#number_clusters"){
							getline (myfile,line);
							g_number_clusters = atoi( line.c_str() );
						}else if (line == "#cluster_min_dist"){
							getline (myfile,line);
							g_cluster_min_dist = atof( line.c_str() );
						}
					}
					myfile.close();
				}catch (...){
					std::cerr<<"Exception thrown when read from configuration file"<<std::endl;
				}
			}else{
				std::cerr<<"Unable to open configuration file: "<<filename<<std::endl;
			}

			std::cout<<"shape_penalty: "<<g_shape_penalty<<std::endl;
			std::cout<<"fitting_error_threshold: "<<g_fitting_error_threshold<<std::endl;
			std::cout<<"diable_abnormal: "<<g_disable_abnormal<<std::endl;
			std::cout<<"number_clusters: "<<g_number_clusters<<std::endl;
		}
	};
}

#endif