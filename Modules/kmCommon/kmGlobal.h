#ifndef KM_GLOBAL_H
#define KM_GLOBAL_H

#define KM_DEBUG_ERROR(X) std::cout<<"[DEBUG ERROR]"<<(X)<<std::endl;
#define KM_DEBUG_INFO(X) std::cout<<"[DEBUG INFO] "<<(X)<<std::endl;
#define KM_DEBUG_PRINT(X,Y) std::cout<<"[DEBUG PRINT] "<<X<<":"<<(Y)<<std::endl;
#define KM_DEBUG_PRINT_VALUE(X) std::cout<<"[DEBUG PRINT] "<<"X:"<<(X)<<std::endl;
#define KM_PRINT_EXCEPTION(E) std::cout<<"[EXCEPTION] "<<(E)<<std::endl;
#define KM_ASSERT(X) if(!(X)) std::cout<<"[ASSERT FAIL] "<<std::endl;

#define DEFINE_STATIC_GLOBAL_VARIABLE(varName, varType, varDefaultVal) \
	static varType g_##varName = varDefaultVal;

//Only support double, int, bool
#define ASSIGN_GLOBAL_VARIABLE(varName, varType, varVal) \
	g_##varName = (varType)varVal;

#define READ_CONFIG_LINE(varName, varType, file, line) \
	if (line == #varName) \
{ \
	getline (file,line); \
	ASSIGN_GLOBAL_VARIABLE(varName, varType, atof(line.c_str()) ); \
}

#include "itkVector.h"
#include "itkPoint.h"
#include <vector>
#include <map>
#include <fstream>
#include <string>

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
			}else{
				return "DEFAUTL";
			}
		}
	};

	static char   g_output_dir[1024];
	static PROCESS_PHASE g_phase;

	static double g_resample_spacing = 2.0; //Spacing used to re-sample input image.
	static double g_sigma = 1.0; //Sigma used for calculate gradient
	static int    g_profile_dim = 9; //Number of sample points on each profile.
	static double g_profile_spacing = 1.0; //Distance between each sample point on each profile.
	static double g_shift_inside = 3.0; //Distance shifted toward inside along normal direction during training.
	static double g_shift_outside = 3.0; //Distance shifted toward outside along normal direction during training.
	static double g_shift_boundary = 1.0; //Random distance shifted near true boundary during training.
	static int    g_number_of_inside_per_point = 3; //Number of inside sample points for each landmark.
	static int    g_number_of_boundary_per_point = 6; //Number of true boundary points for each landmark.
	static int    g_number_of_outside_per_point = 9; //Number of outside sample points for each landmark.

	static double g_shape_penalty = 0.0;
	static int    g_number_clusters = 5;
	static double g_cluster_min_dist = 5.0;
	static int    g_number_principle_components = 17;
	static double g_alpha = 0.3;
	static double g_beta  = 0.3;
	static double g_kappa = 0.3;
	static double g_gamma = 0.1;
	static int    g_rigidity = 1;
	static bool   g_landmark_status_evalution = false;

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
						}else if (line == "#number_clusters"){
							getline (myfile,line);
							g_number_clusters = atoi( line.c_str() );
						}else if (line == "#cluster_min_dist"){
							getline (myfile,line);
							g_cluster_min_dist = atof( line.c_str() );
						}else if (line == "#number_principle_components"){
							getline (myfile,line);
							g_number_principle_components = atoi( line.c_str() );
						}else if (line == "#alpha"){
							getline (myfile,line);
							g_alpha = atof( line.c_str() );
						}else if (line == "#beta"){
							getline (myfile,line);
							g_beta = atof( line.c_str() );
						}else if (line == "#kappa"){
							getline (myfile,line);
							g_kappa = atof( line.c_str() );
						}else if (line == "#gamma"){
							getline (myfile,line);
							g_gamma = atof( line.c_str() );
						}else if (line == "#rigidity"){
							getline (myfile,line);
							g_rigidity = atoi( line.c_str() );
						}else if (line == "#resample_spacing"){
							getline (myfile,line);
							g_resample_spacing = atof( line.c_str() );
						}else if (line == "#sigma"){
							getline (myfile,line);
							g_sigma = atof( line.c_str() );
						}else if (line == "#profile_dim"){
							getline (myfile,line);
							g_profile_dim = atoi( line.c_str() );
						}else if (line == "#profile_spacing"){
							getline (myfile,line);
							g_profile_spacing = atof( line.c_str() );
						}else if (line == "#shift_inside"){
							getline (myfile,line);
							g_shift_inside = atof( line.c_str() );
						}else if (line == "#shift_outside"){
							getline (myfile,line);
							g_shift_outside = atof( line.c_str() );
						}else if (line == "#shift_boundary"){
							getline (myfile,line);
							g_shift_boundary = atof( line.c_str() );
						}else if (line == "#number_of_inside_per_point"){
							getline (myfile,line);
							g_number_of_inside_per_point = atoi( line.c_str() );
						}else if (line == "#number_of_boundary_per_point"){
							getline (myfile,line);
							g_number_of_boundary_per_point = atoi( line.c_str() );
						}else if (line == "#number_of_outside_per_point"){
							getline (myfile,line);
							g_number_of_outside_per_point = atoi( line.c_str() );
						}
					}
					myfile.close();
				}catch (...){
					std::cerr<<"Exception thrown when read from configuration file"<<std::endl;
				}
			}else{
				std::cerr<<"Unable to open configuration file: "<<filename<<std::endl;
			}

			std::cout<<"[Global Config] shape_penalty: "<<g_shape_penalty<<std::endl;
			std::cout<<"[Global Config] number_clusters: "<<g_number_clusters<<std::endl;
			std::cout<<"[Global Config] landmark_status_evalution: "<<g_landmark_status_evalution<<std::endl;
			std::cout<<"[Global Config] alpha: "<<g_alpha<<std::endl;
			std::cout<<"[Global Config] beta: "<<g_beta<<std::endl;
			std::cout<<"[Global Config] kappa: "<<g_kappa<<std::endl;
			std::cout<<"[Global Config] gamma: "<<g_gamma<<std::endl;
			std::cout<<"[Global Config] rigidity: "<<g_rigidity<<std::endl;
			std::cout<<"[Global Config] number_clusters: "<<g_number_clusters<<std::endl;
			std::cout<<"[Global Config] cluster_min_dist: "<<g_cluster_min_dist<<std::endl;
			std::cout<<"[Global Config] number_principle_components: "<<g_number_principle_components<<std::endl;
			std::cout<<"[Global Config] resample_spacing: "<<g_resample_spacing<<std::endl;
			std::cout<<"[Global Config] sigma: "<<g_sigma<<std::endl;
			std::cout<<"[Global Config] profile_dim: "<<g_profile_dim<<std::endl;
			std::cout<<"[Global Config] profile_spacing: "<<g_profile_spacing<<std::endl;
			std::cout<<"[Global Config] shift_inside: "<<g_shift_inside<<std::endl;
			std::cout<<"[Global Config] shift_outside: "<<g_shift_outside<<std::endl;
			std::cout<<"[Global Config] shift_boundary: "<<g_shift_boundary<<std::endl;
			std::cout<<"[Global Config] number_of_inside_per_point: "<<g_number_of_inside_per_point<<std::endl;
			std::cout<<"[Global Config] number_of_boundary_per_point: "<<g_number_of_boundary_per_point<<std::endl;
			std::cout<<"[Global Config] number_of_outside_per_point: "<<g_number_of_outside_per_point<<std::endl;
		}
	};
}

#endif