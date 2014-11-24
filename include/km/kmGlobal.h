#ifndef KM_GLOBAL_H
#define KM_GLOBAL_H

#include "itkVector.h"
#include "itkPoint.h"
#include <vector>
#include <map>

static double SIGMA = 2.0; //�����ݶ�profileʱ�������ݶȳ������Ĳ���sigma
static int NEIGHBOR_RADIUS = 1;
static int PROFILE_DIM = 9;
static double PROFILE_SPACING = 1.5; //���� PROFILE_LENTH/(PROFILE_DIM-1)
static double PROFILE_INSIDE_RATIO = 0.8;
static double SHIFT_INSIDE = 3.0; //����������profileʱ�������ӱ������������ƫ�Ƶľ���
static double SHIFT_OUTSIDE = 3.0; //����������profileʱ�������ӱ������������ƫ�Ƶľ���
static int NUMBER_OF_INSIDE_PER_POINT = 3; // ��ÿ��MESH�㸽������inside profile��������
static int NUMBER_OF_BOUNDARY_PER_POINT = 1; // ��ÿ��MESH�㸽������boundary profile��������
static int NUMBER_OF_OUTSIDE_PER_POINT = 5; // ��ÿ��MESH�㸽������outside profile��������
//#define NUMBER_OF_PROFILE_PER_POINT (NUMBER_OF_INSIDE_PER_POINT+NUMBER_OF_BOUNDARY_PER_POINT+NUMBER_OF_OUTSIDE_PER_POINT)

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
	static PROCESS_PHASE g_phase;
	static std::vector<std::pair<double, double>> g_liverThresholds;
	static itk::Point<double> g_liverCentroid;
	static double g_shape_penalty = 0.0;

	std::string category2string(PROFILE_CATEGORY category)
	{
		if (category == BOUNDARY)
		{
			return "BOUNDARY";
		}
		else if (category == PLAIN)
		{
			return "PLAIN";
		}
		else if (category == LIVER)
		{
			return "LIVER";
		}
		else if (category == COORDINATE)
		{
			return "COORDINATE";
		}
		else
		{
			return "DEFAUTL";
		}
	}

	static void loadConfig(const char* filename)
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