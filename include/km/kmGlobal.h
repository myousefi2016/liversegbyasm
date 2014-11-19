#ifndef KM_GLOBAL_H
#define KM_GLOBAL_H

#include "itkVector.h"
#include "itkPoint.h"
#include <vector>
#include <map>

#define KM_DEBUG_ERROR(X) std::cout<<"[DEBUG ERROR]"<<(X)<<std::endl;

static double SIGMA = 2.0; //采样梯度profile时，计算梯度场采样的参数sigma
static int NEIGHBOR_RADIUS = 1;
static int PROFILE_DIM = 7;
static double PROFILE_SPACING = 2.0; //毫米 PROFILE_LENTH/(PROFILE_DIM-1)
static double PROFILE_INSIDE_RATIO = 0.8;
static double SHIFT_INSIDE = 2.0; //采样表面内profile时，样本从表面上向表面内偏移的距离
static double SHIFT_OUTSIDE = 2.0; //采样表面外profile时，样本从表面上向表面外偏移的距离
static int NUMBER_OF_INSIDE_PER_POINT = 3; // 在每个MESH点附近采样inside profile的数量。
static int NUMBER_OF_BOUNDARY_PER_POINT = 1; // 在每个MESH点附近采样boundary profile的数量。
static int NUMBER_OF_OUTSIDE_PER_POINT = 3; // 在每个MESH点附近采样outside profile的数量。
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

	static string category2string(PROFILE_CATEGORY category);
	static void loadConfig(const char* filename);
}

#include "kmGlobal.hxx"

#endif