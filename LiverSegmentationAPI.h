#ifndef _LiverSegmentation_H
#define _LiverSegmentation_H

#include <fstream>						// file I/O
#include <map>

#include <itkImage.h>
#include <itkMesh.h>

#include "itkTimeProbesCollectorBase.h"
#include "itkMemoryProbesCollectorBase.h"

#include "itkDefaultDynamicMeshTraits.h"
#include "itkSimplexMesh.h"
#include "Representers/ITK/itkSimplexMeshRepresenter.h"
#include "statismo_ITK/itkStatisticalModel.h"
#include "statismo_ITK/itkStatisticalShapeModelTransform.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkVersorRigid3DTransform.h"
#include <itkEuler3DTransform.h>
#include <itkCenteredEuler3DTransform.h>
#include <itkCompositeTransform.h>
#include <itkAffineTransform.h>
#include <itkCenteredAffineTransform.h>
#include <itkScaleTransform.h>
#include "itkRigid3DTransform.h"
#include "itkSimilarity3DTransform.h"
#include <itkFixedCenterOfRotationAffineTransform.h>
#include "itkBSplineTransform.h"
#include <itkScaleVersor3DTransform.h>
#include <itkScaleSkewVersor3DTransform.h>
#include "itkSimilarity3DTransform.h"

#include "itkSimpleFilterWatcher.h"
#include "itkRelabelComponentImageFilter.h"

#include "itkLinearInterpolateImageFunction.h"
#include "itkSimpleFilterWatcher.h"

#include "itkListSample.h"
#include "itkKdTreeGenerator.h"
#include "itkWeightedCentroidKdTreeGenerator.h"

#include "itkDeformableSimplexMesh3DAdaboostClassifierForceFilter.h"

#include "LiverSegmentationAPI_Export.h"
#include "kmGlobal.h"
#include "AdaSegmentAPI.h"
#include "kmNotifierBase.h"
#include "kmProfileClassifier.h"
#include "kmUtility.h"
#include "kmVtkItkUtility.h"
#include "kmProcessing.h"
#include "kmModelFitting.h"

extern char* outputdir;

const unsigned int Dimension = 3;
typedef itk::Image<float, Dimension> ShortImageType;
typedef itk::Image<unsigned char, Dimension> UCharImageType;
typedef itk::Image<float, Dimension> FloatImageType;
typedef itk::CovariantVector<double, Dimension> GradientVectorType;
typedef itk::Image<GradientVectorType, Dimension> GradientImageType;

typedef itk::LinearInterpolateImageFunction<ShortImageType, double> ShortInterpolatorType;
typedef itk::LinearInterpolateImageFunction<FloatImageType, double> FloatInterpolatorType;
typedef itk::LinearInterpolateImageFunction<GradientImageType, double> GradientInterpolatorType;

typedef ShortImageType::PixelType PixelType;         
typedef ShortImageType::IndexType IndexType;         
typedef ShortImageType::PointType PointType;         
typedef ShortImageType::SizeType  SizeType;          
typedef ShortImageType::SpacingType SpacingType;     
typedef ShortImageType::RegionType RegionType;       
typedef ShortImageType::DirectionType DirectionType;

typedef double MeshPixelType;
typedef itk::SimplexMeshRepresenter<MeshPixelType, Dimension> RepresenterType;
typedef itk::StatisticalModel<RepresenterType>                StatisticalModelType;
typedef RepresenterType::MeshType                             MeshType;

typedef itk::Mesh<MeshPixelType, 3> TriangleMeshType;

typedef km::NotifierBase<UCharImageType, MeshType> NotifierType;

bool WRITE_MIDDLE_RESULT = true;

char* outputdir = NULL;

using namespace km;

LiverSegmentationAPI_EXPORT 
void LiverSeg( const NotifierType* notifier,
			  const char* inputImageFile,
			  const char* SSMFile,
			  const char* boundaryClassifierFile,
			  const char* liverClassifierFile,
			  const char* adaboostSegmentFile,
			  const char* geoFile,
			  const char* atlasImageFile,
			  const char* configFile,
			  const char* varianceMap);
	
#endif //_LiverSegmentation_H