#include <iostream>
#include <string>

#include "itkImage.h"
#include "itkSimplexMesh.h"
#include "itkDeformableSimplexMesh3DWithShapePriorFilter.h"
#include "itkAffineTransform.h"
#include "Representers/ITK/itkSimplexMeshRepresenter.h"
#include "statismo_ITK/itkStatisticalModel.h"
#include "statismo_ITK/itkStatisticalShapeModelTransform.h"

#include "kmGlobal.h"
#include "kmUtility.h"
#include "kmVtkItkUtility.h"
#include "kmProfileClassifier.h"

using namespace km;

int main(int argc, char* argv[])
{
	const char* inputImageFile = "";
	const char* geoImageFile = "";
	const char* inputMeshFile = "";
	const char* liverClassifierFile = "";
	const char* boundaryClassifierFile = "";
	const char* ssmFile = "";

	const int Dimension = 3;
	typedef itk::Image<unsigned short, Dimension> InputImageType;
	InputImageType::Pointer inputImage = km::readImage<InputImageType>( inputImageFile ); //Input image
	
	typedef itk::SimplexMeshRepresenter<double, Dimension> RepresenterType;
	typedef itk::StatisticalModel<RepresenterType>         StatisticalModelType;
	StatisticalModelType::Pointer model = StatisticalModelType::New();
	model->Load( ssmFile ); //Model
	
	GeometryImageType::Pointer geoImage = km::readImage<GeometryImageType>( geoImageFile );
	
	typedef RepresenterType::MeshType MeshType;
	MeshType::Pointer inputMesh = km::readMesh<MeshType>( inputMeshFile ); //Input mesh
	km::loadSimplexMeshGeometryData<MeshType>(geoImage, inputMesh);
	g_liverCentroid.CastFrom(km::getMeshCentroid<MeshType>(inputMesh)); //Liver centroid
	
	typedef itk::AffineTransform<double, Dimension> RigidTransformType;
	RigidTransformType::Pointer rigidTransform = RigidTransformType::New();
	rigidTransform->SetCenter( g_liverCentroid );
	rigidTransform->SetIdentity(); //Rigid transform
		
	typedef itk::StatisticalShapeModelTransform<RepresenterType, double, Dimension> ShapeTransformType;
	ShapeTransformType::Pointer shapeTransform = ShapeTransformType::New();
	shapeTransform->SetStatisticalModel( model );
	shapeTransform->SetUsedNumberOfCoefficients( shapeTransform->GetNumberOfParameters() );
	shapeTransform->SetIdentity(); //Shape transform
	
	typedef km::ProfileClassifier ProfileClassifierType;
	ProfileClassifierType ProfileClassifier_Boundary;
	ProfileClassifier_Boundary.load(boundaryClassifierFile); //Boundary classifier
	ProfileClassifierType ProfileClassifier_Liver;
	ProfileClassifier_Liver.load(liverClassifierFile); //Liver classifier

	typedef itk::DeformableSimplexMesh3DWithShapePriorFilter<MeshType,
															 MeshType,
															 InputImageType,
															 StatisticalModelType,
															 RigidTransformType,
															 ShapeTransformType> DeformableFilterType;
	DeformableFilterType::Pointer deformaFilter = DeformableFilterType::New();
	deformaFilter->SetAlpha(0.3);
	deformaFilter->SetKappa(0.1);
	deformaFilter->SetIterations(500);
	deformaFilter->SetInput(inputMesh);
	deformaFilter->SetInputImage(inputImage);
	deformaFilter->SetStatisticalModel(model);
	deformaFilter->SetRigidTransform(rigidTransform);
	deformaFilter->SetShapeTransform(shapeTransform);
	deformaFilter->SetBoundaryClassifier(&ProfileClassifier_Boundary);
	deformaFilter->SetLiverClassifier(&ProfileClassifier_Liver);
	deformaFilter->Update();
	
	system("pause");
	return 0;
}