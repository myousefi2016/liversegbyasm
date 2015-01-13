#include <iostream>
#include <string>

#include "itkImage.h"
#include "itkSimplexMesh.h"
#include "itkTimeProbesCollectorBase.h"
#include "itkMemoryProbesCollectorBase.h"
#include "itkDeformableSimplexMesh3DWithShapePriorFilter.h"
#include "itkAffineTransform.h"
#include "Representers/ITK/itkSimplexMeshRepresenter.h"
#include "statismo_ITK/itkStatisticalModel.h"
#include "statismo_ITK/itkStatisticalShapeModelTransform.h"

#include "kmCommon.h"

using namespace km;

int main(int argc, char* argv[])
{
	const char* inputImageFile = "D:\\Workspace\\ASM\\projects\\LiverSegbyASM\\experiments\\modelFitting\\output_20140815\\inputImageSmoothed.nii.gz";
	const char* geoImageFile = "D:\\Workspace\\ASM\\projects\\LiverSegbyASM\\experiments\\training\\alignment\\output_20140815\\geoImage.mha";
	const char* inputMeshFile = "D:\\Workspace\\ASM\\projects\\LiverSegbyASM\\experiments\\modelFitting\\output_20140815\\initializedMesh.vtk";
	const char* liverClassifierFile = "D:\\Workspace\\ASM\\projects\\LiverSegbyASM\\experiments\\training\\appearance\\output_20140815\\AdaboostClassifier_LIVER.h5";
	const char* boundaryClassifierFile = "D:\\Workspace\\ASM\\projects\\LiverSegbyASM\\experiments\\training\\appearance\\output_20140815\\AdaboostClassifier_BOUNDARY.h5";
	const char* ssmFile = "D:\\Workspace\\ASM\\projects\\LiverSegbyASM\\experiments\\buidModel\\shape\\output_20140815\\StatisticalShapeModel_20140815.h5";

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
	km::assigneMesh<MeshType>(inputMesh, 0.0);
	km::loadSimplexMeshGeometryData<MeshType>(geoImage, inputMesh);
	MeshType::PointType liverCentroid = km::getMeshCentroid<MeshType>(inputMesh); //Liver centroid
	
	typedef itk::AffineTransform<double, Dimension> RigidTransformType;
	RigidTransformType::Pointer rigidTransform = RigidTransformType::New();
	rigidTransform->SetCenter( liverCentroid );
	rigidTransform->SetIdentity(); //Rigid transform
		
	typedef itk::StatisticalShapeModelTransform<RepresenterType, double, Dimension> ShapeTransformType;
	ShapeTransformType::Pointer shapeTransform = ShapeTransformType::New();
	shapeTransform->SetStatisticalModel( model );
	shapeTransform->SetUsedNumberOfCoefficients( shapeTransform->GetNumberOfParameters() );
	shapeTransform->SetIdentity(); //Shape transform
	
	typedef km::ProfileClassifier ProfileClassifierType;
	ProfileClassifierType ProfileClassifier_Boundary;
	ProfileClassifier_Boundary.load(boundaryClassifierFile); //Boundary classifier
	ProfileClassifier_Boundary.print();
	ProfileClassifierType ProfileClassifier_Liver;
	ProfileClassifier_Liver.load(liverClassifierFile); //Liver classifier
	ProfileClassifier_Liver.print();

	itk::TimeProbesCollectorBase chronometer;
	itk::MemoryProbesCollectorBase memorymeter;
	memorymeter.Start( "segmentation" );
	chronometer.Start( "segmentation" );

	typedef itk::DeformableSimplexMesh3DWithShapePriorFilter<MeshType,
															 MeshType,
															 InputImageType,
															 StatisticalModelType,
															 RigidTransformType,
															 ShapeTransformType> DeformableFilterType;
	DeformableFilterType::Pointer deformaFilter = DeformableFilterType::New();
	deformaFilter->SetAlpha(0.0);
	deformaFilter->SetKappa(1.0);
	deformaFilter->SetIterations(15);
	deformaFilter->SetInput(inputMesh);
	deformaFilter->SetInputImage(inputImage);
	deformaFilter->SetStatisticalModel(model);
	deformaFilter->SetRigidTransform(rigidTransform);
	deformaFilter->SetShapeTransform(shapeTransform);
	deformaFilter->SetBoundaryClassifier(&ProfileClassifier_Boundary);
	deformaFilter->SetLiverClassifier(&ProfileClassifier_Liver);
	deformaFilter->Update();

	memorymeter.Stop( "segmentation" );
	chronometer.Stop( "segmentation" );
	chronometer.Report( std::cout );
	memorymeter.Report( std::cout );

	km::writeMesh<MeshType>("deformed.vtk", inputMesh);
	
	system("pause");
	return 0;
}