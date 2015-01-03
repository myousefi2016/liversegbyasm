#include <iostream>
#include <string>

#include "itkDefaultDynamicMeshTraits.h"
#include "itkSimplexMesh.h"
#include "Representers/ITK/itkSimplexMeshRepresenter.h"
#include "statismo_ITK/itkStatisticalModel.h"
#include "statismo_ITK/itkStatisticalShapeModelTransform.h"
#include "itkAffineTransform.h"

#include "kmCommon.h"

int main(int argc, char* argv[])
{
	//const char* SSMFile = argv[1];
	const char* SSMFile = "D:\\Workspace\\ASM\\projects\\LiverSegbyASM\\experiments\\buidModel\\shape\\output_20140815\\StatisticalShapeModel_20140815.h5";
	const char* targetMeshFile = "D:\\Workspace\\ASM\\projects\\LiverSegbyASM\\experiments\\training\\shape\\output_20140815\\liverMesh.46.vtk";
	const char* configFile = "D:\\Workspace\\LiverSegByASM\\liversegbyasm-v2\\Data\\config.txt";

	KM_DEBUG_INFO("Load config file...");
	km::Config::loadConfig(configFile);

	const int Dimension = 3;
	typedef double MeshPixelType;
	typedef itk::SimplexMeshRepresenter<MeshPixelType, Dimension> RepresenterType;
	typedef itk::StatisticalModel<RepresenterType>                StatisticalModelType;
	typedef RepresenterType::MeshType                             MeshType;
	
	StatisticalModelType::Pointer model = StatisticalModelType::New();
	model->Load( SSMFile );

	MeshType::Pointer meanShape = model->DrawMean();
	MeshType::PointType centroid = km::getMeshCentroid<MeshType>(meanShape);

	typedef itk::AffineTransform<double, Dimension> RigidTransformType;
	typedef itk::StatisticalShapeModelTransform<RepresenterType, double, Dimension> ShapeTransformType;
	ShapeTransformType::Pointer shapeTransform = ShapeTransformType::New();
	shapeTransform->SetStatisticalModel(model);
	shapeTransform->SetIdentity();

	RigidTransformType::Pointer rigidTransform = RigidTransformType::New();
	rigidTransform->SetCenter(centroid);
	rigidTransform->SetIdentity();

	typedef km::SSMUtils<MeshType, StatisticalModelType, RigidTransformType, ShapeTransformType> SSMUtilsType;
	SSMUtilsType ssmUtils;
	ssmUtils.SetSSM(model);
	ssmUtils.SetRigidTransform(rigidTransform);
	ssmUtils.SetShapeTransform(shapeTransform);

	ssmUtils.cluster(km::g_number_clusters);

	MeshType::Pointer targetMesh = km::readMesh<MeshType>(targetMeshFile);
	MeshType::Pointer outputMesh = km::cloneMesh<MeshType, MeshType>(meanShape);

	km::writeMesh<MeshType>("targetMesh.vtk", targetMesh);
	km::writeMesh<MeshType>("unfittedMesh.vtk", meanShape);

	for (int i=0;i<1;i++)
	{
		std::cout<<"****************iter: "<<i<<"**************"<<std::endl;
		ssmUtils.update(targetMesh, outputMesh);
		char filename[1024];
		sprintf(filename, "fittedMesh-%d.vtk", i);
		km::writeMesh<MeshType>(filename, outputMesh);

		ssmUtils.printTransform();
	}

	//km::assigneMesh<MeshType>(meanShape, 0);
	//ssmUtils.cluster(1);

	//for (int i=0;i<meanShape->GetNumberOfPoints();i++)
	//{
	//	meanShape->SetPointData(i, ssmUtils.getShapeCluster(i)->clusterId);
	//}
	//km::writeMesh<MeshType>("clusteredMesh.vtk", meanShape);

	system("pause");
	
	return 0;
}