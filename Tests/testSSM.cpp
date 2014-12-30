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

	const int Dimension = 3;
	typedef double MeshPixelType;
	typedef itk::SimplexMeshRepresenter<MeshPixelType, Dimension> RepresenterType;
	typedef itk::StatisticalModel<RepresenterType>                StatisticalModelType;
	typedef RepresenterType::MeshType                             MeshType;
	
	StatisticalModelType::Pointer model = StatisticalModelType::New();
	model->Load( SSMFile );

	typedef itk::AffineTransform<double, Dimension> RigidTransformType;
	typedef itk::StatisticalShapeModelTransform<RepresenterType, double, Dimension> ShapeTransformType;
	typedef km::SSMUtils<MeshType, StatisticalModelType, RigidTransformType, ShapeTransformType> SSMUtilsType;
	SSMUtilsType ssmUtils;
	ssmUtils.SetSSM(model);

	MeshType::Pointer meanShape = model->DrawMean();
	km::assigneMesh<MeshType>(meanShape, 0);
	ssmUtils.cluster(1);

	for (int i=0;i<meanShape->GetNumberOfPoints();i++)
	{
		meanShape->SetPointData(i, ssmUtils.getShapeCluster(i)->clusterId);
	}
	km::writeMesh<MeshType>("clusteredMesh.vtk", meanShape);

	system("pause");
	
	return 0;
}