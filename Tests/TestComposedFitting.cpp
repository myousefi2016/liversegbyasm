#include <cmath>
#include <ctime>

#include "Representers/ITK/itkSimplexMeshRepresenter.h"
#include "statismo_ITK/itkStatisticalModel.h"
#include "statismo_ITK/itkStatisticalShapeModelTransform.h"
#include "itkSimilarity3DTransform.h"
#include "itkAffineTransform.h"

#include "kmUtility.h"
#include "kmModelFitting.h"

const unsigned int Dimension = 3;
typedef double MeshPixelType;
typedef itk::SimplexMeshRepresenter<MeshPixelType, Dimension> RepresenterType;
typedef itk::StatisticalModel<RepresenterType> StatisticalModelType;
typedef RepresenterType::MeshType MeshType;
typedef itk::StatisticalShapeModelTransform<RepresenterType, double, Dimension> ShapeTransformType;
//typedef itk::Similarity3DTransform<double> RigidTransformType;
typedef itk::AffineTransform<double, Dimension> RigidTransformType;
typedef itk::CompositeTransform<double, Dimension> CompositeTransformType;

typedef statismo::MatrixType StatismoMatrixType;
typedef statismo::VectorType StatismoVectorType;

void fillMatrix(StatismoMatrixType& matrix, const MeshType* mesh)
{
	typedef MeshType::PointType PointType;
	typedef MeshType::PointIdentifier PointIdentifier;
	typedef MeshType::PointsContainerConstIterator PointsContainerConstIterator;

	PointsContainerConstIterator ptIt = mesh->GetPoints()->Begin();
	PointsContainerConstIterator ptItEnd = mesh->GetPoints()->End();
	while(ptIt != ptItEnd)
	{
		PointIdentifier idx = ptIt.Index();
		PointType pt = ptIt.Value();
		for (int d=0;d<Dimension;d++)
		{
			matrix(idx, d) = pt[d];
		}
		ptIt++;
	}
}

void FillRigidTransform(const StatismoMatrixType & mat, RigidTransformType::Pointer & rigidTransform)
{
	RigidTransformType::MatrixType rigidMatrix = rigidTransform->GetMatrix();
	RigidTransformType::OffsetType rigidOffset = rigidTransform->GetOffset();

	for (int i=0;i<3;i++)
	{
		for (int j=0;j<3;j++)
		{
			rigidMatrix[i][j] = mat(i, j);
		}
	}
	std::cout<<rigidMatrix<<std::endl;

	for (int i=0;i<Dimension;i++)
	{
		rigidOffset[i] = mat(i, Dimension);
	}
	std::cout<<rigidOffset<<std::endl;

	rigidTransform->SetMatrix( rigidMatrix );
	rigidTransform->SetOffset( rigidOffset );
}

void FillShapeTransform(const StatismoMatrixType & mat, ShapeTransformType::ParametersType & params)
{
	
}

int main(int argc, char * argv[] )
{
	if(argc<5)
	{
		std::cout<<"Usage:"<<std::endl;
		std::cout<<"       ShapeModel(.h5) MeanMesh TestMesh FittedMesh";
		return EXIT_FAILURE;
	}
	
	const char* ssmFile = argv[1];
	const char* meanMeshFile = argv[2];
	const char* testMeshFile = argv[3];
	const char* fittedMeshFile = argv[4];
	
	srand(time(0));
	//(float)rand()/RAND_MAX
	
	//ÔØÈëSSM
	KM_DEBUG_PRINT( "Loading statistical shape model..", ssmFile );
	StatisticalModelType::Pointer model = StatisticalModelType::New();
	model->Load( ssmFile );
	
	ShapeTransformType::Pointer shapeTransform = ShapeTransformType::New();
	shapeTransform->SetStatisticalModel( model );
	shapeTransform->SetIdentity();
	
	MeshType::Pointer referenceShapeMesh = model->GetRepresenter()->GetReference();
	MeshType::Pointer meanShapeMesh = model->DrawMean();
	km::writeMesh<MeshType>( meanMeshFile, meanShapeMesh ); //Output mean mesh
	MeshType::PointType meshCentroid = km::getMeshCentroid<MeshType>( meanShapeMesh );
	
	ShapeTransformType::ParametersType shapeParam = shapeTransform->GetParameters();
	
	//Generate test shape parameters.
	for(int i=0;i<shapeTransform->GetNumberOfParameters();i++)
	{
		shapeParam[i] = ((float)rand()/RAND_MAX)*3.0-1.5;
	}
	std::cout<<"Test shape parameters: "<<shapeParam<<std::endl;
	shapeTransform->SetParameters(shapeParam);
	
	RigidTransformType::Pointer rigidTransform = RigidTransformType::New();
	rigidTransform->SetIdentity();
	rigidTransform->SetCenter( meshCentroid );

	typedef RigidTransformType::OutputVectorType OutputVectorType;
	OutputVectorType axisVec;
	axisVec[0] = 1.0; axisVec[1] = 0.0; axisVec[2] = 0.0;
	rigidTransform->Rotate3D(axisVec, 15.0);

	axisVec[0] = 0.0; axisVec[1] = 1.0; axisVec[2] = 0.0;
	rigidTransform->Rotate3D(axisVec, 15.0);

	axisVec[0] = 0.0; axisVec[1] = 0.0; axisVec[2] = 1.0;
	rigidTransform->Rotate3D(axisVec, 45.0);

	OutputVectorType scalesVec;
	scalesVec[0] = 0.9; scalesVec[1] = 1.1; scalesVec[2] = 1.3;
	rigidTransform->Scale(scalesVec);

	OutputVectorType transVec;
	transVec[0] = 50; transVec[1] = 50; transVec[2] = 100;
	rigidTransform->Translate(transVec);
 	
	std::cout<<"Test rigid parameters: \n"<<rigidTransform->GetMatrix()<<std::endl;
	std::cout<<"Test rigid parameters: \n"<<rigidTransform->GetOffset()<<std::endl;
	
	CompositeTransformType::Pointer compositeTransform = CompositeTransformType::New();
	compositeTransform->AddTransform( rigidTransform );
	compositeTransform->AddTransform( shapeTransform );
	
	MeshType::Pointer testMesh = km::transformMesh<MeshType, CompositeTransformType>( referenceShapeMesh, compositeTransform );
	km::writeMesh<MeshType>( testMeshFile, testMesh ); //Output mean mesh
	
	/**************************Start to fit*************************/
	km::compositeTransformFitting<MeshType, StatisticalModelType, RigidTransformType, ShapeTransformType>( testMesh, model, rigidTransform, shapeTransform );

	std::cout<<"Fitted rigid parameters: \n"<<rigidTransform->GetMatrix()<<std::endl;
	std::cout<<"Fitted rigid parameters: \n"<<rigidTransform->GetOffset()<<std::endl;
	std::cout<<"Fitted shape parameters: \n"<<shapeTransform->GetParameters()<<std::endl;
	
	std::cout<<"Transform fitted mesh..."<<std::endl;
	MeshType::Pointer fittedMesh = km::transformMesh<MeshType, CompositeTransformType>( referenceShapeMesh, compositeTransform );
	
	std::cout<<"Write fitted mesh..."<<std::endl;
	km::writeMesh<MeshType>( fittedMeshFile, fittedMesh ); //Output mean mesh
	
	return EXIT_SUCCESS;
}