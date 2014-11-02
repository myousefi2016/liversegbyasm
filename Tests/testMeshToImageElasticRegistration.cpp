#include "itkImage.h"
#include "itkMesh.h"

#include "kmUtility.h"
#include "kmRegistration.h"

int main(int argc, char* argv[])
{
	if (argc < 4)
	{
		std::cout<<"Usage: FixedImage MovingImage MovingMesh OutputMesh"<<std::endl;
		return EXIT_FAILURE;
	}
	
	const unsigned int Dimension = 3;
	typedef itk::Image< unsigned char, Dimension > ImageType;
	typedef itk::Mesh<double, Dimension>           MeshType;
	
	ImageType::Pointer fixedImage  = km::readImage<ImageType>( argv[1] );
	ImageType::Pointer movingImage = km::readImage<ImageType>( argv[2] );

	MeshType::Pointer  mesh  = km::readMesh<MeshType>( argv[3] );

	typedef itk::BSplineTransform<double,Dimension> TransformType;
	TransformType::Pointer transform = TransformType::New();

	ImageType::SpacingType spacing;
	spacing.Fill( 2.0 );
	fixedImage = km::resampleImage<ImageType>( fixedImage, spacing, 0 );

	km::elasticRegistration<ImageType, TransformType>( fixedImage, movingImage, transform,  12 );
	movingImage = km::transformImage<ImageType, TransformType>( movingImage, movingImage, transform, 0 );
	km::writeImage<ImageType>( "alignedImage-coarse.nii.gz", movingImage );

	mesh  = km::transformMesh<MeshType, TransformType>( mesh, transform );
	km::writeMesh<MeshType>( "alignedMesh-coarse.vtk", mesh );


	
	return EXIT_SUCCESS;
}