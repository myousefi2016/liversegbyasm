#include "itkImage.h"
#include "itkMesh.h"

#include "kmUtility.h"
#include "kmRegistration.h"

int main(int argc, char* argv[])
{
	if (argc < 4)
	{
		std::cout<<"Usage: FixedImage MovingImage OutputImage"<<std::endl;
		return EXIT_FAILURE;
	}
	
	const unsigned int Dimension = 3;
	typedef itk::Image< signed short, Dimension > ImageType;
	
	ImageType::Pointer fixedImage = km::readImage<ImageType>( argv[1] );
	ImageType::Pointer movingImage = km::readImage<ImageType>( argv[2] );

	typedef itk::BSplineTransform<double,Dimension> TransformType;
	TransformType::Pointer transform = TransformType::New();

	ImageType::SpacingType spacing;
	spacing.Fill( 2.0 );
	fixedImage = km::resampleImage<ImageType>( fixedImage, spacing, 0 );

	km::elasticRegistration<ImageType, TransformType>( fixedImage, movingImage, transform,  6 );
	movingImage = km::transformImageByReference<ImageType, ImageType, TransformType>( movingImage, movingImage, transform, 0 );
	km::writeImage<ImageType>( "aligned-coarse.nii.gz", movingImage );
	km::elasticRegistration<ImageType, TransformType>( fixedImage, movingImage, transform,  8 );
	movingImage = km::transformImageByReference<ImageType, ImageType, TransformType>( movingImage, movingImage, transform, 0 );
	km::writeImage<ImageType>( "aligned-fine.nii.gz", movingImage );
	
	return EXIT_SUCCESS;
}