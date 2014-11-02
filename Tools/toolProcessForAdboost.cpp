#include <iostream>
#include "kmUtility.h"
#include "kmProcessing.h"
#include "kmRegistration.h"
#include <itkNormalizeImageFilter.h>
#include "itkBinaryContourImageFilter.h"

int main(int argc, char* argv[])
{	
	if(argc<6)
	{
		std::cout<<"Usage: origimage livermask orig_output livermask_output skinmask_output livercontour_output"<<std::endl;
		return EXIT_FAILURE;
	}

	typedef itk::Image<signed short, 3> ShortImageType;
	typedef itk::Image<unsigned char, 3> UCharImageType;
	typedef itk::Image<float, 3> FloatImageType;

	ShortImageType::Pointer origImage = km::readImage<ShortImageType>( argv[1] );
	UCharImageType::Pointer livermask = km::readImage<UCharImageType>( argv[2] );

	double minval, maxval;
	km::calculateMinAndMax<ShortImageType>(origImage, minval, maxval);

	ShortImageType::SpacingType downspac;
	downspac.Fill( 3.0 );

	ShortImageType::Pointer origimageOutput = km::resampleImage<ShortImageType>( origImage, downspac, minval, 0 );
	UCharImageType::Pointer livermaskOutput = km::resampleImage<UCharImageType>( livermask, downspac, 0, 1 );

	UCharImageType::Pointer skinmaskOutput = UCharImageType::New();
	km::removeAir<ShortImageType, UCharImageType>( origImage,skinmaskOutput, minval, true );

	//typedef itk::NormalizeImageFilter< ShortImageType, FloatImageType >
	//	NormalizeFilterType;
	//NormalizeFilterType::Pointer normalizeFilter = NormalizeFilterType::New();
	//normalizeFilter->SetInput( origdown );
	//normalizeFilter->Update();

	//FloatImageType::Pointer origimageOutput = normalizeFilter->GetOutput();

	typedef itk::BinaryContourImageFilter<UCharImageType, UCharImageType> BinaryContourFilterType;
	BinaryContourFilterType::Pointer contourFileter = BinaryContourFilterType::New();
	contourFileter->SetInput( livermaskOutput );
	contourFileter->SetForegroundValue(1);
	contourFileter->Update();

	UCharImageType::Pointer livercontourOutput = contourFileter->GetOutput();

	km::writeImage<ShortImageType>( argv[3], origimageOutput );
	km::writeImage<UCharImageType>( argv[4], livermaskOutput );
	km::writeImage<UCharImageType>( argv[5], skinmaskOutput );
	km::writeImage<UCharImageType>( argv[6], livercontourOutput );

	return EXIT_SUCCESS;
}