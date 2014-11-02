#include "itkHessianRecursiveGaussianImageFilter.h"
#include "itkSimpleFilterWatcher.h"

int main(int argc, char* argv[] )
{

	// Define the dimension of the images
	const unsigned int myDimension = 3;

	// Declare the types of the images
	typedef itk::Image<float, myDimension>           myImageType;

	// Declare the type of the index to access images
	typedef itk::Index<myDimension>             myIndexType;

	// Declare the type of the size
	typedef itk::Size<myDimension>              mySizeType;

	// Declare the type of the Region
	typedef itk::ImageRegion<myDimension>        myRegionType;

	// Create the image
	myImageType::Pointer inputImage  = myImageType::New();


	// Define their size, and start index
	mySizeType size;
	size[0] = 8;
	size[1] = 8;
	size[2] = 8;

	myIndexType start;
	start.Fill(0);

	myRegionType region;
	region.SetIndex( start );
	region.SetSize( size );

	// Initialize Image A
	inputImage->SetLargestPossibleRegion( region );
	inputImage->SetBufferedRegion( region );
	inputImage->SetRequestedRegion( region );
	inputImage->Allocate();

	// Declare Iterator type for the input image
	typedef itk::ImageRegionIteratorWithIndex<myImageType>  myIteratorType;

	// Create one iterator for the Input Image A (this is a light object)
	myIteratorType it( inputImage, inputImage->GetRequestedRegion() );

	// Initialize the content of Image A
	while( !it.IsAtEnd() )
	{
		it.Set( 0.0 );
		++it;
	}

	size[0] = 4;
	size[1] = 4;
	size[2] = 4;

	start[0] = 2;
	start[1] = 2;
	start[2] = 2;

	// Create one iterator for an internal region
	region.SetSize( size );
	region.SetIndex( start );
	myIteratorType itb( inputImage, region );

	// Initialize the content the internal region
	while( !itb.IsAtEnd() )
	{
		itb.Set( 100.0 );
		++itb;
	}

	// Declare the type for the
	typedef itk::HessianRecursiveGaussianImageFilter<
		myImageType >  myFilterType;

	typedef myFilterType::OutputImageType myHessianImageType;


	// Create a  Filter
	myFilterType::Pointer filter = myFilterType::New();
	itk::SimpleFilterWatcher watcher(filter);


	// Connect the input images
	filter->SetInput( inputImage );

	// Select the value of Sigma
	filter->SetSigma( 2.5 );


	// Execute the filter
	filter->Update();

	// Get the Smart Pointer to the Filter Output
	// It is important to do it AFTER the filter is Updated
	// Because the object connected to the output may be changed
	// by another during GenerateData() call
	myHessianImageType::Pointer outputImage = filter->GetOutput();

	// Declare Iterator type for the output image
	typedef itk::ImageRegionIteratorWithIndex<
		myHessianImageType>  myOutputIteratorType;

	// Create an iterator for going through the output image
	myOutputIteratorType itg( outputImage,
		outputImage->GetRequestedRegion() );

	//  Print the content of the result image
	std::cout << " Result " << std::endl;
	itg.GoToBegin();
	while( !itg.IsAtEnd() )
	{
		std::cout << itg.Get()<<std::endl;
		++itg;
	}

	// the following just tests for warnings in 2D
	typedef itk::Image<float,2>                                       my2DImageType;
	typedef itk::HessianRecursiveGaussianImageFilter<my2DImageType >  my2DFilterType;
	my2DFilterType::Pointer test = my2DFilterType::New();
	if(test.IsNull())
	{
		return EXIT_FAILURE;
	}

	// All objects should be automatically destroyed at this point
	return EXIT_SUCCESS;

}