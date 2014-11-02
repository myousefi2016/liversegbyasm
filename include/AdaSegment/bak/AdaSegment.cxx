/*
AdaBoost classification based on Dr. Hongzhi Wang's AdaBoost code.
*/

#include "util.h"
#include "AdaBoost.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
#include "itkCommand.h"
#include "itkImageMaskSpatialObject.h"
#include "itkResampleImageFilter.h"
#include "itkMaskImageFilter.h"
#include "itkDivideImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkDerivativeImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkShiftScaleImageFilter.h"
#include "itkImageRandomConstIteratorWithIndex.h"
#include "itkGradientMagnitudeImageFilter.h"
#include "itkBinaryImageToShapeLabelMapFilter.h"
#include "itkComposeImageFilter.h"
#include "itkStreamingImageFilter.h"
#include "itkStreamingImageFilter.h"
#include "itkImageRegionSplitterMultidimensional.h"
#include "itkPipelineMonitorImageFilter.h"
#include "itkHessianRecursiveGaussianImageFilter.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include "itkDiscreteGaussianDerivativeImageFilter.h"


typedef  unsigned char                         MaskPixelType;
typedef itk::Image< MaskPixelType, Dimension > MaskImageType;
typedef itk::Image< PixelType, Dimension >     ImageType;
typedef itk::ImageFileReader< MaskImageType >  MaskReaderType;
typedef itk::ImageFileReader< ImageType >      ReaderType;
typedef itk::ImageFileWriter< ImageType >      WriterType;
typedef itk::VectorImage<float, 3>             VectorImageType;
typedef itk::ImageFileWriter< MaskImageType >  MaskWriterType;
typedef itk::DiscreteGaussianImageFilter<ImageType, ImageType> GaussianFilterType;
typedef itk::DerivativeImageFilter<ImageType,ImageType>        DerivativeFilterType;
typedef itk::ShiftScaleImageFilter<ImageType,ImageType>        ShiftScaleFilterType;
typedef itk::GradientMagnitudeImageFilter<ImageType,ImageType> GradientMagnitudeFilterType;
typedef itk::ImageRegionIteratorWithIndex<VectorImageType>     VectorIteratorType;
typedef itk::ComposeImageFilter<ImageType>                     ImageToVectorImageFilterType;
typedef itk::ImageRegionIteratorWithIndex< MaskImageType >     MaskIndexIteratorType;
typedef itk::ImageRegionConstIteratorWithIndex<ImageType>      ImageIteratorType;



using namespace std;

int main( int argc, char ** argv )
{
	if ( argc < 4 )
	{
		std::cerr << "Missing parameters. " << std::endl;
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0]
		<< " inputImage AdaBoostOutPutPrefix outputSegmentation [-x label image.nii]"<<endl<<endl
			<< " Location of organ based on multi-features AdaBoost classifiers."<<endl<<endl
			<< " Meanings of the parameters:"<<endl 
			<< " inputImage:           the image file" <<endl
			<< " AdaBoostOutPutPrefix: the path and prefix to the learned AdaBoost files used in ./bl"<<endl
			<< " outputSegmentation:   the path and file name of the corrected segmentation."<<endl
			<< " -x label image.nii    Specify an exclusion region for the given label. " << endl
			<< "                       If a voxel has non-zero value in an exclusion image," << endl
			<< "                       the corresponding label is not allowed at that voxel." << endl
			<< std::endl;
		return -1;
	}
	
	ReaderType::Pointer im1 = ReaderType::New();
	im1->SetFileName( argv[1] );
	try
	{
		im1->Update();
	}
	catch ( itk::ExceptionObject &err)
	{
		std::cout << "ExceptionObject caught !" << std::endl;
		std::cout << err << std::endl;
		return -1;
	}



	vector<int> sign;
	vector<FLOATTYPE> alpha;
	vector<FLOATTYPE> threshold;
	vector<int> featureID;


	char tfn[1024];
	sprintf(tfn,"%s-AdaBoostResults",argv[2]);
	std::cout<<tfn<<std::endl;
	ifstream ifs; 
	ifs.open( tfn , ifstream::in );

	FLOATTYPE t;
	int LC;

	while (ifs.good())
	{
		ifs >> t;
		if (!ifs.good())
			break;
		ifs >> t;
		alpha.push_back(t);
		ifs >> t;
		ifs >> t;
		featureID.push_back(int(t));
		ifs >> t;
		sign.push_back(int(t));
		ifs >> t;
		threshold.push_back(t);
		ifs >> t;
	}
	ifs.close();
	LC = alpha.size();
	cout<<"# weak learners: "<<LC<<endl;


	MaskImageType::Pointer nseg = MaskImageType::New();
	nseg->SetRegions(im1->GetOutput()->GetRequestedRegion());
	nseg->SetSpacing( im1->GetOutput()->GetSpacing() );
	nseg->SetOrigin( im1->GetOutput()->GetOrigin() );
	nseg->SetDirection(im1->GetOutput()->GetDirection());
	nseg->Allocate();




	double maxError = 0.02;
	double maxKernelWidth = 100;

	unsigned int orderX[3],orderY[3],orderZ[3];
	orderX[0] = 1; orderX[1] = 0; orderX[2] = 0;
	orderY[0] = 0; orderY[1] = 1; orderY[2] = 0;
	orderZ[0] = 0; orderY[1] = 0; orderZ[2] = 1;
	

	float sigma[3]={0.7071,1.4142,2.8284};

	std::vector<ImageType::Pointer> imagefeatures;


	for(int scale=0;scale<3;scale++)
	{

		typedef itk::DiscreteGaussianDerivativeImageFilter<ImageType,ImageType> GaussianDerivativeFilterType;

		GaussianDerivativeFilterType::Pointer derivativeX =	GaussianDerivativeFilterType::New();
		derivativeX->SetInput( im1->GetOutput() );
		derivativeX->SetMaximumError( maxError );
		derivativeX->SetMaximumKernelWidth( maxKernelWidth );
		derivativeX->SetOrder( orderX );
		derivativeX->SetNormalizeAcrossScale( false );
		derivativeX->SetUseImageSpacing( true );
		derivativeX->SetVariance( sigma[scale]*sigma[scale] );
		
		GaussianDerivativeFilterType::Pointer derivativeY =	GaussianDerivativeFilterType::New();
		derivativeY->SetInput( im1->GetOutput() );
		derivativeY->SetMaximumError( maxError );
		derivativeY->SetMaximumKernelWidth( maxKernelWidth );
		derivativeY->SetOrder( orderY );
		derivativeY->SetNormalizeAcrossScale( false );
		derivativeY->SetUseImageSpacing( true );
		derivativeY->SetVariance( sigma[scale]*sigma[scale] );
		
		GaussianDerivativeFilterType::Pointer derivativeZ =	GaussianDerivativeFilterType::New();
		derivativeZ->SetInput( im1->GetOutput() );
		derivativeZ->SetMaximumError( maxError );
		derivativeZ->SetMaximumKernelWidth( maxKernelWidth );
		derivativeZ->SetOrder( orderZ );
		derivativeZ->SetNormalizeAcrossScale( false );
		derivativeZ->SetUseImageSpacing( true );
		derivativeZ->SetVariance( sigma[scale]*sigma[scale] );




		typedef itk::RecursiveGaussianImageFilter<ImageType >  RecursiveGaussianImageFilterType;

		// derivative of x^2
		RecursiveGaussianImageFilterType::Pointer gaussianFilter_x_x = RecursiveGaussianImageFilterType::New();
		gaussianFilter_x_x->SetInput( im1->GetOutput() );
		gaussianFilter_x_x->SetSigma(sigma[scale]);
		gaussianFilter_x_x->SetDirection(0); // "x" axis
		gaussianFilter_x_x->SetSecondOrder();


		// derivative of y^2
		RecursiveGaussianImageFilterType::Pointer gaussianFilter_y_y = RecursiveGaussianImageFilterType::New();
		gaussianFilter_y_y->SetInput( im1->GetOutput() );
		gaussianFilter_y_y->SetSigma(sigma[scale]);
		gaussianFilter_y_y->SetDirection(1); // "y" axis
		gaussianFilter_y_y->SetSecondOrder();

		// derivative of z^2
		RecursiveGaussianImageFilterType::Pointer gaussianFilter_z_z = RecursiveGaussianImageFilterType::New();
		gaussianFilter_z_z->SetInput( im1->GetOutput() );
		gaussianFilter_z_z->SetSigma(sigma[scale]);
		gaussianFilter_z_z->SetDirection(2); // "z" axis
		gaussianFilter_z_z->SetSecondOrder();


	    // derivative of x_y
		RecursiveGaussianImageFilterType::Pointer gaussianFilter_x_y = RecursiveGaussianImageFilterType::New();
		gaussianFilter_x_y->SetInput( derivativeX->GetOutput() );
		gaussianFilter_x_y->SetSigma(sigma[scale]);
		gaussianFilter_x_y->SetDirection(1); // "y" axis
		gaussianFilter_x_y->SetFirstOrder();


		// derivative of y_z
		RecursiveGaussianImageFilterType::Pointer gaussianFilter_y_z = RecursiveGaussianImageFilterType::New();
		gaussianFilter_y_z->SetInput( derivativeY->GetOutput() );
		gaussianFilter_y_z->SetSigma(sigma[scale]);
		gaussianFilter_y_z->SetDirection(2); // "z" axis
		gaussianFilter_y_z->SetFirstOrder();

		// derivative of x_z
		RecursiveGaussianImageFilterType::Pointer gaussianFilter_x_z = RecursiveGaussianImageFilterType::New();
		gaussianFilter_x_z->SetInput( derivativeX->GetOutput() );
		gaussianFilter_x_z->SetSigma(sigma[scale]);
		gaussianFilter_x_z->SetDirection(2); // "z" axis
		gaussianFilter_x_z->SetFirstOrder();

	
		typedef itk::GradientMagnitudeRecursiveGaussianImageFilter<ImageType> GradientMagnitudeFilterType;
		GradientMagnitudeFilterType::Pointer gradientmagnitudefunction = GradientMagnitudeFilterType::New();
		gradientmagnitudefunction->SetInput( im1->GetOutput() );
		gradientmagnitudefunction->SetNormalizeAcrossScale( false );
		gradientmagnitudefunction->SetSigma( sigma[scale] );

		typedef itk::SmoothingRecursiveGaussianImageFilter<ImageType> GaussianBlurFilterType;
		GaussianBlurFilterType::Pointer gaussianblur = GaussianBlurFilterType::New();

		gaussianblur->SetInput(im1->GetOutput());
		gaussianblur->SetNormalizeAcrossScale( false);
		gaussianblur->SetSigma(sigma[scale]);
	    
		derivativeX->Update();
		derivativeY->Update();
		derivativeZ->Update();
		gaussianFilter_x_x->Update();
		gaussianFilter_x_y->Update();
		gaussianFilter_x_z->Update();
		gaussianFilter_y_y->Update();
		gaussianFilter_y_z->Update();
		gaussianFilter_z_z->Update();
        
		gradientmagnitudefunction->Update();
		gaussianblur->Update();


		ImageType::Pointer image_derivativeX = derivativeX->GetOutput();
		image_derivativeX->DisconnectPipeline();
		imagefeatures.push_back(image_derivativeX);

		ImageType::Pointer image_derivativeY = derivativeY->GetOutput();
		image_derivativeY->DisconnectPipeline();
		imagefeatures.push_back(image_derivativeY);


		ImageType::Pointer image_derivativeZ = derivativeZ->GetOutput();
		image_derivativeZ->DisconnectPipeline();
		imagefeatures.push_back(image_derivativeZ);

		ImageType::Pointer image_derivativeXX = gaussianFilter_x_x->GetOutput();
		image_derivativeXX->DisconnectPipeline();
		imagefeatures.push_back(image_derivativeXX);

		ImageType::Pointer image_derivativeXY = gaussianFilter_x_y->GetOutput();
		image_derivativeXY->DisconnectPipeline();
		imagefeatures.push_back(image_derivativeXY);

		ImageType::Pointer image_derivativeXZ = gaussianFilter_x_z->GetOutput();
		image_derivativeXZ->DisconnectPipeline();
		imagefeatures.push_back(image_derivativeXZ);

		ImageType::Pointer image_derivativeYY = gaussianFilter_y_y->GetOutput();
		image_derivativeYY->DisconnectPipeline();
		imagefeatures.push_back(image_derivativeYY);

		ImageType::Pointer image_derivativeYZ = gaussianFilter_y_z->GetOutput();
		image_derivativeYZ->DisconnectPipeline();
		imagefeatures.push_back(image_derivativeYZ);

		ImageType::Pointer image_derivativeZZ = gaussianFilter_z_z->GetOutput();
		image_derivativeZZ->DisconnectPipeline();
		imagefeatures.push_back(image_derivativeZZ);

		ImageType::Pointer image_gradientmagnitude = gradientmagnitudefunction->GetOutput();
		image_gradientmagnitude->DisconnectPipeline();
		imagefeatures.push_back(image_gradientmagnitude);

		ImageType::Pointer image_gaussianblur = gaussianblur->GetOutput();
		image_gaussianblur->DisconnectPipeline();
		imagefeatures.push_back(image_gaussianblur);

	
	}


	const int max_order = 2;
	const int no_scales = 3;
	

	// location features
	ImageType::Pointer outputImageX = ImageType::New();
	ImageType::Pointer outputImageY = ImageType::New();
	ImageType::Pointer outputImageZ = ImageType::New();
	outputImageX->SetRegions( im1->GetOutput()->GetRequestedRegion() );
	outputImageY->SetRegions( im1->GetOutput()->GetRequestedRegion() );
	outputImageZ->SetRegions( im1->GetOutput()->GetRequestedRegion() );
	outputImageX->CopyInformation( im1->GetOutput() );
	outputImageY->CopyInformation( im1->GetOutput() );
	outputImageZ->CopyInformation( im1->GetOutput() );
	outputImageX->Allocate();
	outputImageY->Allocate();
	outputImageZ->Allocate();

	IteratorType outputItX( outputImageX, outputImageX->GetRequestedRegion() );
	IteratorType outputItY( outputImageY, outputImageY->GetRequestedRegion() );
	IteratorType outputItZ( outputImageZ, outputImageZ->GetRequestedRegion() );
	ImageType::IndexType requestedIndexX =outputImageX->GetRequestedRegion().GetIndex();
	ImageType::IndexType requestedIndexY =outputImageY->GetRequestedRegion().GetIndex();
	ImageType::IndexType requestedIndexZ =outputImageZ->GetRequestedRegion().GetIndex();
	ImageType::SizeType requestedSize = outputImageX->GetRequestedRegion().GetSize();
	for ( outputItX.GoToBegin(),outputItY.GoToBegin(),outputItZ.GoToBegin() ; !outputItX.IsAtEnd(); ++outputItX,++outputItY,++outputItZ)
	{
		ImageType::IndexType idx = outputItX.GetIndex();
		outputItX.Set(idx[0]);
		outputItY.Set(idx[1]);
		outputItZ.Set(idx[2]);
	}

	imagefeatures.push_back(outputImageX);
	imagefeatures.push_back(outputImageY);
	imagefeatures.push_back(outputImageZ);


	MaskIndexIteratorType nsegit(nseg, nseg->GetRequestedRegion());
	

	for (nsegit.GoToBegin(); !nsegit.IsAtEnd(); ++nsegit)
	{
		nsegit.Set(0);
	}

	FLOATTYPE H,cH;
	
	int NFeature =  36;

	std::vector<ImageIteratorType> imagefeaturesiterator;

	for(int i=0;i<imagefeatures.size();i++)
	{
		ImageIteratorType it(imagefeatures[i],imagefeatures[i]->GetRequestedRegion());
		it.GoToBegin();
		imagefeaturesiterator.push_back(it);
	}

	nsegit.GoToBegin();
	while(!nsegit.IsAtEnd())
	{
		vector<FLOATTYPE> X;
		for(int i=0;i<NFeature;i++)
		{
			X.push_back(imagefeaturesiterator[i].Get());
		}

		//Apply AdaBoost classifier on the features to obtain classifier response in probability, stored in cH.
		AdaBoostClassify(&X[0], 1, NFeature, LC, &featureID[0], &alpha[0], &sign[0], &threshold[0], &H, &cH);
		if (H>0.5)
			nsegit.Set(1);

		++nsegit;	
		for(int i=0;i<NFeature;i++)
			++imagefeaturesiterator[i];

	}




	cout<<"saving the segmentation to "<<argv[3]<<endl;
	// write the segmentation into the output file
	MaskWriterType::Pointer writer = MaskWriterType::New();
	writer->SetFileName( argv[3] );
	writer->SetInput(nseg);
	try
	{
		writer->Update();
	}
	catch ( itk::ExceptionObject &err)
	{
		std::cout << "ExceptionObject caught !" << std::endl;
		std::cout << err << std::endl;
		return -1;
	}

	return 0;
}
