/*
AdaBoost classification based on Dr. Hongzhi Wang's AdaBoost code.
*/


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

const int Dimension = 2;
typedef float PixelType;
typedef itk::Image< PixelType, Dimension >  ImageType;
typedef itk::ImageFileReader< ImageType > ReaderType;
typedef itk::ImageFileWriter< ImageType >  WriterType;
typedef float WritePixelType;
typedef itk::Image< WritePixelType, Dimension > WriteImageType;
typedef itk::ImageLinearConstIteratorWithIndex< ImageType >  ConstIteratorType;
typedef itk::NeighborhoodIterator< ImageType > NeighborhoodIteratorType;
typedef itk::ImageRegionIteratorWithIndex< ImageType > IndexIteratorType;
typedef itk::ImageRegionIterator< ImageType>        IteratorType;

typedef  unsigned char                         MaskPixelType;
typedef itk::Image< MaskPixelType, Dimension > MaskImageType;
typedef itk::ImageFileReader< MaskImageType >  MaskReaderType;
typedef itk::VectorImage<float, Dimension>     VectorImageType;
typedef itk::ImageFileWriter< MaskImageType >  MaskWriterType;
typedef itk::DiscreteGaussianImageFilter<ImageType, ImageType> GaussianFilterType;
typedef itk::DerivativeImageFilter<ImageType,ImageType>        DerivativeFilterType;
typedef itk::ShiftScaleImageFilter<ImageType,ImageType>        ShiftScaleFilterType;
typedef itk::GradientMagnitudeImageFilter<ImageType,ImageType> GradientMagnitudeFilterType;
typedef itk::ImageRegionIteratorWithIndex<VectorImageType>     VectorIteratorType;
typedef itk::ComposeImageFilter<ImageType>                     ImageToVectorImageFilterType;
typedef itk::ImageRegionIteratorWithIndex< MaskImageType >     MaskIndexIteratorType;



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

	GaussianFilterType::Pointer filter = GaussianFilterType::New();
	DerivativeFilterType::Pointer filterX = DerivativeFilterType::New();
	DerivativeFilterType::Pointer filterY = DerivativeFilterType::New();
	ShiftScaleFilterType::Pointer shiftScale = ShiftScaleFilterType::New();
	filterX->SetDirection(0);
	filterY->SetDirection(1);
	
	filter->SetInput(im1->GetOutput());
	filter->SetUseImageSpacingOn();

	char temp_dir[]="./";


	
	ImageToVectorImageFilterType::Pointer imageToVectorImageFilter = ImageToVectorImageFilterType::New();

	const int max_order = 2;
	const int no_scales = 3;
	float sigmas2[3]={0.25,1.0,4.0};


	int feat_n=0;
	for(int i=0; i<no_scales; i++)
	{
		for(int j=0; j<=max_order; j++)
		{
			// order j
			for(int ox=0; ox<=j; ox++)
			{
				for(int oy=0; oy<=j; oy++)
				{
					
						if(ox+oy==j)
						{
							filter->SetVariance(sigmas2[i]);
							filter->SetMaximumKernelWidth( 12*sigmas2[i] );
							filter->SetMaximumError(0.02);

							filter->Update();
							if(ox>0)
							{
								filterX->SetInput(filter->GetOutput());
								filterX->SetOrder(ox);
								filterX->Update();
							}
							if(oy>0)
							{
								if(ox>0)
								{
									filterY->SetInput(filterX->GetOutput());
								}
								else
								{
									filterY->SetInput(filter->GetOutput());
								}
								filterY->SetOrder(oy);
								filterY->Update();
							}
							


							ShiftScaleFilterType::Pointer shiftScale = ShiftScaleFilterType::New();

							ImageType::Pointer image = ImageType::New(); 


							
							if(oy>0)
							{
								filterY->Update();
								//		shiftScale->SetInput(filterY->GetOutput());
								image = filterY->GetOutput();
								image->DisconnectPipeline();
							} 
							else if(ox>0)
							{
								filterX->Update();
								//	shiftScale->SetInput(filterX->GetOutput());
								image = filterX->GetOutput();
								image->DisconnectPipeline();
							}
							else
							{
								filter->Update();
								//	shiftScale->SetInput(filter->GetOutput());
								image = filter->GetOutput();
								image->DisconnectPipeline();
							}

							imageToVectorImageFilter->SetInput(feat_n, image);
							feat_n++;
						
					}
				}
			}
		}
	}


	GaussianFilterType::Pointer filter1 = GaussianFilterType::New();
	GaussianFilterType::Pointer filter2 = GaussianFilterType::New();
	GaussianFilterType::Pointer filter3 = GaussianFilterType::New();

	filter1->SetInput(im1->GetOutput());
	filter2->SetInput(im1->GetOutput());
	filter3->SetInput(im1->GetOutput());


	filter1->SetVariance(sigmas2[0]);
	filter1->SetMaximumKernelWidth( 12*sigmas2[0] );
	filter1->SetMaximumError(0.02);

	filter2->SetVariance(sigmas2[1]);
	filter2->SetMaximumKernelWidth( 12*sigmas2[1] );
	filter2->SetMaximumError(0.02);

	filter3->SetVariance(sigmas2[2]);
	filter3->SetMaximumKernelWidth( 12*sigmas2[2] );
	filter3->SetMaximumError(0.02);


	GradientMagnitudeFilterType::Pointer gradientfilter1 = GradientMagnitudeFilterType::New();
	GradientMagnitudeFilterType::Pointer gradientfilter2 = GradientMagnitudeFilterType::New();
	GradientMagnitudeFilterType::Pointer gradientfilter3 = GradientMagnitudeFilterType::New();

	gradientfilter1->SetInput(filter1->GetOutput());
	gradientfilter1->Update();

	gradientfilter2->SetInput(filter2->GetOutput());
	gradientfilter2->Update();

	gradientfilter3->SetInput(filter3->GetOutput());
	gradientfilter3->Update();


	ImageType::Pointer gradient1 = gradientfilter1->GetOutput();
	ImageType::Pointer gradient2 = gradientfilter2->GetOutput();
	ImageType::Pointer gradient3 = gradientfilter3->GetOutput();


	gradient1->DisconnectPipeline();
	gradient2->DisconnectPipeline();
	gradient3->DisconnectPipeline();


	imageToVectorImageFilter->SetInput(feat_n++,gradient1);
	imageToVectorImageFilter->SetInput(feat_n++,gradient2);
	imageToVectorImageFilter->SetInput(feat_n++,gradient3);


	// location features
	ImageType::Pointer outputImageX = ImageType::New();
	ImageType::Pointer outputImageY = ImageType::New();
	outputImageX->SetRegions( im1->GetOutput()->GetRequestedRegion() );
	outputImageY->SetRegions( im1->GetOutput()->GetRequestedRegion() );
	outputImageX->CopyInformation( im1->GetOutput() );
	outputImageY->CopyInformation( im1->GetOutput() );
	outputImageX->Allocate();
	outputImageY->Allocate();

	IteratorType outputItX( outputImageX, outputImageX->GetRequestedRegion() );
	IteratorType outputItY( outputImageY, outputImageY->GetRequestedRegion() );
	ImageType::IndexType requestedIndexX =outputImageX->GetRequestedRegion().GetIndex();
	ImageType::IndexType requestedIndexY =outputImageY->GetRequestedRegion().GetIndex();
	ImageType::SizeType requestedSize = outputImageX->GetRequestedRegion().GetSize();
	for ( outputItX.GoToBegin(),outputItY.GoToBegin() ; !outputItX.IsAtEnd(); ++outputItX,++outputItY)
	{
		ImageType::IndexType idx = outputItX.GetIndex();
		outputItX.Set(idx[0]+1);
		outputItY.Set(idx[1]+1);
	}


	imageToVectorImageFilter->SetInput(feat_n++,outputImageX);
	imageToVectorImageFilter->SetInput(feat_n++,outputImageY);
	imageToVectorImageFilter->Update();

	MaskIndexIteratorType nsegit(nseg, nseg->GetRequestedRegion());
	

	for (nsegit.GoToBegin(); !nsegit.IsAtEnd(); ++nsegit)
	{
		nsegit.Set(0);
	}

	FLOATTYPE H,cH;
	
	int NFeature =  23;

	VectorIteratorType vit(imageToVectorImageFilter->GetOutput(),imageToVectorImageFilter->GetOutput()->GetRequestedRegion());


	VectorImageType::PixelType pixel;
	VectorImageType::IndexType index;
	int fgcounter=0;


	vit.GoToBegin();
	nsegit.GoToBegin();
	while(!vit.IsAtEnd())
	{
		vector<FLOATTYPE> X;
		pixel = vit.Get();
		index = vit.GetIndex();
		for(int i=0;i<NFeature;i++)
		{
			X.push_back(pixel[i]);
		}

		//Apply AdaBoost classifier on the features to obtain classifier response in probability, stored in cH.
		AdaBoostClassify(&X[0], 1, NFeature, LC, &featureID[0], &alpha[0], &sign[0], &threshold[0], &H, &cH);
		if (H>0.5)
			nsegit.Set(1);

		++vit;	
		++nsegit;
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
