
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


const unsigned int Dimension = 3;
typedef float PixelType;
typedef itk::Image< PixelType, Dimension >  ImageType;
typedef itk::ImageFileReader< ImageType > ReaderType;
typedef itk::ImageFileWriter< ImageType >  WriterType;
typedef float WritePixelType;
typedef itk::Image< WritePixelType, 3 > WriteImageType;
typedef itk::ImageLinearConstIteratorWithIndex< ImageType >  ConstIteratorType;
typedef itk::NeighborhoodIterator< ImageType > NeighborhoodIteratorType;
typedef itk::ImageRegionIteratorWithIndex< ImageType > IndexIteratorType;
typedef itk::ImageRegionIterator< ImageType>        IteratorType;

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

using namespace std;



int main(int argc, char *argv[])
{

	if(argc<2)
	{
		std::cout<<"Usage: "<<argv[0]<<" input3dImage"<<std::endl;
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

	GaussianFilterType::Pointer filter = GaussianFilterType::New();
	DerivativeFilterType::Pointer filterX = DerivativeFilterType::New();
	DerivativeFilterType::Pointer filterY = DerivativeFilterType::New();
	DerivativeFilterType::Pointer filterZ = DerivativeFilterType::New();
	ShiftScaleFilterType::Pointer shiftScale = ShiftScaleFilterType::New();
	filterX->SetDirection(0);
	filterY->SetDirection(1);
	filterZ->SetDirection(2);

	filter->SetInput(im1->GetOutput());
	filter->SetUseImageSpacingOn();

	char temp_dir[]="./";

	char newname[1000];

	//ImageToVectorImageFilterType::Pointer imageToVectorImageFilter = ImageToVectorImageFilterType::New();

//	VectorImageType::Pointer vectorImage = imageToVectorImageFilter->GetOutput();

	const int max_order = 2;
	const int no_scales = 3;
	float sigmas2[3]={0.25,1.0,4.0};


	int feat_n=0;
	for(int i=0; i<no_scales; i++){
		for(int j=0; j<=max_order; j++){
			// order j
			for(int ox=0; ox<=j; ox++){
				for(int oy=0; oy<=j; oy++){
					for(int oz=0; oz<=j; oz++){
						if(ox+oy+oz==j){
							filter->SetVariance(sigmas2[i]);
							filter->SetMaximumKernelWidth( 12*sigmas2[i] );
							filter->SetMaximumError(0.02);

							filter->Update();
							if(ox>0){
								filterX->SetInput(filter->GetOutput());
								filterX->SetOrder(ox);
								filterX->Update();
							}
							if(oy>0){
								if(ox>0){
									filterY->SetInput(filterX->GetOutput());
								}
								else{
									filterY->SetInput(filter->GetOutput());
								}
								filterY->SetOrder(oy);
								filterY->Update();
							}
							if(oz>0){
								if(oy>0){
									filterZ->SetInput(filterY->GetOutput());
								}
								else if(ox>0){
									filterZ->SetInput(filterY->GetOutput());
								}
								else{
									filterZ->SetInput(filter->GetOutput());
								}
								filterZ->SetOrder(oz);
								filterZ->Update();
							}


							ShiftScaleFilterType::Pointer shiftScale = ShiftScaleFilterType::New();

							ImageType::Pointer image = ImageType::New(); 


							if(oz>0)
							{
								filterZ->Update();
								//shiftScale->SetInput(filterZ->GetOutput());
								image = filterZ->GetOutput();
								image->DisconnectPipeline();
								WriterType::Pointer writer = WriterType::New();
								sprintf(newname,"%s_%02d.mhd",argv[1],feat_n);
								writer->SetFileName(newname);
								writer->SetInput(image);
								writer->Update();
								std::cout<<"The "<<feat_n+1<<"th feature image has been written!"<<std::endl;

							} 
							else if(oy>0)
							{
								filterY->Update();
								//		shiftScale->SetInput(filterY->GetOutput());
								image = filterY->GetOutput();
								image->DisconnectPipeline();
								WriterType::Pointer writer = WriterType::New();
								sprintf(newname,"%s_%02d.mhd",argv[1],feat_n);
								writer->SetFileName(newname);
								writer->SetInput(image);
								writer->Update();
								std::cout<<"The "<<feat_n+1<<"th feature image has been written!"<<std::endl;
							} 
							else if(ox>0)
							{
								filterX->Update();
								//	shiftScale->SetInput(filterX->GetOutput());
								image = filterX->GetOutput();
								image->DisconnectPipeline();
								WriterType::Pointer writer = WriterType::New();
								sprintf(newname,"%s_%02d.mhd",argv[1],feat_n);
								writer->SetFileName(newname);
								writer->SetInput(image);
								writer->Update();
								std::cout<<"The "<<feat_n+1<<"th feature image has been written!"<<std::endl;
							}
							else
							{
								filter->Update();
								//	shiftScale->SetInput(filter->GetOutput());
								image = filter->GetOutput();
								image->DisconnectPipeline();
								WriterType::Pointer writer = WriterType::New();
								sprintf(newname,"%s_%02d.mhd",argv[1],feat_n);
								writer->SetFileName(newname);
								writer->SetInput(image);
								writer->Update();
								std::cout<<"The "<<feat_n+1<<"th feature image has been written!"<<std::endl;
							}

						//	imageToVectorImageFilter->SetInput(feat_n, image);

							feat_n++;
						}
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

	WriterType::Pointer writer1 = WriterType::New();
	sprintf(newname,"%s_%02d.mhd",argv[1],30);
	writer1->SetFileName(newname);
	writer1->SetInput(gradient1);
	writer1->Update();
    std::cout<<"The 31th feature image has been written!"<<std::endl;

	WriterType::Pointer writer2 = WriterType::New();
	sprintf(newname,"%s_%02d.mhd",argv[1],31);
	writer2->SetFileName(newname);
	writer2->SetInput(gradient2);
	writer2->Update();
	std::cout<<"The 32th feature image has been written!"<<std::endl;

	WriterType::Pointer writer3 = WriterType::New();
	sprintf(newname,"%s_%02d.mhd",argv[1],32);
	writer3->SetFileName(newname);
	writer3->SetInput(gradient3);
	writer3->Update();
    std::cout<<"The 33th feature image has been written!"<<std::endl;

//	imageToVectorImageFilter->SetInput(30,gradient1);
//	imageToVectorImageFilter->SetInput(31,gradient2);
//	imageToVectorImageFilter->SetInput(32,gradient3);


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
		outputItX.Set(idx[0]+1);
		outputItY.Set(idx[1]+1);
		outputItZ.Set(idx[2]+1);
	}


	//imageToVectorImageFilter->SetInput(33,outputImageX);
	//imageToVectorImageFilter->SetInput(34,outputImageY);
	//imageToVectorImageFilter->SetInput(35,outputImageZ);
	//imageToVectorImageFilter->Update();

	WriterType::Pointer writer4 = WriterType::New();
	sprintf(newname,"%s_%02d.mhd",argv[1],33);
	writer4->SetFileName(newname);
	writer4->SetInput(outputImageX);
	writer4->Update();
    std::cout<<"The 34th feature image has been written!"<<std::endl;

	WriterType::Pointer writer5 = WriterType::New();
	sprintf(newname,"%s_%02d.mhd",argv[1],34);
	writer5->SetFileName(newname);
	writer5->SetInput(outputImageY);
	writer5->Update();
    std::cout<<"The 35th feature image has been written!"<<std::endl;

	WriterType::Pointer writer6 = WriterType::New();
	sprintf(newname,"%s_%02d.mhd",argv[1],35);
	writer6->SetFileName(newname);
	writer6->SetInput(outputImageZ);
	writer6->Update();
	std::cout<<"The 36th feature image has been written!"<<std::endl;


	return 0;
}