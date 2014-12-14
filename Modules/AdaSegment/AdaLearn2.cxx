#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkMeanImageFunction.h"
#include "itkMedianImageFunction.h"
#include "itkVarianceImageFunction.h"
#include "itkSumOfSquaresImageFunction.h"
#include "itkCentralDifferenceImageFunction.h"
#include "itkGaussianBlurImageFunction.h"
#include "itkGaussianDerivativeImageFunction.h"
#include "itkDiscreteGaussianDerivativeImageFunction.h"
#include "itkDiscreteGradientMagnitudeGaussianImageFunction.h"
#include "itkDiscreteHessianGaussianImageFunction.h"
#include "itkSobelEdgeDetectionImageFilter.h"
#include "itkZeroCrossingImageFilter.h"
#include "itkBoxMeanImageFilter.h"
#include "itkBoxSigmaImageFilter.h"
#include "itkBox2SumImageFilter.h"
#include "itkLaplacianRecursiveGaussianImageFilter.h"
#include "itkHessian3DToVesselnessMeasureImageFilter.h"
#include "itkHessianRecursiveGaussianImageFilter.h"
#include "itkDerivativeImageFilter.h"
#include "itkDifferenceOfGaussiansGradientImageFilter.h"
#include "itkGradientMagnitudeImageFilter.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkSubtractImageFilter.h"
#include "itkRankImageFilter.h"
#include "AdaBoost.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkImageLinearConstIteratorWithIndex.h"
#include "itkImageLinearIteratorWithIndex.h"

typedef itk::Image<unsigned char, 3 > UCImageType;
typedef itk::Image<signed short, 3> ImageType;
typedef itk::Image<float, 3> FloatImageType;
typedef itk::MeanImageFunction< ImageType >   MeanImageFunctionType;
typedef itk::MedianImageFunction< ImageType > MedianImageFunctionType;
typedef itk::VarianceImageFunction<ImageType> VarianceImageFunctionType;
typedef itk::SumOfSquaresImageFunction<ImageType> SumOfSquaresImageFunctionType;
typedef itk::CentralDifferenceImageFunction<ImageType> CenteralDifferenceImageFunctionType;
typedef itk::DiscreteGaussianDerivativeImageFunction<ImageType> GaussianDerivativeImageFunctionType;
typedef itk::DiscreteGradientMagnitudeGaussianImageFunction<ImageType> GradientMagnitudeImageFunctionType;
typedef itk::DiscreteHessianGaussianImageFunction<ImageType> HessianGaussianImageFunctionType;
typedef itk::BoxMeanImageFilter<ImageType,FloatImageType> MeanType;
typedef itk::BoxSigmaImageFilter<ImageType,FloatImageType> SigmaType;
typedef itk::Box2SumImageFilter<FloatImageType, FloatImageType> SumType;
typedef itk::RankImageFilter<ImageType, FloatImageType > RankType;
typedef itk::ConstNeighborhoodIterator< FloatImageType > NeighborhoodIteratorType;


typedef itk::CastImageFilter<ImageType, FloatImageType> CastType;

typedef itk::SobelEdgeDetectionImageFilter<ImageType, FloatImageType> SobelType;
typedef itk::ZeroCrossingImageFilter<ImageType, FloatImageType> ZeroCrossingType;
typedef itk::LaplacianRecursiveGaussianImageFilter<FloatImageType,FloatImageType >    LoGFilter;
typedef itk::HessianRecursiveGaussianImageFilter< ImageType > HessianFilterType;
typedef itk::Hessian3DToVesselnessMeasureImageFilter< float > VesselnessMeasureFilterType;
typedef itk::DerivativeImageFilter< ImageType, FloatImageType >  derivativeFilterType;

typedef itk::DifferenceOfGaussiansGradientImageFilter<FloatImageType,float> DOGFilterType;
typedef itk::GradientMagnitudeImageFilter<ImageType,FloatImageType>  GradientMagnitudeFilterType;
typedef itk::GradientMagnitudeRecursiveGaussianImageFilter<ImageType,FloatImageType> GMRGFilterType;

typedef itk::ImageRegionConstIterator< UCImageType > UCImageIteratorType;
typedef itk::ImageRegionConstIterator< ImageType > ShortImageIteratorType;
typedef itk::ImageRegionConstIterator< FloatImageType > FloatImageIteratorType;
typedef itk::ImageRegionIterator< FloatImageType> FloatImageWriteIteratorType;

typedef itk::ImageRegionConstIteratorWithIndex< UCImageType > UCIndexIteratorType;
typedef itk::ImageRegionIteratorWithIndex< UCImageType > UCIndexWriteIteratorType;

typedef itk::ImageRegionConstIteratorWithIndex< ImageType   > ShortIndexIteratorType;
typedef itk::ImageRegionConstIteratorWithIndex< FloatImageType > FloatIndexIteratorType;



typedef itk::ImageFileWriter< FloatImageType > WriterType;
typedef itk::ImageFileReader< ImageType > ReaderType;
typedef itk::ImageFileReader< UCImageType > UCReaderType;

typedef itk::SubtractImageFilter<UCImageType,UCImageType,UCImageType> SubtractFilterType;

typedef itk::SignedMaurerDistanceMapImageFilter<UCImageType, FloatImageType>  DistanceMapFilterType;
typedef itk::ImageLinearIteratorWithIndex< FloatImageType > FloatLineIteratorType;
typedef itk::ImageLinearConstIteratorWithIndex< UCImageType > UCLineIteratorType;

template<class ImageType>
void WriteImage(typename ImageType::Pointer image,int no)
{
	typedef itk::ImageFileWriter<ImageType> WriterType;
	WriterType::Pointer writer = WriterType::New();
	char filename[256];
	sprintf(filename,"./feature_%03d.nii.gz",no);
	writer->SetInput( image);
	writer->SetFileName(filename);
	//try
	//{
	//	writer->Update();
	//}
	//catch(itk::ExceptionObject &e)
	//{
	//	std::cout<<e.GetDescription()<<std::endl;
	//	return ;
	//}
}

int main( int argc, char * argv[] )
{
	
	if(argc<6)
	{
		std::cout<<"Usage: "<<argv[0]<<" input3dImage livermask3d skin_air_mask3d learnoutput iternum"<<std::endl;
		return -1;
	}

	
/*	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(argv[1]);
	reader->Update();

	UCReaderType::Pointer livermaskreader = UCReaderType::New();

	livermaskreader->SetFileName(argv[2]);
	livermaskreader->Update();


	UCReaderType::Pointer airskinreader = UCReaderType::New();
	airskinreader->SetFileName(argv[3]);
	airskinreader->Update();
*/



	ifstream imfile;
	std::string imstr;
	std::string livermaskstr;
	std::string abdominalstr;

	imfile.open(argv[1]);
	int subjectN=0;
	string t;

	while (!imfile.eof())
	{
		getline(imfile,t);
		if (t.length()>0)
			subjectN++;
	}

	imfile.clear();
	imfile.seekg(0, ios::beg);

	cout<<"subject #: "<<subjectN<<endl;

	ifstream livermaskfile;
	livermaskfile.open(argv[2]);

	ifstream abdominalmaskfile;
	abdominalmaskfile.open(argv[3]);


	std::vector<FLOATTYPE> X;
	std::vector<char> Y;
	int NFeature;

	int n=0;

	ofstream AdaFeatureFile("./features.txt", ios::out);

	int featureno=0;

	while (!imfile.eof())
	{
		n++;
		

		getline(imfile,imstr);
		std::cout<<"image file: "<<imstr<<std::endl;
		getline(livermaskfile,livermaskstr);
		std::cout<<"liver mask seg: "<<livermaskstr<<std::endl;
		getline(abdominalmaskfile,abdominalstr);
		std::cout<<"abdominal mask: "<<abdominalstr<<std::endl;

		ReaderType::Pointer orgimage = ReaderType::New();
		orgimage->SetFileName(imstr.c_str());
		try
		{
			orgimage->Update();
		}
		catch(itk::ExceptionObject &e)
		{
			std::cout<<e.GetDescription()<<std::endl;
			return -1;
		}

		UCReaderType::Pointer liverSeg = UCReaderType::New();
		liverSeg->SetFileName(livermaskstr.c_str());
		try
		{
			liverSeg->Update();
		}
		catch(itk::ExceptionObject &e)
		{
			std::cout<<e.GetDescription()<<std::endl;
			return -1;
		}

		UCReaderType::Pointer abdominalmaskSeg = UCReaderType::New();
		abdominalmaskSeg->SetFileName(abdominalstr.c_str());
		try
		{
			abdominalmaskSeg->Update();
		}
		catch(itk::ExceptionObject &e)
		{
			std::cout<<e.GetDescription()<<std::endl;
			return -1;
		}

		ImageType::Pointer image = orgimage->GetOutput();

		HessianFilterType::Pointer hessianFilter = HessianFilterType::New();
		VesselnessMeasureFilterType::Pointer vesselnessFilter = VesselnessMeasureFilterType::New();

		CastType::Pointer caster = CastType::New();
		caster->SetInput(image);
		caster->Update();

		std::vector<FloatImageType::Pointer> imagefeatures;

		double sigma = 4.0;


		for(int ii=0;ii<3;ii++)
		{
			AdaFeatureFile<<++featureno<<" mean at neighborhood radius "<<ii+1<<" x -1"<<std::endl;
			AdaFeatureFile<<++featureno<<" mean at neighborhood radius "<<ii+1<<" x +1"<<std::endl;
			AdaFeatureFile<<++featureno<<" mean at neighborhood radius "<<ii+1<<" y -1"<<std::endl;
			AdaFeatureFile<<++featureno<<" mean at neighborhood radius "<<ii+1<<" y +1"<<std::endl;
			AdaFeatureFile<<++featureno<<" mean at neighborhood radius "<<ii+1<<" z -1"<<std::endl;
			AdaFeatureFile<<++featureno<<" mean at neighborhood radius "<<ii+1<<" z +1"<<std::endl;

			AdaFeatureFile<<++featureno<<" sigma at neighborhood radius "<<ii+1<<" x -1"<<std::endl;
			AdaFeatureFile<<++featureno<<" sigma at neighborhood radius "<<ii+1<<" x +1"<<std::endl;
			AdaFeatureFile<<++featureno<<" sigma at neighborhood radius "<<ii+1<<" y -1"<<std::endl;
			AdaFeatureFile<<++featureno<<" sigma at neighborhood radius "<<ii+1<<" y +1"<<std::endl;
			AdaFeatureFile<<++featureno<<" sigma at neighborhood radius "<<ii+1<<" z -1"<<std::endl;
			AdaFeatureFile<<++featureno<<" sigma at neighborhood radius "<<ii+1<<" z +1"<<std::endl;

			AdaFeatureFile<<++featureno<<" minimum at neighborhood radius "<<ii+1<<" x -1"<<std::endl;
			AdaFeatureFile<<++featureno<<" minimum at neighborhood radius "<<ii+1<<" x +1"<<std::endl;
			AdaFeatureFile<<++featureno<<" minimum at neighborhood radius "<<ii+1<<" y -1"<<std::endl;
			AdaFeatureFile<<++featureno<<" minimum at neighborhood radius "<<ii+1<<" y +1"<<std::endl;
			AdaFeatureFile<<++featureno<<" minimum at neighborhood radius "<<ii+1<<" z -1"<<std::endl;
			AdaFeatureFile<<++featureno<<" minimum at neighborhood radius "<<ii+1<<" z +1"<<std::endl;

			AdaFeatureFile<<++featureno<<" median at neighborhood radius "<<ii+1<<" x -1"<<std::endl;
			AdaFeatureFile<<++featureno<<" median at neighborhood radius "<<ii+1<<" x +1"<<std::endl;
			AdaFeatureFile<<++featureno<<" median at neighborhood radius "<<ii+1<<" y -1"<<std::endl;
			AdaFeatureFile<<++featureno<<" median at neighborhood radius "<<ii+1<<" y +1"<<std::endl;
			AdaFeatureFile<<++featureno<<" median at neighborhood radius "<<ii+1<<" z -1"<<std::endl;
			AdaFeatureFile<<++featureno<<" median at neighborhood radius "<<ii+1<<" z +1"<<std::endl;

			AdaFeatureFile<<++featureno<<" maximum at neighborhood radius "<<ii+1<<" x -1"<<std::endl;
			AdaFeatureFile<<++featureno<<" maximum at neighborhood radius "<<ii+1<<" x +1"<<std::endl;
			AdaFeatureFile<<++featureno<<" maximum at neighborhood radius "<<ii+1<<" y -1"<<std::endl;
			AdaFeatureFile<<++featureno<<" maximum at neighborhood radius "<<ii+1<<" y +1"<<std::endl;
			AdaFeatureFile<<++featureno<<" maximum at neighborhood radius "<<ii+1<<" z -1"<<std::endl;
			AdaFeatureFile<<++featureno<<" maximum at neighborhood radius "<<ii+1<<" z +1"<<std::endl;

		}



		imagefeatures.push_back(caster->GetOutput());  // intensity 
		WriteImage<FloatImageType>(caster->GetOutput(),featureno);
		AdaFeatureFile<<++featureno<<" intensity"<<std::endl;

		for(int ii=0;ii<3;ii++)
		{

			MeanType::Pointer meanfilter = MeanType::New();
			SigmaType::Pointer sigmafilter = SigmaType::New();
			RankType::Pointer rankfilter0 = RankType::New(); // minimum
			RankType::Pointer rankfilter1 = RankType::New(); // median
			RankType::Pointer rankfilter2 = RankType::New(); // maximum

			MeanType::RadiusType radius;
			radius.Fill( ii+1 );

			meanfilter->SetInput(image);
			meanfilter->SetRadius(radius);
			meanfilter->Update();

			imagefeatures.push_back(meanfilter->GetOutput());
			WriteImage<FloatImageType>(meanfilter->GetOutput(),featureno);

			AdaFeatureFile<<++featureno<<" mean at radius "<<ii+1<<std::endl;


			sigmafilter->SetInput(image);
			sigmafilter->SetRadius(radius);
			sigmafilter->Update();

			imagefeatures.push_back(sigmafilter->GetOutput());
			WriteImage<FloatImageType>(sigmafilter->GetOutput(),featureno);

			AdaFeatureFile<<++featureno<<" sigma at radius "<<ii+1<<std::endl;

		
			rankfilter0->SetInput(image);
			rankfilter0->SetRadius(radius);
			rankfilter0->SetRank(0);    // minimum value
			rankfilter0->Update();
			imagefeatures.push_back(rankfilter0->GetOutput());
			WriteImage<FloatImageType>(rankfilter0->GetOutput(),featureno);
			AdaFeatureFile<<++featureno<<" minimum at radius "<<ii+1<<std::endl;
		

			rankfilter1->SetInput(image);
			rankfilter1->SetRadius(radius);
			rankfilter1->SetRank(0.5);  // median value
			rankfilter1->Update();
			imagefeatures.push_back(rankfilter1->GetOutput());
			WriteImage<FloatImageType>(rankfilter1->GetOutput(),featureno);
			AdaFeatureFile<<++featureno<<" median at radius "<<ii+1<<std::endl;
			

			rankfilter2->SetInput(image);
			rankfilter2->SetRadius(radius);
			rankfilter2->SetRank(1);    // maximum value
			rankfilter2->Update();
			imagefeatures.push_back(rankfilter2->GetOutput());
			WriteImage<FloatImageType>(rankfilter2->GetOutput(),featureno);

			AdaFeatureFile<<++featureno<<" maximum at radius "<<ii+1<<std::endl;

			
		}

	

		//SumType::Pointer sumfilter = SumType::New();
//		sumfilter->SetInput(caster->GetOutput());
//		sumfilter->Update();

		hessianFilter->SetInput( image );
		hessianFilter->SetSigma( sigma );


		// Create and setup a derivative filter
		derivativeFilterType::Pointer derivativeFilter_x = derivativeFilterType::New();
		derivativeFilterType::Pointer derivativeFilter_y = derivativeFilterType::New();
		derivativeFilterType::Pointer derivativeFilter_z = derivativeFilterType::New();

		derivativeFilter_x->SetInput( image );
		derivativeFilter_y->SetInput( image );
		derivativeFilter_z->SetInput( image );

		derivativeFilter_x->SetDirection(0); // "x" axis
		derivativeFilter_y->SetDirection(1); // "y" axis
		derivativeFilter_z->SetDirection(2); // "z" axis

		derivativeFilter_x->Update();
		derivativeFilter_y->Update();
		derivativeFilter_z->Update();

		imagefeatures.push_back(derivativeFilter_x->GetOutput());
		WriteImage<FloatImageType>(derivativeFilter_x->GetOutput(),featureno);
		AdaFeatureFile<<++featureno<<" derivative x "<<std::endl;
		imagefeatures.push_back(derivativeFilter_y->GetOutput());
		WriteImage<FloatImageType>(derivativeFilter_y->GetOutput(),featureno);
		AdaFeatureFile<<++featureno<<" derivative y "<<std::endl;
		imagefeatures.push_back(derivativeFilter_z->GetOutput());
		WriteImage<FloatImageType>(derivativeFilter_z->GetOutput(),featureno);
		AdaFeatureFile<<++featureno<<" derivative z "<<std::endl;

		DOGFilterType::Pointer DOGFilter = DOGFilterType::New();
		DOGFilter->SetInput(caster->GetOutput());
		DOGFilter->SetWidth(4);
		DOGFilter->Update();

		//imagefeatures.push_back(DOGFilter->GetOutput());

		DOGFilterType::TOutputImage::Pointer gradientImage = DOGFilter->GetOutput();
		GradientMagnitudeFilterType::Pointer gradientmagnitudefilter = GradientMagnitudeFilterType::New();

		gradientmagnitudefilter->SetInput(image);
		gradientmagnitudefilter->Update();
		imagefeatures.push_back(gradientmagnitudefilter->GetOutput());
		WriteImage<FloatImageType>(gradientmagnitudefilter->GetOutput(),featureno);
		AdaFeatureFile<<++featureno<<" gradient magnitude"<<std::endl;

		GMRGFilterType::Pointer gmrgfilter = GMRGFilterType::New();
		gmrgfilter->SetInput(image);
		gmrgfilter->Update();

		imagefeatures.push_back(gmrgfilter->GetOutput());
		WriteImage<FloatImageType>(gmrgfilter->GetOutput(),featureno);
		AdaFeatureFile<<++featureno<<" gradient magnitude recursive gaussian"<<std::endl;
	
		UCIndexIteratorType livermaskit(liverSeg->GetOutput(),liverSeg->GetOutput()->GetLargestPossibleRegion());
		
		UCIndexIteratorType air_skin_it(abdominalmaskSeg->GetOutput(),abdominalmaskSeg->GetOutput()->GetLargestPossibleRegion());

		ShortIndexIteratorType imageit(image,image->GetLargestPossibleRegion());


		UCImageType::Pointer im = UCImageType::New();

		im->SetRegions(image->GetLargestPossibleRegion());
		im->CopyInformation(image);
		im->Allocate();
		im->FillBuffer(1);

		SubtractFilterType::Pointer subtracter = SubtractFilterType::New();

		subtracter->SetInput(0,im);
		subtracter->SetInput(1,abdominalmaskSeg->GetOutput());
		subtracter->Update();

		UCIndexWriteIteratorType abdominalmask_it(subtracter->GetOutput(),subtracter->GetOutput()->GetLargestPossibleRegion());


		abdominalmask_it.GoToBegin();
		livermaskit.GoToBegin();
		air_skin_it.GoToBegin();
		imageit.GoToBegin();


		// location features
		FloatImageType::Pointer outputImageX = FloatImageType::New();
		FloatImageType::Pointer outputImageY = FloatImageType::New();
		FloatImageType::Pointer outputImageXdist = FloatImageType::New();
		FloatImageType::Pointer outputImageYdist = FloatImageType::New();
		outputImageX->SetRegions( image->GetLargestPossibleRegion() );
		outputImageY->SetRegions( image->GetLargestPossibleRegion() );
		outputImageXdist->SetRegions( image->GetLargestPossibleRegion() );
		outputImageYdist->SetRegions( image->GetLargestPossibleRegion() );

		outputImageX->CopyInformation( image );
		outputImageY->CopyInformation( image );
		outputImageXdist->CopyInformation( image );
		outputImageYdist->CopyInformation( image );
		outputImageX->Allocate();
		outputImageY->Allocate();
		outputImageXdist->Allocate();
		outputImageYdist->Allocate();

		outputImageXdist->FillBuffer(0);
		outputImageYdist->FillBuffer(0);


		FloatImageWriteIteratorType outputItX( outputImageX, outputImageX->GetLargestPossibleRegion() );
		FloatImageWriteIteratorType outputItY( outputImageY, outputImageY->GetLargestPossibleRegion() );
		FloatLineIteratorType outputItXdist( outputImageXdist, outputImageXdist->GetLargestPossibleRegion() );
		FloatLineIteratorType outputItYdist( outputImageYdist, outputImageYdist->GetLargestPossibleRegion() );

		FloatImageType::IndexType requestedIndexX =outputImageX->GetLargestPossibleRegion().GetIndex();
		FloatImageType::IndexType requestedIndexY =outputImageY->GetLargestPossibleRegion().GetIndex();
		FloatImageType::IndexType requestedIndexXdist = outputImageXdist->GetLargestPossibleRegion().GetIndex();
		FloatImageType::IndexType requestedIndexYdist = outputImageYdist->GetLargestPossibleRegion().GetIndex();

		FloatImageType::SizeType requestedSize = outputImageX->GetLargestPossibleRegion().GetSize();
		for ( outputItX.GoToBegin(),outputItY.GoToBegin() ;
			!outputItX.IsAtEnd(); ++outputItX,++outputItY)
		{
			FloatImageType::IndexType idx = outputItX.GetIndex();
			outputItX.Set(idx[0]);
			outputItY.Set(idx[1]);
		}

		//imagefeatures.push_back(outputImageX);
		WriteImage<FloatImageType>(outputImageX,featureno);
		AdaFeatureFile<<++featureno<<" location x coordinate"<<std::endl;
		imagefeatures.push_back(outputImageY);
		WriteImage<FloatImageType>(outputImageY,featureno);
		AdaFeatureFile<<++featureno<<" location y coordinate"<<std::endl;
	//	imagefeatures.push_back(outputImageZ);
	//	AdaFeatureFile<<++featureno<<" location z coordinate"<<std::endl;

		DistanceMapFilterType::Pointer distancemapfilter = DistanceMapFilterType::New();
		distancemapfilter->SetInput( abdominalmaskSeg->GetOutput() );
		distancemapfilter->SetSquaredDistance( false );
		distancemapfilter->SetUseImageSpacing( true  );
		distancemapfilter->SetInsideIsPositive( true );
		distancemapfilter->Update();


		imagefeatures.push_back(distancemapfilter->GetOutput());
		WriteImage<FloatImageType>(distancemapfilter->GetOutput(),featureno);
		AdaFeatureFile<<++featureno<<" distance map"<<std::endl;


		UCLineIteratorType abdominallineit(subtracter->GetOutput(),subtracter->GetOutput()->GetLargestPossibleRegion());
	
		abdominallineit.SetDirection(0);
		outputItXdist.SetDirection(0);
		abdominallineit.GoToBegin();
		outputItXdist.GoToBegin();

		UCImageType::SpacingType spacing = subtracter->GetOutput()->GetSpacing();

		while(!abdominallineit.IsAtEnd())
		{
			int abdominalleft = 9999;
			abdominallineit.GoToBeginOfLine();
			
			if(abdominallineit.Get())
				abdominalleft = 0;
			else
				while(!abdominallineit.IsAtEndOfLine())
				{
					if(abdominallineit.Get())
					{
						abdominalleft = abdominallineit.GetIndex()[0];
						break;
					}
					++abdominallineit;
				}
		

			outputItXdist.GoToBeginOfLine();

			if(abdominalleft!=9999)
			while(!outputItXdist.IsAtEndOfLine())
			{

			    float tmpdist = (outputItXdist.GetIndex()[0]-abdominalleft)*spacing[0];
			    outputItXdist.Set(tmpdist);
				++outputItXdist;
			}

			abdominallineit.NextLine();
			outputItXdist.NextLine();
		}
		
		// Y dist

		abdominallineit.SetDirection(1);
		outputItYdist.SetDirection(1);
		abdominallineit.GoToBegin();
		outputItYdist.GoToBegin();
		
		while(!abdominallineit.IsAtEnd())
		{
			int abdominaltop = 9999;
			abdominallineit.GoToBeginOfLine();

			if(abdominallineit.Get())
				abdominaltop = 0;
			else
				while(!abdominallineit.IsAtEndOfLine())
				{
					if(abdominallineit.Get())
					{
						abdominaltop = abdominallineit.GetIndex()[1];
						break;
					}
					++abdominallineit;
				}


				outputItYdist.GoToBeginOfLine();
				if(abdominaltop!=9999)
				while(!outputItYdist.IsAtEndOfLine())
				{

					float tmpdist = (outputItYdist.GetIndex()[1]-abdominaltop)*spacing[1];
					outputItYdist.Set(tmpdist);
					++outputItYdist;
				}

				abdominallineit.NextLine();
				outputItYdist.NextLine();
		}


		imagefeatures.push_back(outputImageXdist);
		WriteImage<FloatImageType>(outputImageXdist,featureno);

		AdaFeatureFile<<++featureno<<" distance of X"<<std::endl;
		
		imagefeatures.push_back(outputImageYdist);
		WriteImage<FloatImageType>(outputImageYdist,featureno);

		AdaFeatureFile<<++featureno<<" distance of Y"<<std::endl;

		

		std::vector<FloatImageIteratorType> imagefeaturesiterator;

	    NFeature = imagefeatures.size();

		for(int i=0;i<NFeature;i++)
		{
			FloatImageIteratorType it(imagefeatures[i],imagefeatures[i]->GetRequestedRegion());
			it.GoToBegin();
			imagefeaturesiterator.push_back(it);
		}

		std::cout<<"adding context features..."<<std::endl;

		NeighborhoodIteratorType::RadiusType neighborhood_radius;
		neighborhood_radius.Fill(3);

		//NeighborhoodIteratorType it(neighborhood_radius,mean,mean);

		std::vector<NeighborhoodIteratorType> mean_context_iterators;
		std::vector<NeighborhoodIteratorType> sigma_context_iterators;
		std::vector<NeighborhoodIteratorType> minimum_context_iterators;
		std::vector<NeighborhoodIteratorType> maximum_context_iterators;
		std::vector<NeighborhoodIteratorType> median_context_iterators;

		for(int i=0;i<3;i++)
		{
			int t = i*5;
			NeighborhoodIteratorType it1(neighborhood_radius,imagefeatures[t],imagefeatures[t]->GetLargestPossibleRegion());
			mean_context_iterators.push_back(it1);

			NeighborhoodIteratorType it2(neighborhood_radius,imagefeatures[t+1],imagefeatures[t+1]->GetLargestPossibleRegion());
			sigma_context_iterators.push_back(it2);

			NeighborhoodIteratorType it3(neighborhood_radius,imagefeatures[t+2],imagefeatures[t+2]->GetLargestPossibleRegion());
			minimum_context_iterators.push_back(it3);

			NeighborhoodIteratorType it4(neighborhood_radius,imagefeatures[t+3],imagefeatures[t+3]->GetLargestPossibleRegion());
			median_context_iterators.push_back(it4);

			NeighborhoodIteratorType it5(neighborhood_radius,imagefeatures[t+4],imagefeatures[t+4]->GetLargestPossibleRegion());
			maximum_context_iterators.push_back(it5);

		}
		NeighborhoodIteratorType::OffsetType offset1 = { {-3, 0, 0} };
		NeighborhoodIteratorType::OffsetType offset2 = { { 3, 0, 0} };
		NeighborhoodIteratorType::OffsetType offset3 = { { 0,-3, 0} };
		NeighborhoodIteratorType::OffsetType offset4 = { { 0, 3, 0} };
		NeighborhoodIteratorType::OffsetType offset5 = { { 0, 0,-3} };
		NeighborhoodIteratorType::OffsetType offset6 = { { 0, 0, 3} };


		while(!imageit.IsAtEnd())
		{

			if(abdominalmask_it.Get())
			{
				for(int i=0;i<3;i++)
				{
					X.push_back(mean_context_iterators[i].GetPixel(offset1));
					X.push_back(mean_context_iterators[i].GetPixel(offset2));
					X.push_back(mean_context_iterators[i].GetPixel(offset3));
					X.push_back(mean_context_iterators[i].GetPixel(offset4));
					X.push_back(mean_context_iterators[i].GetPixel(offset5));
					X.push_back(mean_context_iterators[i].GetPixel(offset6));
					
					X.push_back(sigma_context_iterators[i].GetPixel(offset1));
					X.push_back(sigma_context_iterators[i].GetPixel(offset2));
					X.push_back(sigma_context_iterators[i].GetPixel(offset3));
					X.push_back(sigma_context_iterators[i].GetPixel(offset4));
					X.push_back(sigma_context_iterators[i].GetPixel(offset5));
					X.push_back(sigma_context_iterators[i].GetPixel(offset6));
					
					X.push_back(minimum_context_iterators[i].GetPixel(offset1));
					X.push_back(minimum_context_iterators[i].GetPixel(offset2));
					X.push_back(minimum_context_iterators[i].GetPixel(offset3));
					X.push_back(minimum_context_iterators[i].GetPixel(offset4));
					X.push_back(minimum_context_iterators[i].GetPixel(offset5));
					X.push_back(minimum_context_iterators[i].GetPixel(offset6));
					
					X.push_back(median_context_iterators[i].GetPixel(offset1));
					X.push_back(median_context_iterators[i].GetPixel(offset2));
					X.push_back(median_context_iterators[i].GetPixel(offset3));
					X.push_back(median_context_iterators[i].GetPixel(offset4));
					X.push_back(median_context_iterators[i].GetPixel(offset5));
					X.push_back(median_context_iterators[i].GetPixel(offset6));

					
					X.push_back(maximum_context_iterators[i].GetPixel(offset1));
					X.push_back(maximum_context_iterators[i].GetPixel(offset2));
					X.push_back(maximum_context_iterators[i].GetPixel(offset3));
					X.push_back(maximum_context_iterators[i].GetPixel(offset4));
					X.push_back(maximum_context_iterators[i].GetPixel(offset5));
					X.push_back(maximum_context_iterators[i].GetPixel(offset6));

				}

				for(int i=0;i<NFeature;i++)
				{
					X.push_back(imagefeaturesiterator[i].Get());
				}


				if(livermaskit.Get())
					Y.push_back(1);
				else
					Y.push_back(-1);

			}

			for(int i=0;i<NFeature;i++)
				++imagefeaturesiterator[i];

			for(int i=0;i<3;i++)
			{
				++mean_context_iterators[i];
				++sigma_context_iterators[i];
				++minimum_context_iterators[i];
				++maximum_context_iterators[i];
				++median_context_iterators[i];
			}

			++imageit;
			++air_skin_it;
			++livermaskit;
			++abdominalmask_it;
		}
	}

	std::cout<<"Features finished..."<<std::endl;
	std::cout<<X.size()<<" "<<Y.size()<<" "<<X.size()/Y.size()<<std::endl;
	std::cout<<"Number of feature: "<<NFeature<<std::endl;

	char fileName[1024];
	sprintf (fileName, "%s-AdaBoostResult",argv[4]);

	try
	{
	AdaBoostTrain(&X[0],&Y[0],Y.size(),X.size()/Y.size(),atoi(argv[5]),fileName);
	}
	catch(const exception &e)
	{
		std::cout<<e.what()<<std::endl;
		return -1;
	}


	AdaFeatureFile.close();


	return EXIT_SUCCESS;
}