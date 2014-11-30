#include <iostream>
#include <sstream>

#include "itkMeshFileWriter.h"
#include "itkImageFileWriter.h"

#include "LiverSegmentationAPI.h"

namespace km
{
	template<class TPointSet>
	class NotifierImp:public NotifierBase<TPointSet>
	{
	public:
		typedef TPointSet PointSetType;

		NotifierImp()
		{
			iter_num = 0;
		}

		~NotifierImp()
		{

		}

		void setOutputDir(const char* dir)
		{
			strcpy (outputdir,dir);
		}

		void notify( const PointSetType * pointSet )
		{
			/*
			std::stringstream ss;
			if (outputdir != NULL)
			{
				ss  <<  outputdir  <<  "/"  <<  "meshUpdated"   <<  "."  <<  iter_num++  <<  ".vtk";
			}
			else
			{
				ss  <<  "meshUpdated"   <<  "."  <<  iter_num++  <<  ".vtk";
			}
			typedef itk::MeshFileWriter<PointSetType> MeshWriterType;
			MeshWriterType::Pointer writer = MeshWriterType::New();
			writer->SetInput( pointSet );
			writer->SetFileName( ss.str().c_str() );
			try
			{
				writer->Update();
			}
			catch( itk::ExceptionObject & err )
			{
				std::cerr << "ExceptionObject caught !" << std::endl;
				std::cerr << err << std::endl;
			}
			catch(...)
			{
				std::cerr << "Shit happens!"<<std::endl;
			}
			*/
		}

	private:
		int iter_num;
		char outputdir[1024];
	};
}

int main(int argc, char* argv[])
{
	if(argc<9)
	{
		std::cerr<<"Usage: "<<std::endl;
		std::cout<<" inputImage"
		         <<" SSMFile"
                 <<" outputDir"
				 <<" boundaryClassifierFile"
				 <<" liverClassifierFile"
				 <<" adaboostSegmentFile"
				 <<" geoFile"
				 <<" atlasImageFile"
				 <<" configFile"
				 <<std::endl;
		return -1;
	}

	int paramidx = 1;

	const char* inputImageFile = argv[paramidx++];
	const char* SSMFile = argv[paramidx++];
	const char* outputdir = argv[paramidx++];
	const char* boundaryClassifierFile = argv[paramidx++];
	const char* liverClassifierFile = argv[paramidx++];
	const char* adaboostSegmentFile = argv[paramidx++];
	const char* geoFile = argv[paramidx++];
	const char* atlasImageFile = argv[paramidx++];
	const char* configFile = argv[paramidx++];

	std::cout<<"** input image                 : " << inputImageFile << std::endl;
	std::cout<<"** SSM file                    : " << SSMFile << std::endl;
	std::cout<<"** output dir                  : " << outputdir << std::endl;
	std::cout<<"** boundary classifier file    : " << boundaryClassifierFile << std::endl;
	std::cout<<"** liver classifier file       : " << liverClassifierFile << std::endl;
	std::cout<<"** adaboost segmentation file  : " << adaboostSegmentFile << std::endl;
	std::cout<<"** geometry file               : " << geoFile << std::endl;
	std::cout<<"** atlas image                 : " << atlasImageFile << std::endl;
	std::cout<<"** config file                 : " <<configFile<<std::endl;

	typedef km::MeshType MeshType;
	typedef km::NotifierImp<MeshType> NotifierType;
	NotifierType* notifier = new NotifierType;
	notifier->setOutputDir(outputdir);

	typedef km::UCharImageType BinaryImageType;
	BinaryImageType::Pointer segmentationResult = BinaryImageType::New();

	std::cout<<"Start to call LiverSeg(..)"<<std::endl;

	LiverSeg( segmentationResult,
			  notifier,
			  outputdir,
	          inputImageFile,
			  SSMFile,
			  boundaryClassifierFile,
			  liverClassifierFile,
			  adaboostSegmentFile,
			  geoFile,
			  atlasImageFile,
			  configFile);

	/*
	std::stringstream ss;
	ss << outputdir << "/" <<"segmentationResult.nii.gz";

	typedef itk::ImageFileWriter<BinaryImageType> ImageFileWriterType;
	typedef ImageFileWriterType::Pointer ImageFileWriterPointer;
	ImageFileWriterPointer writer = ImageFileWriterType::New();
	writer->SetFileName( ss.str().c_str() );
	writer->SetInput( segmentationResult );
	try
	{
		writer->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
	}
	*/

	return 0;
}