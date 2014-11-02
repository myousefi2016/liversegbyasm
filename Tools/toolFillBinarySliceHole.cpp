#include "itkTimeProbesCollectorBase.h"
#include "itkMemoryProbesCollectorBase.h"

#include "kmUtility.h"
#include "kmProcessing.h"

using namespace std;
using namespace km;

const unsigned Dimensions = 3;
typedef itk::Image<unsigned char, Dimensions> BinaryImageType;

int main(int argc, char* argv[]) {

  if (argc < 3) {
		std::cout << "usage " << argv[0] << " inputBinary outputBinary" << std::endl;
    exit(-1);
  }

  itk::TimeProbesCollectorBase chronometer;
  itk::MemoryProbesCollectorBase memorymeter;
  memorymeter.Start( "filling" );
  chronometer.Start( "filling" );

	int type = 0;
	if (argc>3)
	{
		type = atoi(argv[3]);
	}

	BinaryImageType::Pointer inputBinaryImage = km::readImage<BinaryImageType>( argv[1] );

	km::fillSliceHole<BinaryImageType>( inputBinaryImage );
	
	km::writeImage<BinaryImageType>( argv[2], inputBinaryImage );

  chronometer.Stop( "filling" );
  memorymeter.Stop( "filling" );
  chronometer.Report( std::cout );
  memorymeter.Report( std::cout );

  return EXIT_SUCCESS;
}

