#include "vtkMaskPolyData.h"
#include <vtkSurfaceReconstructionFilter.h>
#include <vtkContourFilter.h>
#include <vtkReverseSense.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkDecimatePro.h>

#include "vtkTriangleFilter.h"
#include "vtkQuadricClustering.h"
#include <vtkFillHolesFilter.h>

#include "itkTimeProbesCollectorBase.h"
#include "itkMemoryProbesCollectorBase.h"

#include "kmUtility.h"

using namespace std;
using namespace km;

int main(int argc, char* argv[]) {

	if (argc < 4) {
		std::cout << "usage " << argv[0] << " inputVtkPolyData outputVtkPolyData -ratio newRatio -number number [-smooth iterations]" << std::endl;
		exit(-1);
	}

	itk::TimeProbesCollectorBase chronometer;
	itk::MemoryProbesCollectorBase memorymeter;
	memorymeter.Start( "converting" );
	chronometer.Start( "converting" );

	bool ratioFlag = true;
	bool smoothFlag = false;
	double ratio = 1.0;
	unsigned int pointsNumber = 0;
	unsigned int smoothIterations = 0;

	for(int i=1;i<argc;i++)
	{
		if(strcmp(argv[i], "-ratio") == 0)
		{
			ratioFlag = true;
			ratio=atof(argv[++i]);
		}
		else if(strcmp(argv[i], "-number") == 0)
		{
			ratioFlag = false;
			pointsNumber = atoi(argv[++i]);
		}
		else if(strcmp(argv[i], "-smooth") == 0)
		{
			smoothFlag = true;
			smoothIterations = atoi(argv[++i]);
		}
	}

	vtkSmartPointer<vtkPolyData> inputPolyData = readPolyData( argv[1] );


	//vtkTriangleFilter *tri = vtkTriangleFilter::New();
	//tri->SetInput( inputPolyData );

	//vtkQuadricClustering *decimate = vtkQuadricClustering::New();
	////decimate->SetNumberOfXDivisions(32);
	////decimate->SetNumberOfYDivisions(32);
	////decimate->SetNumberOfZDivisions(32);
	//decimate->SetAutoAdjustNumberOfDivisions(1);
	//decimate->SetInput(tri->GetOutput());

	//decimate->Update();

	////tri->Update();

	//writePolyData( argv[2],decimate->GetOutput() );

	//return 0;

	vtkSmartPointer<vtkPolyData> outputPolyData = vtkSmartPointer<vtkPolyData>::New();

	if(!ratioFlag)
	{
		int origPointsNumber = inputPolyData->GetNumberOfPoints();
		std::cout<<"original points number: "<<origPointsNumber<<std::endl;

		ratio = static_cast<double>(pointsNumber) / origPointsNumber;
	}

	outputPolyData = km::decimatePolydata(inputPolyData, ratio);

	if( smoothFlag )
	{
		std::cout<<"smooth iterations: "<<smoothIterations<<std::endl;
		outputPolyData = km::smoothPolyData(outputPolyData, smoothIterations);
	}

	//vtkSmartPointer<vtkFillHolesFilter> holeFiller = vtkSmartPointer<vtkFillHolesFilter>::New();
	//holeFiller->SetInput( outputPolyData );
	//holeFiller->Update();

	//vtkSmartPointer<vtkPolyData> filledPolydata = holeFiller->GetOutput();
	//writePolyData( "filledPolydata.vtk",outputPolyData );

	writePolyData( argv[2],outputPolyData );

	chronometer.Stop( "converting" );
	memorymeter.Stop( "converting" );
	chronometer.Report( std::cout );
	memorymeter.Report( std::cout );

	return EXIT_SUCCESS;
}