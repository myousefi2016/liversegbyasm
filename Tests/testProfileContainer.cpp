
#include "kmProfileContainer.h"

int main(int argc, char* argv[])
{
	const char* profileFilename = "D:\\Workspace\\ASM\\projects\\LiverSegbyASM\\experiments\\training\\appearance\\output_20140815\\samples_PLAIN.txt";
	const char* hdfFilename = "Samples.h5";

	//std::map<int, km::ProfileContainer*> tmpmap;
	//if (tmpmap[999] == NULL)
	//{
	//	std::cout<<"No found."<<std::endl;
	//}
	//else
	//{
	//	std::cout<<"Found."<<std::endl;
	//}

	km::ProfileContainer profileContainer;
	profileContainer.setShapeNumber(10);
	profileContainer.setShapePointsNumber(3996);
	profileContainer.loadSamples(profileFilename);
	profileContainer.print();
	profileContainer.save(hdfFilename);
	profileContainer.cluster(9);
	profileContainer.print();
	profileContainer.save("Samples_Clustered.h5");

	km::ProfileContainer profileContainerCopied;
	profileContainerCopied.load(hdfFilename);
	profileContainerCopied.copyCluster(profileContainer);
	profileContainerCopied.save("Samples_Copied.h5");

	return 0;
}