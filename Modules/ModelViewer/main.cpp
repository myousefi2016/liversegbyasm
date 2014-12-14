
#include "QVTKApplication.h"
#include "ModelViewer.h"



using namespace statismo;


int main(int argc, char** argv)
{
  QVTKApplication app(argc, argv);
  ModelViewer widget;
  
  widget.show();

  return app.exec();
//  MeshStatisticalModelPointer meshModel = MeshStatisticalModelType::New();
//  meshModel->Load( "k:/model.h5" );
//  std::cout << "loaded model with " << meshModel->GetNumberOfPrincipalComponents() << " Principal Components" << std::endl;
//  
//  system("pause");
 
 //return 0;
}
