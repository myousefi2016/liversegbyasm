#ifndef __kmProfileContainer_h
#define __kmProfileContainer_h

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <math.h>
#include <set>
#include <fstream>
#include <map>
#include <set>

#include <itkVariableLengthVector.h>
#include "itkListSample.h"
#include "itkKdTree.h"
#include "itkWeightedCentroidKdTreeGenerator.h"
#include "itkMinimumDecisionRule.h"
#include "itkSampleClassifierFilter.h"
#include "itkKdTreeBasedKmeansEstimator.h"
#include "itkVariableLengthVector.h"

#include "statismo/CommonTypes.h"
#include "statismo/HDF5Utils.h"

#include "kmGlobal.h"

using namespace std;
using namespace statismo;

namespace km
{
	typedef unsigned long long LONGTYPE;
	/************************************************************************************************/
	// Profile container
	/************************************************************************************************/	
	class ProfileContainer
	{
	public:
		int numberOfShapes;
		int numberOfPointsEachShape;
		int numberOfSamplesPerLandmark;
		int profileDimension;

		MatrixType profiles;
		VectorType classLabels;
		VectorType pointIDs;
		VectorType clusterLabels;
		std::set<int> clusterLabelSet;

		enum SampleOffset
		{
			Offset_PointId = 0,
			Offset_ClassLabel,
			Offset_Profile
		};
	private:
		bool initilized;
		static const int defaultClusterLabel = 0;
	public:
		ProfileContainer( )
		{
			initilized = false;
			numberOfShapes = 0;
			numberOfPointsEachShape = 0;
			numberOfSamplesPerLandmark = 1;
			profileDimension = 0;
		}

		~ProfileContainer()
		{
		}

		void setShapeNumber(int n)
		{
			this->numberOfShapes = n;	
		}

		void setShapePointsNumber(int n)
		{
			this->numberOfPointsEachShape = n;
		}

		void push_back( int ptid, int label, VectorType & single_profile, const LONGTYPE rowid )
		{
			this->pointIDs[rowid] = ptid;
			this->classLabels[rowid] = label;
			for (int dim=0;dim<this->profileDimension;dim++){
				this->profiles(rowid, dim) = single_profile[dim];
			}
		}

		void updateCluseterLabelSet()
		{
			clusterLabelSet.clear();
			for (int i=0;i<clusterLabels.size();i++)
			{
				clusterLabelSet.insert(clusterLabels[i]);
			}
		}

		int getClusterLabelByRowId(LONGTYPE rowid)
		{
			//LONGTYPE shapeBlockSize = numberOfPointsEachShape*numberOfSamplesPerLandmark;
			//int pointBlockSize = numberOfSamplesPerLandmark;
			//int pointId = (rowid%shapeBlockSize)/pointBlockSize;
			//return this->clusterLabels[pointId];
			int pointId = this->pointIDs[rowid];
			return this->clusterLabels[pointId];
		}

		void allocate(const char* sampleFilename)
		{
			LONGTYPE numberOfRows = 0;
			int      numberOfColumns = 0;
			//Read profile count and profile dimension
			
			string line;
			ifstream myfile (sampleFilename);
			if (myfile.is_open()){
				while ( getline (myfile,line) ){
					if ( line.find_first_not_of(' ') != std::string::npos ){
						if (numberOfRows == 0){
							istringstream s;
							s.str(line);
							double t;
							while (s >> t){
								if (!s.good()){
									break;
								}else{
									numberOfColumns ++;
								}
							}
						}
						numberOfRows++;
					}
				}
				myfile.close();
			}else{
				std::cout << "Unable to open file: " << sampleFilename << std::endl;
				return;
			}
			this->profileDimension = numberOfColumns - Offset_Profile;
			this->numberOfSamplesPerLandmark = (numberOfRows/this->numberOfShapes)/this->numberOfPointsEachShape;
			this->clusterLabels.setConstant(this->numberOfPointsEachShape, defaultClusterLabel);
			this->profiles.setZero(numberOfRows, this->profileDimension);
			this->classLabels.setZero(numberOfRows);
			this->pointIDs.setZero(numberOfRows);
			std::cout<<"Profile dimension: "<<this->profileDimension<<std::endl;
			std::cout<<"Profile count: "<<numberOfRows<<std::endl;
			std::cout<<"Number of samples per landmark: "<<this->numberOfSamplesPerLandmark<<std::endl;
		}

		void loadSamples( const char* filename )
		{
			std::cout<<"Load from sample file: "<<filename<<std::endl;
			
			this->allocate(filename);

			int ptid;
			int label;
			VectorType profile_single = VectorType::Zero(this->profileDimension);

			//Read data and label
			string line;
			ifstream myfile (filename);

			LONGTYPE row_num = 0;  //useless
			int      col_offset = 0;

			if (myfile.is_open()){
				while ( getline (myfile,line) ){
					if ( line.find_first_not_of(' ') != std::string::npos ){
						istringstream s;
						s.str(line);
						float t = 0.0;
						col_offset = 0;
						while (s >> t)
						{
							if (!s.good()){
								break;
							}else{
								if (col_offset == Offset_PointId)
								{
									ptid = int(t);
								}
								else 
								if (col_offset == Offset_ClassLabel){
									label = int(t);
								}else{
									if (col_offset>=this->profileDimension+Offset_Profile){
										std::cout<<"Wrong data line: "<<col_offset<<" "<<row_num<<" "<<t<<std::endl;
									}else{
										profile_single[col_offset-Offset_Profile] = t;
									}
								}
								col_offset++;
							}
						}
						this->push_back(ptid, label, profile_single, row_num++);
					}
				}
				myfile.close();
			}
			else
			{
				std::cout << "Unable to open file: " << filename << std::endl;
				return;
			}

			this->updateCluseterLabelSet();
			std::cout<<"Load from sample file done"<<std::endl;
		}

		void save( const std::string filename )
		{
			using namespace H5;
			H5File file;
			try {
				file = H5::H5File( filename.c_str(), H5F_ACC_TRUNC);
			} catch (FileIException& e) {
				std::string msg(std::string("Could not open HDF5 file for writing \n") + e.getCDetailMsg());
				throw StatisticalModelException(msg.c_str());
			}
			H5::Group groupRoot = file.openGroup("/");
			try {
				// create the group structure
				HDF5Utils::writeInt(groupRoot, "numberOfPointsEachShape", this->numberOfPointsEachShape);
				HDF5Utils::writeInt(groupRoot, "numberOfShapes", this->numberOfShapes);
				HDF5Utils::writeInt(groupRoot, "profileDimension", this->profileDimension);
				HDF5Utils::writeMatrix(groupRoot, "profiles", this->profiles);
				HDF5Utils::writeVector(groupRoot, "classLabels", this->classLabels);
				HDF5Utils::writeVector(groupRoot, "pointIDs", this->pointIDs);
				HDF5Utils::writeVector(groupRoot, "clusterLabels", this->clusterLabels);
			} catch (H5::Exception& e) {
				std::string msg(std::string("an exception occurred while writing HDF5 file \n") + e.getCDetailMsg());
				std::cerr<<msg<<std::endl;
				throw StatisticalModelException(msg.c_str());
			}
			groupRoot.close();
			file.close();	

			std::cout<<"Save profile container successfully."<<std::endl;
		}

		void load( const std::string filename )
		{
			using namespace H5;
			H5::H5File file;
			try {
				file = H5File(filename.c_str(), H5F_ACC_RDONLY);
			}
			catch (H5::Exception& e) {
				std::string msg(std::string("could not open HDF5 file \n") + e.getCDetailMsg());
				throw StatisticalModelException(msg.c_str());
			}
			Group groupRoot = file.openGroup("/");
			try{
				this->numberOfPointsEachShape = HDF5Utils::readInt(groupRoot, "numberOfPointsEachShape");
				this->numberOfShapes = HDF5Utils::readInt(groupRoot, "numberOfShapes");
				this->profileDimension = HDF5Utils::readInt(groupRoot, "profileDimension");
				HDF5Utils::readMatrix(groupRoot, "profiles", this->profiles);
				HDF5Utils::readVector(groupRoot, "classLabels", this->classLabels);
				HDF5Utils::readVector(groupRoot, "pointIDs", this->pointIDs);
				HDF5Utils::readVector(groupRoot, "clusterLabels", this->clusterLabels);

			}catch (H5::Exception& e){
				std::string msg(std::string("an exception occurred while reading HDF5 file \n") + e.getCDetailMsg());
				std::cerr<<msg<<std::endl;
				throw StatisticalModelException(msg.c_str());
			}
			groupRoot.close();
			file.close();

			this->updateCluseterLabelSet();

			std::cout<<"Load profile container successfully."<<std::endl;
		}

		LONGTYPE getRowsOfProfiles(int cluster_label = defaultClusterLabel)
		{
			LONGTYPE cnt = 0;
			for (int ptid=0;ptid<this->clusterLabels.size();ptid++){
				if (cluster_label == static_cast<int>(clusterLabels[ptid])) {
					cnt++;
				}
			}
			return cnt*this->numberOfShapes*this->numberOfSamplesPerLandmark;
		}

		void copyCluster(const ProfileContainer & source)
		{
			if (clusterLabels.size() != source.clusterLabels.size()){
				clusterLabels.setZero(source.clusterLabels.size());
			}

			for (int ptid=0;ptid<clusterLabels.size();ptid++){
				this->clusterLabels[ptid] = source.clusterLabels[ptid];
			}
			this->updateCluseterLabelSet();
		}

		void cluster( int cluster_number_request )
		{
			std::cout<<"Start to cluster."<<std::endl;
			if (cluster_number_request<2)
			{
				std::cerr<<"Cluster number requested < 2. No need to cluster."<<std::endl;
				return;
			}
			const unsigned int measureLength = this->profileDimension*this->numberOfShapes*numberOfSamplesPerLandmark;
			const unsigned int measureCount = this->numberOfPointsEachShape;

			std::cout<<"measure length: "<<measureLength<<std::endl;
			std::cout<<"measure count:"<<measureCount<<std::endl;

			typedef itk::VariableLengthVector< double> MeasurementVectorType;
			typedef itk::Statistics::ListSample< MeasurementVectorType > SampleType;
			SampleType::Pointer sample = SampleType::New();
			sample->SetMeasurementVectorSize( measureLength );
			MeasurementVectorType mv;
			mv.SetSize( measureLength );
			mv.Fill( 0 );

			for (int ptid=0;ptid<this->numberOfPointsEachShape;ptid++)
			{
				for (int shapeid=0;shapeid<this->numberOfShapes;shapeid++)
				{
					for (int k=0;k<numberOfSamplesPerLandmark;k++)
					{
						for (int dim=0;dim<this->profileDimension;dim++)
						{
							int measureIdx = shapeid*numberOfSamplesPerLandmark*profileDimension + k*profileDimension + dim;
							int profileRowId = shapeid*numberOfPointsEachShape*numberOfSamplesPerLandmark + ptid*numberOfSamplesPerLandmark + k;
							mv[measureIdx] = this->profiles(profileRowId, dim);
						}
					}
				}
				sample->PushBack( mv );
			}
			typedef itk::Statistics::WeightedCentroidKdTreeGenerator< SampleType > TreeGeneratorType;
			TreeGeneratorType::Pointer treeGenerator = TreeGeneratorType::New();
			treeGenerator->SetSample( sample );
			treeGenerator->SetBucketSize( 16 );
			treeGenerator->Update();

			typedef TreeGeneratorType::KdTreeType TreeType;
			typedef itk::Statistics::KdTreeBasedKmeansEstimator<TreeType> EstimatorType;
			EstimatorType::Pointer estimator = EstimatorType::New();
			EstimatorType::ParametersType initialMeans( measureLength * cluster_number_request );
			initialMeans.Fill(0.0);
			estimator->SetParameters( initialMeans );
			estimator->SetKdTree( treeGenerator->GetOutput() );
			estimator->SetMaximumIteration( 500 );
			estimator->SetCentroidPositionChangesThreshold(1.0);
			estimator->StartOptimization();

			EstimatorType::ParametersType estimatedMeans = estimator->GetParameters();
			typedef itk::Statistics::DistanceToCentroidMembershipFunction< MeasurementVectorType > MembershipFunctionType;
			typedef itk::Statistics::MinimumDecisionRule DecisionRuleType;
			DecisionRuleType::Pointer decisionRule = DecisionRuleType::New();
			typedef itk::Statistics::SampleClassifierFilter< SampleType > ClassifierType;
			ClassifierType::Pointer classifier = ClassifierType::New();
			classifier->SetDecisionRule( decisionRule );
			classifier->SetInput( sample );
			classifier->SetNumberOfClasses( cluster_number_request );

			typedef ClassifierType::ClassLabelVectorObjectType
				ClassLabelVectorObjectType;
			typedef ClassifierType::ClassLabelVectorType ClassLabelVectorType;
			typedef ClassifierType::ClassLabelType ClassLabelType;
			ClassLabelVectorObjectType::Pointer classLabelsObject =
				ClassLabelVectorObjectType::New();
			ClassLabelVectorType& classLabelsVector = classLabelsObject->Get();
			for (int k=0;k<cluster_number_request;k++)
			{
				classLabelsVector.push_back( k );
			}
			classifier->SetClassLabels( classLabelsObject );

			typedef ClassifierType::MembershipFunctionVectorObjectType MembershipFunctionVectorObjectType;
			typedef ClassifierType::MembershipFunctionVectorType MembershipFunctionVectorType;
			MembershipFunctionVectorObjectType::Pointer membershipFunctionVectorObject =
				MembershipFunctionVectorObjectType::New();
			MembershipFunctionVectorType& membershipFunctionVector =
				membershipFunctionVectorObject->Get();

			int index = 0;
			for ( unsigned int i = 0 ; i < cluster_number_request ; i++ ){
				MembershipFunctionType::Pointer membershipFunction = MembershipFunctionType::New();
				MembershipFunctionType::CentroidType centroid( sample->GetMeasurementVectorSize() );
				for ( unsigned int j = 0 ; j < sample->GetMeasurementVectorSize(); j++ ){
					centroid[j] = estimatedMeans[index++];
				}
				membershipFunction->SetCentroid( centroid );
				membershipFunctionVector.push_back( membershipFunction.GetPointer() );
			}
			classifier->SetMembershipFunctions( membershipFunctionVectorObject );
			classifier->Update();

			const ClassifierType::MembershipSampleType* membershipSample = classifier->GetOutput();
			ClassifierType::MembershipSampleType::ConstIterator labelIt = membershipSample->Begin();
			ClassifierType::MembershipSampleType::ConstIterator labelItEnd = membershipSample->End();
			index = 0;
			while(labelIt!=labelItEnd)
			{
				this->clusterLabels[index++] = static_cast<int>(labelIt.GetClassLabel());
				++labelIt;
			}

			this->updateCluseterLabelSet();

			std::cout<<"Cluster done."<<std::endl;
		}

		void print()
		{
			std::cout<<"*********************ProfileClassifier***************************"<<std::endl;
			std::cout<<"Number of shape: "<<this->numberOfShapes<<std::endl;
			std::cout<<"Number of points on each shape: "<<this->numberOfPointsEachShape<<std::endl;
			std::cout<<"Number of profiles: "<<this->getRowsOfProfiles()<<std::endl;
			std::cout<<"Number of clusters: "<<this->clusterLabelSet.size()<<std::endl;
			std::cout<<"********************************************************************"<<std::endl;
		}
	};

	template<class TMesh>
	class ProfileContainerUtils
	{
	public:
		static void assignClusterLabels(ProfileContainer& profileContainer, typename TMesh * mesh)
		{
			mesh->GetPointData()->Reserve(mesh->GetNumberOfPoints());
			for (int pid=0;pid<mesh->GetNumberOfPoints();pid++)
			{
				mesh->SetPointData( pid, profileContainer.clusterLabels[pid] );
			}
		}
	};
}
#endif