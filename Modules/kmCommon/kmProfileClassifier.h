#ifndef __kmProfileClassifier_h
#define __kmProfileClassifier_h

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <math.h>
#include <set>
#include <fstream>

#include "itkMacro.h"
#include "itkCovariantVector.h"
#include "itkSimplexMeshGeometry.h"
#include "itkVector.h"
#include "itkPoint.h"
#include "itkMesh.h"

#include "statismo/CommonTypes.h"
#include "statismo/HDF5Utils.h"

#include "AdaBoost.h"
#include "kmGlobal.h"
#include "kmProfileContainer.h"

using namespace std;
using namespace itk;
using namespace statismo;

namespace km
{
	typedef unsigned long long LONGTYPE;
	
	typedef std::map<BClassType, double> ProbabilityContainer;
	typedef std::map<BClassType, double> DistanceContainer;
	/************************************************************************************************/
	// Adaboost profile classifier
	/************************************************************************************************/
	class ProfileClassifierUnit
	{
	public:
		int clusterLabel; //label for k means clustering
		PROFILE_CATEGORY profileCategory;

		std::vector<FLOATTYPE> X;
		std::vector<char>   Y;

		std::vector<int> sign;
		std::vector<FLOATTYPE> alpha;
		std::vector<FLOATTYPE> threshold;
		std::vector<int> featureID;
		std::vector<FLOATTYPE> error;

		VectorType signVec;
		VectorType alphaVec;
		VectorType thresholdVec;
		VectorType featureIDVec;
		VectorType errorVec;

		ProfileClassifierUnit( PROFILE_CATEGORY profile_category_, int clusterLabel_ )
		{
			this->clusterLabel = clusterLabel_;
			this->profileCategory = profile_category_;
		}

		~ProfileClassifierUnit()
		{
		}

		void save(const H5::CommonFG& modelRoot)
		{
			this->stdVectorToStatismoVector();

			using namespace H5;
			try {
				// create the group structure
				char groupname[128];
				sprintf(groupname, "%s%d", "./cluster_", clusterLabel);
				Group clusterGroup = modelRoot.createGroup(groupname);
				HDF5Utils::writeVector(clusterGroup, "sign", this->signVec);
				HDF5Utils::writeVector(clusterGroup, "alpha", this->alphaVec);
				HDF5Utils::writeVector(clusterGroup, "threshold", this->thresholdVec);
				HDF5Utils::writeVector(clusterGroup, "featureID", this->featureIDVec);
				HDF5Utils::writeVector(clusterGroup, "error", this->errorVec);
				clusterGroup.close();
			} catch (H5::Exception& e) {
				std::string msg(std::string("an exception occurred while writing HDF5 file \n") + e.getCDetailMsg());
				throw StatisticalModelException(msg.c_str());
			}
		}
		void load(const H5::CommonFG& modelRoot)
		{
			using namespace H5;
			try {
				// create the group structure
				char groupname[128];
				sprintf(groupname, "%s%d", "./cluster_", clusterLabel);
				Group clusterGroup = modelRoot.openGroup(groupname);
				HDF5Utils::readVector(clusterGroup, "sign", this->signVec);
				HDF5Utils::readVector(clusterGroup, "alpha", this->alphaVec);
				HDF5Utils::readVector(clusterGroup, "threshold", this->thresholdVec);
				HDF5Utils::readVector(clusterGroup, "featureID", this->featureIDVec);
				HDF5Utils::readVector(clusterGroup, "error", this->errorVec);
				clusterGroup.close();
			} catch (H5::Exception& e) {
				std::string msg(std::string("an exception occurred while writing HDF5 file \n") + e.getCDetailMsg());
				throw StatisticalModelException(msg.c_str());
			}

			this->statismoVectorToStdVector();
		}

		void pushBack(VectorType classLabels, MatrixType& profiles, LONGTYPE rowid)
		{
			Y.push_back(classLabels[rowid]);
			for (int i=0;i<profiles.cols();i++)
			{
				X.push_back(profiles(rowid, i));
			}
		}

		void train(const char* outputdir)
		{
			char adaboostFile[1024];
			if (outputdir != NULL){
				sprintf (adaboostFile, "%s/AdaBoostResult_%d",outputdir, clusterLabel);
			}else{
				sprintf (adaboostFile, "AdaBoostResult_%d", clusterLabel);
			}
			//std::cout<<adaboostFile<<std::endl;
			try{
				AdaBoostTrain(&X[0],&Y[0],Y.size(),X.size()/Y.size(), 100, adaboostFile);
				this->readTrainFile(adaboostFile, 100);
			}catch(const exception &e){
				std::cout<<e.what()<<std::endl;
			}
		}

		void readTrainFile( const char* adaboostFile, int maxWeakClassifier = 100)
		{
			ifstream ifs; 
			ifs.open( adaboostFile , ifstream::in );

			FLOATTYPE t;
			int LC;
			int weak_classifier_number = 0;
			while (ifs.good() && weak_classifier_number <= maxWeakClassifier)
			{
				ifs >> t;
				if (!ifs.good())
					break;
				ifs >> t;
				this->alpha.push_back(t);
				ifs >> t;
				ifs >> t;
				this->featureID.push_back(int(t));
				ifs >> t;
				this->sign.push_back(int(t));
				ifs >> t;
				this->threshold.push_back(t);
				ifs >> t;
				this->error.push_back(t);

				weak_classifier_number++;
			}
			ifs.close();
			LC = this->alpha.size();
			cout<<"# weak learners: "<<LC<<endl;
		}

		void stdVectorToStatismoVector()
		{
			int LC = sign.size();
			std::cout<<"Call stdVectorToStatismoVector.. LC: "<<LC<<std::endl;
			signVec.setOnes(LC);
			alphaVec.setOnes(LC);
			thresholdVec.setOnes(LC);
			featureIDVec.setOnes(LC);
			errorVec.setOnes(LC);
			for (int i=0;i<LC;i++)
			{
				signVec[i]      = sign[i];
				alphaVec[i]     = alpha[i];
				thresholdVec[i] = threshold[i];
				featureIDVec[i] = featureID[i];
				errorVec[i]     = error[i];
			}
		}

		void statismoVectorToStdVector()
		{
			int LC = signVec.size();
			//std::cout<<"Call statismoVectorToStdVector.. LC: "<<LC<<std::endl;
			for (int i=0;i<LC;i++)
			{
				sign.push_back(static_cast<int>(signVec[i]));
				alpha.push_back(static_cast<FLOATTYPE>(alphaVec[i]));
				threshold.push_back(static_cast<FLOATTYPE>(thresholdVec[i]));
				featureID.push_back(static_cast<int>(featureIDVec[i]));
				error.push_back(static_cast<FLOATTYPE>(errorVec[i]));
			}
		}

		int label2Binary( BClassType c )
		{
			if (profileCategory == BOUNDARY){
				if (c == BPClass){
					return 1;
				}else{
					return 0;
				}
			}else if(profileCategory == LIVER){
				if (c == IPClass){
					return 1;
				}else{
					return 0;
				}
			}else{
				if (c == OPClass){
					return 1;
				}else{
					return 0;
				}
			}
		}

		double test( ProfileContainer & profileContainer )
		{
			LONGTYPE NSample = 0;
			LONGTYPE NSuccess = 0;
			for (LONGTYPE row_id=0;row_id<profileContainer.profiles.rows();row_id++)
			{
				int cl = profileContainer.getClusterLabelByRowId(row_id);
				if (cl == this->clusterLabel){
					BClassType label = (BClassType)profileContainer.classLabels[row_id];
					FLOATTYPE p = this->classify(profileContainer.profiles, row_id);

					NSample++;
					if ((p>=0.5 && label2Binary(label)==1) || (p<0.5 && label2Binary(label)==0) ){
						NSuccess++;
					}
				}else{
					//Ignore this one.
				}
			}
			double successRatio = static_cast<double>(NSuccess)/NSample;
			std::cout<<"Cluster "<<this->clusterLabel<<". NSample="<<NSample<<", NSuccess="<<NSuccess<<", NSuccess/NSample="<<successRatio<<std::endl;
			return successRatio;
		}

		FLOATTYPE classify(const MatrixType & profiles, LONGTYPE rowid)
		{
			const int NFeature = profiles.cols();
			std::vector<FLOATTYPE> sample;
			for (int dim=0;dim<NFeature;dim++)
			{
				sample.push_back(profiles(rowid, dim));
			}
			return this->classify(sample);
		}

		FLOATTYPE classify( const std::vector<FLOATTYPE> & sample )
		{
			const int NFeature = sample.size();
			int LC = this->alpha.size();

			FLOATTYPE* X;
			X = const_cast<FLOATTYPE*>(&sample[0]);

			FLOATTYPE H,cH;
			AdaBoostClassify(X, 1, NFeature, LC, &featureID[0], &alpha[0], &sign[0], &threshold[0], &H, &cH);

			return H;
		}
	};

	class ProfileClassifier
	{
	public:
		typedef std::map<int, ProfileClassifierUnit*> ClassifierUnitMapType;
		typedef std::map<int, float> IntFloatMapType;
		PROFILE_CATEGORY profileCategory;
		int profileDimension;
		VectorType clusterLabels;
		std::set<int> clusterLabelSet;
		ClassifierUnitMapType classifierUnitMap;

		ProfileClassifier()
		{
			profileCategory = DEFAULT;
			this->profileDimension = 1;
		}

		~ProfileClassifier()
		{
		}

		void updateClusterLabelSet()
		{
			clusterLabelSet.clear();
			for (int i=0;i<clusterLabels.size();i++)
			{
				clusterLabelSet.insert(clusterLabels[i]);
			}
		}

		ProfileClassifierUnit* getAdaboostUnitByPointIndex( int ptidx )
		{
			int cl=0;
			if (ptidx>=this->clusterLabels.size()){
				std::cerr<<"Cannot find classifier unit for ptidx "<<ptidx<<std::endl;
				return NULL;
			}else{
				cl = this->clusterLabels[ptidx];
			}
			return this->getClassifierUnitByClusterLabel(cl);
		}

		ProfileClassifierUnit* getClassifierUnitByClusterLabel( int cluster_label )
		{
			ProfileClassifierUnit* unit = classifierUnitMap[cluster_label];
			if ( unit == NULL ){
				std::cerr<<"Cannot find classifier unit for label "<<cluster_label<<std::endl;
				return NULL;
			}else{
				return unit;
			}
		}

		double classify( const std::vector<FLOATTYPE> & sample, int ptidx = 0 )
		{
			if (sample.size() != this->profileDimension)
			{
				std::cerr<<"[Error] Sample dimension: "<<sample.size()<<". Require: "<<this->profileDimension<<std::endl;
				return 0;
			}

			ProfileClassifierUnit* adaboostUnit = this->getAdaboostUnitByPointIndex(ptidx);
			if (adaboostUnit==NULL){
				std::cerr<<"Cannot find any profile unit for point: "<<ptidx<<std::endl;
				return 0;
			}
			return adaboostUnit->classify(sample);
		}

		void train(ProfileContainer & profileContainer, PROFILE_CATEGORY category, const char* outputdir = NULL)
		{
			std::cout<<"Start to train profile classifier.."<<std::endl;

			this->profileCategory = category;
			this->profileDimension = profileContainer.profileDimension;

			std::cout<<"Copy cluster labels.."<<std::endl;
			this->clusterLabels.setZero(profileContainer.clusterLabels.size());
			for (int ptid=0;ptid<profileContainer.clusterLabels.size();ptid++)
			{
				this->clusterLabels[ptid] = profileContainer.clusterLabels[ptid];
			}

			std::cout<<"Get cluster set."<<std::endl;
			this->updateClusterLabelSet();
			std::cout<<"Number of clusters: "<<clusterLabelSet.size()<<std::endl;
			
			std::cout<<"Initialize classifier map."<<std::endl;
			for (std::set<int>::iterator clusterLabelIt = clusterLabelSet.begin(); clusterLabelIt!=clusterLabelSet.end(); clusterLabelIt++)
			{
				int cl = *clusterLabelIt;
				ProfileClassifierUnit * unit = new ProfileClassifierUnit(profileCategory, cl);
				classifierUnitMap[cl] = unit;
			}

			std::cout<<"Push all profiles to classifier units."<<std::endl;
			for (LONGTYPE rowid=0;rowid<profileContainer.profiles.rows();rowid++)
			{
				int cl = profileContainer.getClusterLabelByRowId(rowid);
				ProfileClassifierUnit * unit = this->getClassifierUnitByClusterLabel(cl);
				unit->Y.push_back(unit->label2Binary(profileContainer.classLabels[rowid]));
				for (int dim=0;dim<profileContainer.profiles.cols();dim++)
				{
					unit->X.push_back(static_cast<FLOATTYPE>(profileContainer.profiles(rowid, dim)));
				}
			}

			for (std::set<int>::iterator clusterLabelIt = clusterLabelSet.begin(); clusterLabelIt!=clusterLabelSet.end(); clusterLabelIt++)
			{
				int cl = *clusterLabelIt;
				ProfileClassifierUnit * unit = this->getClassifierUnitByClusterLabel(cl);
				unit->train(outputdir);
			}
			std::cout<<"Train profile classifier successfully.."<<std::endl;
		}

		void save(const std::string filename)
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
				//Store basic information.
				HDF5Utils::writeVector(groupRoot, "clusterLabels", this->clusterLabels);
				HDF5Utils::writeInt(groupRoot, "profileDimension", this->profileDimension);
				HDF5Utils::writeString(groupRoot, "profileCategory", ProfileCategoryUtils::category2string(this->profileCategory)); 
				
				//Store cluster information.
				for (std::set<int>::iterator clusterLabelIt = clusterLabelSet.begin(); clusterLabelIt!=clusterLabelSet.end(); clusterLabelIt++)
				{
					int cl = *clusterLabelIt;
					classifierUnitMap[cl]->save(groupRoot);
				}
			} catch (H5::Exception& e) {
				std::string msg(std::string("an exception occurred while writing HDF5 file \n") + e.getCDetailMsg());
				std::cerr<<msg<<std::endl;
				throw StatisticalModelException(msg.c_str());
			}

			groupRoot.close();
			file.close();	

			std::cout<<"Save profile classifier successfully."<<std::endl;
		}

		void load(const std::string filename)
		{
			using namespace H5;
			H5File file;
			try {
				file = H5::H5File( filename.c_str(), H5F_ACC_RDONLY);
			} catch (FileIException& e) {
				std::string msg(std::string("Could not open HDF5 file for writing \n") + e.getCDetailMsg());
				throw StatisticalModelException(msg.c_str());
			}
			H5::Group groupRoot = file.openGroup("/");
			try {
				//Load basic information.
				HDF5Utils::readVector(groupRoot, "clusterLabels", this->clusterLabels);
				this->profileDimension = HDF5Utils::readInt(groupRoot, "profileDimension");
				const std::string categorystr = HDF5Utils::readString(groupRoot, "profileCategory");
				this->profileCategory = ProfileCategoryUtils::string2category(categorystr);

				//Load cluster information.
				this->updateClusterLabelSet(); //Need to get cluster number before load cluster information.
				for (std::set<int>::iterator clusterLabelIt = clusterLabelSet.begin(); clusterLabelIt!=clusterLabelSet.end(); clusterLabelIt++)
				{
					int cl = *clusterLabelIt;
					ProfileClassifierUnit * unit = new ProfileClassifierUnit(this->profileCategory, cl);
					unit->load(groupRoot);
					classifierUnitMap[cl] = unit;
				}
			} catch (H5::Exception& e) {
				std::string msg(std::string("an exception occurred while writing HDF5 file \n") + e.getCDetailMsg());
				std::cerr<<msg<<std::endl;
				throw StatisticalModelException(msg.c_str());
			}
			groupRoot.close();
			file.close();	

			std::cout<<"Load profile classifier successfully."<<std::endl;
		}

		void test(ProfileContainer & profileContainer, IntFloatMapType & errorMap)
		{
			double ratio_weighted = 0.0;
			for (std::set<int>::iterator clusterLabelIt = clusterLabelSet.begin(); clusterLabelIt!=clusterLabelSet.end(); clusterLabelIt++)
			{
				int cl = *clusterLabelIt;
				double ratio = classifierUnitMap[cl]->test(profileContainer);
				errorMap[cl] = ratio;

				ratio_weighted += ratio*profileContainer.getRowsOfProfiles(cl);
			}
			ratio_weighted /= profileContainer.profiles.rows();
			std::cout<<"Test done!!! Weighted classified success ratio: "<<ratio_weighted<<std::endl;	
		}
		
		void test(ProfileContainer & classifier)
		{
			IntFloatMapType errorMap;
			test(classifier, errorMap);
		}

		void print()
		{
			std::cout<<"*********************ProfileClassifier***************************"<<std::endl;
			std::cout<<"Dimension of profile: "<<this->profileDimension<<std::endl;
			std::cout<<"Number of clusters: "<<this->clusterLabelSet.size()<<std::endl;
			std::cout<<"Category of profile: "<<ProfileCategoryUtils::category2string(this->profileCategory)<<std::endl;
			std::cout<<"********************************************************************"<<std::endl;
		}
	};
}

#endif