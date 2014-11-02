#ifndef __kmKNNProfileClassifier_h
#define __kmKNNProfileClassifier_h

#include "itkMacro.h"
#include "itkCovariantVector.h"

#include <fstream>						// file I/O
#include "ANN/ANN.h"

using namespace std;
using namespace itk;

#define CLASS_DIM 3 //profile类别，分为表面外，表面内，表面上以及未知
#define NEIGHBOR_DIM 5 //在进行KNN分类的时候，选出TOP N样本作为分类依据
#define SHIFT_INSIDE 5 //采样表面内profile时，样本从表面上向表面内偏移的距离
#define SHIFT_OUTSIDE 5 //采样表面外profile时，样本从表面上向表面外偏移的距离
#define SIGMA 1.0 //采样梯度profile时，计算梯度场采样的参数sigma
#define PROFILE_DIM 11
#define PROFILE_SPACING 1.5

namespace km
{
	enum BClassType
	{
		IPClass = 0,
		BPClass = 1,
		OPClass = 2,
		//UNKNOWN = 3
	};

	template<unsigned int Dimension>
	class Profile
	{
	public:
		typedef Profile<Dimension> Self;
		typedef itk::CovariantVector<float, Dimension> ValueVectorType;

		itkStaticConstMacro(Dimension, unsigned, Dimension);

		double spacing;
		unsigned int index;
		unsigned int kmclass;
		BClassType bclass;
		ValueVectorType value;

		Profile():spacing(1.0), index(0), bclass(OPClass), kmclass(0)
		{

		}

		~Profile()
		{

		}

		double GetNorm()
		{
			return value.GetNorm();
		}

		void Normalize()
		{
			double vv = this->GetNorm();
			if (vv==0)
			{
				this->value.Fill(0);
			}
			else
			{
				value.Normalize();
			}
		}

		void Fill( double c )
		{
			value.Fill(c);
		}

		void Init()
		{
			spacing = 1.0;
			index = 0;
			bclass = OPClass;
			kmclass = 0;
			value.Fill(0.0);
		}

		friend ostream& operator <<(ostream &os,const Self &self)
		{
			os << self.index << '\t' << self.bclass << '\t' << self.kmclass << '\t';
			for (int i=0;i<Dimension-1;i++)
			{
				os << self.value[i] << '\t';
			}
			os << self.value[Dimension-1];

			return os;
		}

		float operator[](short index)
		{
			return value[index]; 
		}

		Self & operator=(const Self & self)
		{
			this->index = self.index;
			this->bclass = self.bclass;
			this->spacing = self.spacing;
			this->kmclass = self.kmclass;

			this->value = self.value;

			return *this;
		}
	};

	template<class ProfileType>
	class
		ProfilePackage
	{
	public:
		typedef ProfilePackage<ProfileType> Self;
		typedef typename ProfileType ProfileType;
		typedef std::pair<unsigned int, double> KmeansPairType;
		typedef std::map<unsigned int, double> KmeansPairContainerType;
		typedef typename KmeansPairContainerType::iterator KmeansPairIterator;

		typedef std::vector<ProfileType> ProfileContainerType;

		KmeansPairContainerType kmeansPairs;
		ProfileContainerType profiles;

		unsigned int numberOfKmeansPairs;
		unsigned int numberOfProfiles;

		bool weighted;

		ProfilePackage():numberOfProfiles(0), weighted(false)
		{
			kmeansPairs.insert( KmeansPairType(0, 1.0));
			numberOfKmeansPairs = 1;
		}

		~ProfilePackage()
		{
			kmeansPairs.clear();
			profiles.clear();
		}

		void Clear()
		{
			profiles.clear();
			numberOfProfiles = 0;

			kmeansPairs.clear();
			kmeansPairs.insert(0, 1.0);
			numberOfKmeansPairs = 1;
		}

		void ClearProfiles()
		{
			profiles.clear();
			numberOfProfiles = 0;
		}

		void Init()
		{
			profiles.clear();
			numberOfProfiles = 0;

			kmeansPairs.clear();
			kmeansPairs.insert( KmeansPairType(0, 1.0) );
			numberOfKmeansPairs = 1;
		}

		bool AddKmeansPair( unsigned int kmclass, double kmweight )
		{

			if (kmeansPairs.find(kmclass) != kmeansPairs.end())
			{
				//Already exists. So just modify the old value.
				kmeansPairs[kmclass] = kmweight;
				return true;
			}

			kmeansPairs.insert( KmeansPairType( kmclass, kmweight ) );

			numberOfKmeansPairs++;

			return true;
		}

		bool AddProfile(const typename ProfileType& profile)
		{
			profiles.push_back(profile);

			numberOfProfiles++;

			return true;
		}

		void SetWeightedTrue()
		{
			weighted = true;
		}

		void SetWeightedFalse()
		{
			weighted = false;
		}

		friend ostream& operator <<(ostream &os,const Self &self)
		{
			os << "Number of kmeans: " << self.numberOfKmeansPairs << std::endl;
			os << "Number of profiles: " << self.numberOfProfiles << std::endl;
			os << "Weighted: " << self.weighted << std::endl;
			return os;
		}
	};

	template<class ProfileType>
	void 
		profile2point( const ProfileType & pro, ANNpoint & pt )
	{
		for (int i=0;i<ProfileType::Dimension;i++)
		{
			pt[i] = pro.value[i];
		}
	}

	template<class ProfileType>
	bool 
		readProfilesFromFile( const char* filename, std::vector<ProfileType>& profiles )
	{
		ifstream in;
		in.open( filename , ios::in );// open data file
		if (!in) {
			cerr << "Cannot open data file: "<<filename<<"\n";
			return false;
		}

		unsigned int count = 0;

		while(true)
		{
			ProfileType p;
			int index = 0;
			if (in >> index)
			{
				p.index = index;
			}
			else
			{
				break;
			}

			int intBClass = 0;
			if (in >> intBClass)
			{
				p.bclass = (BClassType)intBClass;
			}
			else
			{
				break;
			}

			unsigned int kmclass = 0;
			if (in >> kmclass)
			{
				p.kmclass = kmclass;
			}
			else
			{
				break;
			}

			for (int i = 0; i < ProfileType::Dimension; i++)
			{
				if(!(in >> p.value[i])) 
					break;
			}
			profiles.push_back( p );
			count++;
		}

		in.close();
		std::cout<<"[readProfilesFromFile] Load points: "<<count<<std::endl;

		return true;
	}

	template<class ProfileType>
	void writeProfilesToFile(const char* filename, std::vector<ProfileType>& allprofiles)
	{
		std::ofstream os;
		os.open( filename );

		os << "PROFILES" << std::endl;
		for (int i=0;i<allprofiles.size();i++)
		{
			os << allprofiles[i] << std::endl;
		}

		os.close();
	}

	template<class ProfilePackage>
	bool 
		readProfilesFromFile( const char* filename, typename ProfilePackage& profiles )
	{
#define READ_WITH_ERROR_INFO(X) if(!(in>>X)){std::cerr<<"Error occurs when read from file!"<<std::endl; return false;}

		profiles.Init();

		ifstream in;
		in.open( filename , ios::in );// open data file
		if (!in) {
			cerr << "Cannot open data file: "<<filename<<"\n";
			return false;
		}

		unsigned int numberOfKmeansPairs = 0;
		unsigned int numberOfProfiles = 0;

		const unsigned int BUFFER_SIZE = 256;
		char tmpbuffer[BUFFER_SIZE];

		/************************************************************************/
		/* Read K means information                                             */
		/************************************************************************/
		in.getline (tmpbuffer,BUFFER_SIZE);
		if (strcmp(tmpbuffer, "KMEANS") != 0)
		{
			std::cerr << "Cannot find dataset KMEANS! Instead it's " << tmpbuffer << std::endl;
			return false;
		}
		else
		{
			READ_WITH_ERROR_INFO(tmpbuffer);
			if (strcmp(tmpbuffer, "NUMBER") != 0)
			{
				std::cerr << "Cannot find NUMBER of KMEANS! Instead it's " << tmpbuffer << std::endl;
				return false;
			}
			else
			{
				READ_WITH_ERROR_INFO(numberOfKmeansPairs);
				//std::cout << "Number of kmeans pairs: "<< numberOfKmeansPairs << std::endl;

				for (int k=0;k<numberOfKmeansPairs;k++)
				{
					unsigned int kmeanclass = 0;
					double kmeanweight = 0.0;

					READ_WITH_ERROR_INFO(kmeanclass);
					READ_WITH_ERROR_INFO(kmeanweight);

					profiles.AddKmeansPair( kmeanclass, kmeanweight );
				}
			}
		}

		/************************************************************************/
		/* Read Profile list.                                                   */
		/************************************************************************/
		in.getline (tmpbuffer,1); //Need to change line here!!!!
		in.getline (tmpbuffer,BUFFER_SIZE);
		if (strcmp(tmpbuffer, "PROFILES") != 0)
		{
			std::cerr << "Cannot find dataset PROFILES! Instead it's " << tmpbuffer << std::endl;
			return false;
		}
		else
		{
			READ_WITH_ERROR_INFO(tmpbuffer);
			if (strcmp(tmpbuffer, "NUMBER") != 0)
			{
				std::cerr << "Cannot find NUMBER of PROFILES! Instead it's " << tmpbuffer << std::endl;
				return false;
			}
			else
			{
				READ_WITH_ERROR_INFO(numberOfProfiles);
				//std::cout << "Number of kmeans pairs: "<< numberOfProfiles << std::endl;

				for (int k=0;k<numberOfProfiles;k++)
				{
					typedef ProfilePackage::ProfileType ProfileType;
					ProfileType p;
					READ_WITH_ERROR_INFO(p.index);

					int intBClass = 0;
					READ_WITH_ERROR_INFO(intBClass);
					p.bclass = (BClassType)intBClass;

					READ_WITH_ERROR_INFO(p.kmclass);

					for (int i = 0; i < ProfileType::Dimension; i++)
					{
						READ_WITH_ERROR_INFO(p.value[i]);
					}

					profiles.AddProfile( p );
				}
			}
		}

		return true;

#undef READ_WITH_ERROR_INFO(X)
	}

	template<class ProfilePackage>
	void writeProfilesToFile(const char* filename, typename ProfilePackage& allprofiles)
	{
		std::ofstream os;
		os.open( filename );

		os << "KMEANS"<< std::endl;
		os << "NUMBER "<< allprofiles.numberOfKmeansPairs << std::endl;

		ProfilePackage::KmeansPairIterator kmPairIt = allprofiles.kmeansPairs.begin();
		ProfilePackage::KmeansPairIterator kmPairItEnd = allprofiles.kmeansPairs.end();
		while( kmPairIt!=kmPairItEnd )
		{
			os << kmPairIt->first << " " << kmPairIt->second << std::endl;
			kmPairIt++;
		}

		os << "PROFILES" << std::endl;
		os << "NUMBER "<<allprofiles.numberOfProfiles<<std::endl;
		for (int i=0;i<allprofiles.numberOfProfiles;i++)
		{
			os << allprofiles.profiles[i] << std::endl;
		}

		os.close();
	}

		/************************************************************************
	
	提取梯度&灰度profile
	
	*************************************************************************/
	template<class GradientInterpolatorType, class IntensityInterpolatorType, class PointType, class NormalType, class ProfileType>
	void
		extractAdvancedProfile(
		const typename GradientInterpolatorType* gradientInterpolator, 
		const typename IntensityInterpolatorType* intensityInterpolator,
		const typename PointType ipoint, 
		const typename NormalType normal, 
		ProfileType & profile)
	{
		km::extractGradientProfile(gradientInterpolator, ipoint, normal, profile);
	}

	/************************************************************************
	
	提取梯度&灰度profile
	
	*************************************************************************/
	template<class GradientInterpolatorType, class IntensityInterpolatorType, class PointType, class NormalType, class ProfileType>
	void
		extractGradientAndIntensityProfile(
		const typename GradientInterpolatorType* gradientInterpolator, 
		const typename IntensityInterpolatorType* intensityInterpolator,
		const typename PointType ipoint, 
		const typename NormalType normal, 
		ProfileType & profile)
	{
		typedef typename GradientInterpolatorType::InputImageType GradientImageType;
		typedef typename IntensityInterpolatorType::InputImageType IntensityImageType;
		typedef typename GradientImageType::PixelType GradientPixelType;
		typedef typename IntensityImageType::PixelType IntensityPixelType;

		PointType profilePos;
		profile.Fill(0);

		const unsigned PointDimension   = PointType::PointDimension;
		const unsigned ProfileDimension = ProfileType::Dimension;
		const double   ProfileSpacing   = profile.spacing;

		const unsigned InnerProfileDimension = ProfileDimension/2;
		typedef Profile<InnerProfileDimension> InnerProfileType;

		InnerProfileType innerProfile1;
		InnerProfileType innerProfile2;

		//First move the position to the end(inside) of profile. 
		for (int d=0;d<PointDimension;d++)
		{
			profilePos[d] = ipoint[d]-(InnerProfileDimension*0.5)*normal[d]*ProfileSpacing;
		}

		//Then start to search towards outside.
		for ( int i=0;i<InnerProfileDimension;i++ )
		{
			if (gradientInterpolator->IsInsideBuffer(profilePos))
			{
				GradientPixelType gradient = gradientInterpolator->Evaluate( profilePos );
				double intensity = intensityInterpolator->Evaluate( profilePos );
				innerProfile1.value[i] = gradient*normal;
				innerProfile2.value[i] = intensity;
			}
			else  //If current position is outside of image region, then just ignore this pixel.
			{
				innerProfile1.value[i] = 0.0;
				innerProfile2.value[i] = 0.0;
			}

			for (int d=0;d<PointDimension;d++)
			{
				profilePos[d] += normal[d]*ProfileSpacing; //One step along normal direction.
			}
		}
		if (innerProfile1.value.GetNorm()>0)
		{
			innerProfile1.value.Normalize();
		}

		if (innerProfile2.value.GetNorm()>0)
		{
			innerProfile2.value.Normalize();
		}

		for (int i=0;i<InnerProfileDimension;i++)
		{
			profile.value[i]                       = innerProfile1[i];
			profile.value[i+InnerProfileDimension] = innerProfile2[i];
		}
	}

	/************************************************************************

	提取梯度profile

	*************************************************************************/
	template<class GradientInterpolatorType, class PointType, class NormalType, class ProfileType>
	void
		extractGradientProfile(
		const typename GradientInterpolatorType* gradientInterpolator, 
		const typename PointType ipoint, 
		const typename NormalType normal, 
		ProfileType & profile)
	{
		typedef typename GradientInterpolatorType::InputImageType GradientImageType;
		typedef typename GradientImageType::PixelType GradientPixelType;

		PointType profilePos;
		profile.Fill(0);

		const unsigned PointDimension   = PointType::PointDimension;
		const unsigned ProfileDimension = ProfileType::Dimension;
		const double   ProfileSpacing   = profile.spacing;

		//First move the position to the end(inside) of profile. 
		for (int d=0;d<PointDimension;d++)
		{
			profilePos[d] = ipoint[d]-(ProfileDimension*0.7)*normal[d]*ProfileSpacing;
		}

		//Then start to search towards outside.
		for ( int i=0;i<ProfileDimension;i++ )
		{
			if (gradientInterpolator->IsInsideBuffer(profilePos))
			{
				GradientPixelType gradient = gradientInterpolator->Evaluate( profilePos );
				profile.value[i] = gradient.GetNorm();
			}
			else  //If current position is outside of image region, then just ignore this pixel.
			{
				profile.value[i] = -100.0;
			}

			for (int d=0;d<PointDimension;d++)
			{
				profilePos[d] += normal[d]*ProfileSpacing; //One step along normal direction.
			}
		}
		if (profile.value.GetNorm()>0)
		{
			profile.value.Normalize();
		}
	}

	template<class ProfileType>
	class KNNProfileClassifier
	{
	public:
		KNNProfileClassifier();
		~KNNProfileClassifier();

		typedef ProfileType ProfileType;
		typedef BClassType BClassType;
		typedef ProfilePackage<ProfileType> ProfilePackageType;

		itkStaticConstMacro( ProfileDimension, unsigned int, ProfileType::Dimension );

		bool   buildFromPoints( const std::vector<ProfileType> & points);
		double testFromPoints( const std::vector<ProfileType> & testPoints); //Return total success rate.
		void   classify( const ProfileType & queryProfile, BClassType & maximumPossibleType, double & maximumProbability );
		double queryProbablity( const ProfileType & queryProfile, BClassType searchType );
		void   queryProbablity( const ProfileType & queryProfile, double & probabilityOfInside, double & probabilityOfBoundary, double & probabilityOfOutside );

		double              eps;
		int                 maxPts;
		ANNpoint            queryPt;
		int				    nPts;					  // actual number of data points
		ANNpointArray		dataPts;				// data points
		ANNidxArray			nnIdx;					// near neighbor indices
		ANNdistArray		dists;					// near neighbor distances
		ANNkd_tree*			kdTree;					// search structure

		BClassType*         bclasses;

	private:
		std::string  m_SampleList;
		std::string  m_TestSampleList;

	};

	//////////////////////////////////////////////////////////////////////////
	//
	// Implementation of API
	//
	//////////////////////////////////////////////////////////////////////////
	template<class ProfileType>
	KNNProfileClassifier<ProfileType>::KNNProfileClassifier()
		:eps(0),maxPts(500000),m_SampleList(""),m_TestSampleList(""),bclasses(NULL),kdTree(NULL)
	{
	}

	template<class ProfileType>
	KNNProfileClassifier<ProfileType>::~KNNProfileClassifier()
	{
		delete [] nnIdx;							// clean things up
		delete [] dists;
		delete [] bclasses;
		delete kdTree;
		annClose();									// done with ANN
	}

	template<class ProfileType>
	bool KNNProfileClassifier<ProfileType>::buildFromPoints( const std::vector<ProfileType> & points)
	{
		nPts = points.size();

		if (nPts==0)
		{
			return false;
		}

		dataPts = annAllocPts(nPts, ProfileDimension);	// allocate data points
		nnIdx   = new ANNidx[NEIGHBOR_DIM];						// allocate near neigh indices
		dists   = new ANNdist[NEIGHBOR_DIM];						// allocate near neighbor dists

		bclasses = new BClassType[nPts];
		//cout << "Load Points: "<< nPts << "\n";

		for (int i=0;i<nPts;i++)
		{
			ProfileType profile = points[i];
			bclasses[i] = profile.bclass;

			profile2point(profile, dataPts[i]);
		}

		kdTree = new ANNkd_tree(					// build search structure
			dataPts,					// the data points
			nPts,						  // number of points
			ProfileDimension);						  // dimension of space

		return true;
	}

	//Return total successful rate.
	template<class ProfileType>
	double KNNProfileClassifier<ProfileType>::testFromPoints( const std::vector<ProfileType> & testProfiles)
	{
		unsigned numberOfTestPts = testProfiles.size();

		cout << "Load Testing Points: "<< numberOfTestPts << "\n";

		if (numberOfTestPts==0)
		{
			return 0.0;
		}

		unsigned int positiveCaseInside = 0;
		unsigned int positiveCaseOutside = 0;
		unsigned int positiveCaseBoundary = 0;
		unsigned int negativeCaseInside = 0;
		unsigned int negativeCaseOutside = 0;
		unsigned int negativeCaseBoundary = 0;

		for (int i=0;i<numberOfTestPts;i++)
		{
			ProfileType queryProfile = testProfiles[i];

			BClassType maximumPossibleType;
			double maximumProbablity = 0.0;

			classify(queryProfile, maximumPossibleType, maximumProbablity);

			BClassType trueType = queryProfile.bclass;
			if (trueType==BPClass)
			{
				if (maximumPossibleType == trueType){
					positiveCaseBoundary++;
				}else{
					negativeCaseBoundary++;
				}
			}
			else if (trueType==IPClass)
			{
				if (maximumPossibleType == trueType){
					positiveCaseInside++;
				}else{
					negativeCaseInside++;
				}
			}		
			else if (trueType==OPClass)
			{
				if (maximumPossibleType == trueType){
					positiveCaseOutside++;
				}else{
					negativeCaseOutside++;
				}
			}
		}

		std::cout<<"Total cases: "<<numberOfTestPts<<std::endl;
		std::cout<<"Positive test cases[class inside]: "<<positiveCaseInside<<std::endl;
		std::cout<<"Negative test cases[class inside]: "<<negativeCaseInside<<std::endl;
		std::cout<<"Success rate: "<<static_cast<double>(positiveCaseInside)/(positiveCaseInside+negativeCaseInside)<<std::endl;

		std::cout<<"Positive test cases[class outside]: "<<positiveCaseOutside<<std::endl;
		std::cout<<"Negative test cases[class outside]: "<<negativeCaseOutside<<std::endl;
		std::cout<<"Success rate: "<<static_cast<double>(positiveCaseOutside)/(positiveCaseOutside+negativeCaseOutside)<<std::endl;

		std::cout<<"Positive test cases[class boundary]: "<<positiveCaseBoundary<<std::endl;
		std::cout<<"Negative test cases[class boundary]: "<<negativeCaseBoundary<<std::endl;
		std::cout<<"Success rate: "<<static_cast<double>(positiveCaseBoundary)/(positiveCaseBoundary+negativeCaseBoundary)<<std::endl;

		double totalSuccessRate = static_cast<double>( positiveCaseBoundary+positiveCaseInside+positiveCaseOutside )/numberOfTestPts;

		std::cout<<"Total success rate: "<< totalSuccessRate <<std::endl;

		return totalSuccessRate;
	}

	template<class ProfileType>
	double KNNProfileClassifier<ProfileType>::queryProbablity( const ProfileType & queryProfile, BClassType searchType)
	{
		queryPt = annAllocPt(ProfileDimension);
		profile2point( queryProfile, queryPt );

		kdTree->annkSearch(	// search
			queryPt,					// query point
			NEIGHBOR_DIM,				// number of near neighbors
			nnIdx,					  // nearest neighbors (returned)
			dists,					  // distance (returned)
			eps);							// error bound

		unsigned int numberOfInsideP = 0;
		unsigned int numberOfOutsideP = 0;
		unsigned int numberOfBoundaryP = 0;
		for (int k=0;k<NEIGHBOR_DIM;k++)
		{
			if ( bclasses[nnIdx[k]] == BPClass ){
				numberOfBoundaryP++;
			}else if ( bclasses[nnIdx[k]] == IPClass ){
				numberOfInsideP++;
			}else if ( bclasses[nnIdx[k]] == OPClass ){
				numberOfOutsideP++;
			}
		}

		unsigned int numberOfSearchType;
		if (searchType == BPClass){
			numberOfSearchType = numberOfBoundaryP;
		}else if(searchType == IPClass){
			numberOfSearchType = numberOfInsideP;
		}else if (searchType == BPClass){
			numberOfSearchType = numberOfBoundaryP;
		}

		delete queryPt;

		return static_cast<double>(numberOfSearchType)/NEIGHBOR_DIM;
	}

	template<class ProfileType>
	void KNNProfileClassifier<ProfileType>::queryProbablity( const ProfileType & queryProfile, double & probabilityOfInside, double & probabilityOfBoundary, double & probabilityOfOutside)
	{
		queryPt = annAllocPt(ProfileDimension);
		profile2point( queryProfile, queryPt );

		kdTree->annkSearch(	// search
			queryPt,					// query point
			NEIGHBOR_DIM,				// number of near neighbors
			nnIdx,					  // nearest neighbors (returned)
			dists,					  // distance (returned)
			eps);							// error bound

		unsigned int numberOfInsideP = 0;
		unsigned int numberOfOutsideP = 0;
		unsigned int numberOfBoundaryP = 0;
		for (int k=0;k<NEIGHBOR_DIM;k++)
		{
			if ( bclasses[nnIdx[k]] == BPClass ){
				numberOfBoundaryP++;
			}else if ( bclasses[nnIdx[k]] == IPClass ){
				numberOfInsideP++;
			}else if ( bclasses[nnIdx[k]] == OPClass ){
				numberOfOutsideP++;
			}
		}

		probabilityOfInside = static_cast<double>(numberOfInsideP)/static_cast<double>( NEIGHBOR_DIM );
		probabilityOfOutside = static_cast<double>(numberOfOutsideP)/static_cast<double>( NEIGHBOR_DIM );
		probabilityOfBoundary = static_cast<double>(numberOfBoundaryP)/static_cast<double>( NEIGHBOR_DIM );

		delete queryPt;
	}

	template<class ProfileType>
	void KNNProfileClassifier<ProfileType>::classify( const ProfileType & queryProfile, BClassType & maximumPossibleType, double & maximumProbability )
	{
		queryPt = annAllocPt(ProfileDimension);
		profile2point( queryProfile, queryPt );

		kdTree->annkSearch(	// search
			queryPt,					// query point
			NEIGHBOR_DIM,				// number of near neighbors
			nnIdx,					  // nearest neighbors (returned)
			dists,					  // distance (returned)
			eps);							// error bound

		unsigned int numberOfInsideP = 0;
		unsigned int numberOfOutsideP = 0;
		unsigned int numberOfBoundaryP = 0;
		for (int k=0;k<NEIGHBOR_DIM;k++)
		{
			if ( bclasses[nnIdx[k]] == BPClass ){
				numberOfBoundaryP++;
			}else if ( bclasses[nnIdx[k]] == IPClass ){
				numberOfInsideP++;
			}else if ( bclasses[nnIdx[k]] == OPClass ){
				numberOfOutsideP++;
			}
		}

		double insideProbability = static_cast<double>(numberOfInsideP)/static_cast<double>( NEIGHBOR_DIM );
		double outsideProbability = static_cast<double>(numberOfOutsideP)/static_cast<double>( NEIGHBOR_DIM );
		double boundaryProbability = static_cast<double>(numberOfBoundaryP)/static_cast<double>( NEIGHBOR_DIM );

		if (boundaryProbability>=outsideProbability && boundaryProbability>=insideProbability){
			maximumPossibleType = BPClass;
			maximumProbability = boundaryProbability;
		}else if(insideProbability>=outsideProbability && insideProbability>=boundaryProbability){
			maximumPossibleType = IPClass;
			maximumProbability = insideProbability;
		}else{
			maximumPossibleType = OPClass;
			maximumProbability = outsideProbability;
		}

		delete queryPt;
	}
}

#endif