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
#include <itkVariableLengthVector.h>
#include "itkListSample.h"
#include "itkKdTree.h"
#include "itkWeightedCentroidKdTreeGenerator.h"
#include "itkMinimumDecisionRule.h"
#include "itkImageToListSampleAdaptor.h"
#include "itkSampleClassifierFilter.h"
#include "itkKdTreeBasedKmeansEstimator.h"
#include <itkImageMomentsCalculator.h>
#include "itkVariableLengthVector.h"

#include "flann/flann.hpp"
#include "flann/io/hdf5.h"

#include "AdaBoost.h"
#include "kmGlobal.h"

using namespace std;
using namespace itk;
using namespace flann;

#define DELETE_NULLABLE_PRT(X) if(X.ptr()!=NULL){delete[] X.ptr();}
#define MAX_NEIGHBOUR 100

namespace km
{
	typedef unsigned long long LONGTYPE;
	typedef flann::Matrix<float>         DatasetMatrixType;
	typedef flann::Matrix<BClassType>    LabelMatrixType;
	typedef flann::Matrix<float>         PointMatrixType;
	typedef flann::Matrix<int>           IndexMatrixType;
	typedef flann::Matrix<float>         DistanceMatrixType;
	typedef flann::Matrix<int>           IntegerMatrixType;
	typedef flann::Matrix<float>         FloatMatrixType;
	
	typedef std::map<BClassType, double> ProbabilityContainer;
	typedef std::map<BClassType, double> DistanceContainer;

	template<class MatrixType>
	void fillFlannMatrix( typename MatrixType& matrix, typename MatrixType::type val );
	
	double mappingItensity( double val );
	
	template<class T>
	void 
		calculateMeanAndStd(std::vector<T> & list, typename T& mean, typename T& sd);
	
	template<class T>
	void 
		normalizeList(std::vector<T> & list, typename T& mean, typename T& sd);
		
	template<class T>
	void
		normalizeList(std::vector<T> & list);
		
	template<class GradientInterpolatorType, class IntensityInterpolatorType>
	void 
		extractFeature(
			const typename GradientInterpolatorType* gradientInterpolator, 
			const typename IntensityInterpolatorType* intensityInterpolator,
			const itk::SimplexMeshGeometry* geoData,
			const itk::SimplexMeshGeometry::PointType & ipoint,
			std::vector<double> & feature,
			PROFILE_CATEGORY profile_category);

	/************************************************************************************************/
	// Adaboost profile classifier
	/************************************************************************************************/	
	class KNNProfileClassifier
	{
	public:
		class KNNUnit
		{
		public:
			typedef flann::L2<float> DistanceFunctorType;
			typedef flann::Index<DistanceFunctorType> KNNIndexType;
		
			DatasetMatrixType   dataset;
			LabelMatrixType     labels;
			flann::IndexParams  index_params;
			DistanceFunctorType distance_function;
			flann::SearchParams search_params;
			KNNIndexType      * index;
			IndexMatrixType     indices;
			DistanceMatrixType  dists;
			IndexMatrixType     indices_one;
			DistanceMatrixType  dists_one;
			int                 cluster_label; //label for k means clustering.
		private:
			unsigned long profile_count;

		public:
			KNNUnit( DatasetMatrixType& dataset_, LabelMatrixType& labels_, int cluster_label_ = 0 )
			{
				this->dataset = dataset_;
				this->labels = labels_;
				this->cluster_label = cluster_label_;
				this->profile_count = 0;

				index_params = flann::KDTreeIndexParams(4);
				distance_function = DistanceFunctorType();

				indices_one = IndexMatrixType(new int[MAX_NEIGHBOUR], 1, MAX_NEIGHBOUR);
				dists_one   = DistanceMatrixType(new float[MAX_NEIGHBOUR], 1, MAX_NEIGHBOUR);

				search_params = flann::SearchParams(128);
				index = new KNNIndexType( index_params, distance_function );
			}
			~KNNUnit()
			{
				std::cout<<"De-constructor of KNNUnit.."<<std::endl;
				DELETE_NULLABLE_PRT(dataset);
				DELETE_NULLABLE_PRT(labels);
				DELETE_NULLABLE_PRT(indices);
				DELETE_NULLABLE_PRT(dists);
				DELETE_NULLABLE_PRT(indices_one);
				DELETE_NULLABLE_PRT(dists_one);

				delete index;
				std::cout<<"De-constructor end.."<<std::endl;
			}

			DatasetMatrixType& getDataset()
			{
				return this->dataset;
			}

			LabelMatrixType& getLabels()
			{
				return this->labels;
			}

			void buildIndex()
			{
				this->index->buildIndex( this->dataset );
			}

			void push_back( int label, PointMatrixType & point )
			{
				if (profile_count >= this->dataset.rows)
				{
					std::cerr<<"Dataset is full. This should never happen! Profile count: " << profile_count <<", maximum dataset count: "<< this->dataset.rows <<std::endl;
					return;
				}

				this->labels[profile_count][0] = label;

				for (int c=0;c<point.cols;c++)
				{
					this->dataset[profile_count][c] = point[0][c];
				}

				profile_count++;
			}

			int getPointsCount()
			{
				return dataset.rows;
			}

			int getProfileDimension()
			{
				return dataset.cols;
			}

			void generateTestDataset( DatasetMatrixType& dataset_test, LabelMatrixType& labels_test, double ratio_test)
			{
				int numberOfTesting  = this->getPointsCount()*ratio_test;
				int offsetOfTesting = this->getPointsCount() - numberOfTesting;

				int profile_dimension = this->getProfileDimension();

				//std::cout<<"dimension of profile: "<<profile_dimension<<std::endl;
				//std::cout<<"number of testing points: "<<numberOfTesting<<std::endl;

				dataset_test = DatasetMatrixType(new float[numberOfTesting*profile_dimension], numberOfTesting, profile_dimension);
				labels_test = LabelMatrixType(new int[numberOfTesting], numberOfTesting, 1);

				for (int i=0;i<numberOfTesting;i++)
				{
					for (int j=0;j<profile_dimension;j++)
					{
						dataset_test[i][j] = dataset[i+offsetOfTesting][j];
					}
					labels_test[i][0]  = labels[i+offsetOfTesting][0];
				}
			}

			void classify( PointMatrixType& point, ProbabilityContainer & probablities, DistanceContainer & dists, int max_neighbour)
			{
				probablities[IPClass] = 0.0;
				probablities[BPClass] = 0.0;
				probablities[OPClass] = 0.0;

				dists[IPClass] = 0.0;
				dists[BPClass] = 0.0;
				dists[OPClass] = 0.0;

				index->knnSearch(point, indices_one, dists_one, max_neighbour, search_params);

				for (int i=0;i<max_neighbour;i++)
				{
					int label = labels[indices_one[0][i]][0];
					probablities[label] += 1;
					dists[label] += dists_one[0][i];
				}

				dists[IPClass] /= (probablities[IPClass]==0?1:probablities[IPClass]);
				dists[BPClass] /= (probablities[BPClass]==0?1:probablities[BPClass]);
				dists[OPClass] /= (probablities[OPClass]==0?1:probablities[OPClass]);

				probablities[IPClass] /= max_neighbour;
				probablities[BPClass] /= max_neighbour;
				probablities[OPClass] /= max_neighbour;
			}

			double test( double ratio_test, int max_neighbour )
			{
				DatasetMatrixType dataset_test;
				LabelMatrixType   labels_test;

				this->generateTestDataset(dataset_test, labels_test, ratio_test);

				//flann::save_to_file(dataset_test, "testing.h5", "dataset_test");
				//flann::save_to_file(labels_test, "testing.h5", "labels_test");

				indices = IndexMatrixType(new int[dataset_test.rows*max_neighbour], dataset_test.rows, max_neighbour);
				dists   = DistanceMatrixType(new float[labels_test.rows*max_neighbour], labels_test.rows, max_neighbour);

				index->knnSearch(dataset_test, indices, dists, max_neighbour, search_params);

				std::map<BClassType, double> probabilities;

				unsigned int total = 0;
				unsigned int right = 0;
				for (int i=0;i<indices.rows;i++)
				{
					probabilities[IPClass] = 0.0;
					probabilities[BPClass] = 0.0;
					probabilities[OPClass] = 0.0;

					for (int j=1;j<indices.cols;j++)
					{
						BClassType label = labels[indices[i][j]][0];
						probabilities[label]++;
					}

					BClassType classifiedLabel;
					if (probabilities[OPClass]>=probabilities[BPClass] && probabilities[OPClass]>=probabilities[IPClass])
					{
						classifiedLabel = OPClass;
					}
					else if (probabilities[BPClass]>=probabilities[OPClass] && probabilities[BPClass]>=probabilities[IPClass])
					{
						classifiedLabel = BPClass;
					}
					else
					{
						classifiedLabel = IPClass;
					}

					if (classifiedLabel == labels_test[i][0])
					{
						right++;
					}

					total++;
				}

				DELETE_NULLABLE_PRT(dataset_test);
				DELETE_NULLABLE_PRT(labels_test);
				DELETE_NULLABLE_PRT(indices);
				DELETE_NULLABLE_PRT(dists);

				return (double)right/total;
			}
		};
	public:
		int  max_neighbour;
		int  shape_number;
		int  shape_points_number;
		int  profile_dimension;

		IndexMatrixType   cluster_labels;
		IntegerMatrixType shape_number_matrix;
		IntegerMatrixType shape_points_number_matrix;
		IntegerMatrixType cluster_labels_set; //Unique cluster label. Like 0, 1, 2..

		typedef std::pair<int, KNNUnit*> KNNUnitPairType;
		typedef std::map<int, KNNUnit*> KNNUnitMapType;
		
		KNNUnitMapType knnUnitMap;

	private:
		bool is_initilized;

	public:
		KNNProfileClassifier( )
		{
			is_initilized = false;
			max_neighbour = 11;
			shape_number = 0;
			shape_points_number = 0;
			profile_dimension = 0;
		}

		~KNNProfileClassifier()
		{
			std::cout<<"De-constructor of KNNProfileClassifier.."<<std::endl;
			//DELETE_NULLABLE_PRT(cluster_labels);
			//DELETE_NULLABLE_PRT(shape_number_matrix);
			//DELETE_NULLABLE_PRT(shape_points_number_matrix);
			//DELETE_NULLABLE_PRT(cluster_labels_set);
			//delete[] cluster_labels.ptr();
			//delete[] shape_number_matrix.ptr();
			//delete[] shape_points_number_matrix.ptr();
			//delete[] cluster_labels_set.ptr();
			std::cout<<"De-constructor end.."<<std::endl;
		}

		void setMaxNeighbour(int n)
		{
			if (n>MAX_NEIGHBOUR)
			{
				std::cout<<"Number of neighbour cannot exceed "<<MAX_NEIGHBOUR<<std::endl;
			}
			this->max_neighbour = n;
		}


		void setShapeNumber(int n)
		{
			this->shape_number = n;	
		}

		void setShapePointsNumber(int n)
		{
			this->shape_points_number = n;
		}

		bool isClustered()
		{
			if (cluster_labels_set.rows > 1)
			{
				return true;
			}
			else
			{
				return false;
			}
		}

		int getProfileDimension()
		{
			return this->profile_dimension;
		}

		void push_back( BClassType label, PointMatrixType & point, int cluster_label = 0 )
		{
			KNNUnit * unit = this->getKNNUnitByClusterLabel(cluster_label); //We push back every sample into the first clustered classifier unit by default.

			unit->push_back(static_cast<int>(label), point);
		}

		KNNUnit* getKNNUnitByPointIndex( int ptidx )
		{
			if (ptidx<0 || ptidx>=this->shape_points_number)
			{
				std::cerr<<"Cannot find KNN unit for point index: "<<ptidx<<std::endl;
				return NULL;
			}

			int cl = this->cluster_labels[ptidx][0];
			return this->getKNNUnitByClusterLabel(cl);
		}

		KNNUnit* getKNNUnitByClusterLabel( int cluster_label )
		{
			KNNUnitMapType::iterator unitIt = knnUnitMap.find( cluster_label );
			
			if ( unitIt == knnUnitMap.end() )
			{
				std::cerr<<"Cannot find KNNUnit for label "<<cluster_label<<std::endl;
			}

			return unitIt->second;
		}

		//Don't call this.
		void initilize( int profile_count_, int profile_dimension_ )
		{
			if (shape_number<=0)
			{
				std::cout<<"Please set shape number before calling initialization."<<std::endl;
				return;
			}

			if (shape_points_number<=0)
			{
				std::cout<<"Please set shape points number before calling initialization."<<std::endl;
				return;
			}

			this->profile_dimension = profile_dimension_;

			DatasetMatrixType datasetDefault = DatasetMatrixType(new float[profile_count_*profile_dimension], profile_count_, profile_dimension);
			LabelMatrixType   labelsDefault  = LabelMatrixType(new int[profile_count_], profile_count_, 1);
			KNNUnit* knnUnitDefault = new KNNUnit( datasetDefault, labelsDefault, 0 );
			knnUnitMap.insert( KNNUnitPairType( 0, knnUnitDefault ) ); //0 is the default kmeans label;

			cluster_labels        = IndexMatrixType(new int[shape_points_number], shape_points_number, 1);
			fillFlannMatrix<IndexMatrixType>( cluster_labels, 0 );

			cluster_labels_set    = IndexMatrixType(new int[1], 1, 1); //By default, there is only 1 cluster.
			fillFlannMatrix<IndexMatrixType>( cluster_labels_set, 0 );

			shape_number_matrix = IntegerMatrixType(new int[1], 1, 1);
			fillFlannMatrix<IntegerMatrixType>( shape_number_matrix, this->shape_number );

			shape_points_number_matrix = IntegerMatrixType( new int[1], 1, 1 );
			fillFlannMatrix<IntegerMatrixType>( shape_points_number_matrix, this->shape_points_number );

			is_initilized = true;
		}

		void buildIndex()
		{
			KNNUnitMapType::iterator it = knnUnitMap.begin();
			KNNUnitMapType::iterator itEnd = knnUnitMap.end();

			while(it!=itEnd)
			{
				std::cout<<"Start to build index for clustered classifier: "<<it->first<<std::endl;
				it->second->buildIndex();
				it++;
			}
		}

		void getDatasetAndLabels( DatasetMatrixType& dataset, LabelMatrixType& labels, int cluster_label = 0 )
		{
			KNNUnit * unit = this->getKNNUnitByClusterLabel(cluster_label);

			dataset = unit->getDataset();
			labels  = unit->getLabels();
		}

		void getClusterLabelAndLabelSet( IndexMatrixType & _cluster_labels, IntegerMatrixType & _cluster_labels_set )
		{
			_cluster_labels = this->cluster_labels;
			_cluster_labels_set = this->cluster_labels_set;
		}

		void loadSamples( const char* filename )
		{
			std::cout<<"Load from sample file: "<<filename<<std::endl;

			LONGTYPE profile_count_ = 0;
			int      profile_dimension_ = 0;

			//Read profile count and profile dimension
			{
				string line;
				ifstream myfile (filename);
				if (myfile.is_open())
				{
					while ( getline (myfile,line) )
					{
						if ( line.find_first_not_of(' ') != std::string::npos )
						{
							if (profile_count_ == 0)
							{
								istringstream s;
								s.str(line);
								double t;
								while (s >> t)
								{
									if (!s.good())
									{
										break;
									}
									else
									{
										profile_dimension_ ++;
									}
								}
							}
							profile_count_++;
						}
					}
					myfile.close();
				}
				else
				{
					std::cout << "Unable to open file: " << filename << std::endl;
					return;
				}
			}
			profile_dimension_ -= 1;
			std::cout<<"profile_count: "<<profile_count_<<std::endl; 
			std::cout<<"profile_dimension: "<<profile_dimension_<<std::endl;

			profile_dimension = profile_dimension_;

			this->initilize(profile_count_, profile_dimension_);

			 BClassType label;
			 DatasetMatrixType data_single(new float[profile_dimension_], 1, profile_dimension_);

			//Read data and label
			string line;
			ifstream myfile (filename);

			LONGTYPE row_num = 0;  //useless
			int      col_num = 0;

			if (myfile.is_open())
			{
				while ( getline (myfile,line) )
				{
					if ( line.find_first_not_of(' ') != std::string::npos )
					{
						istringstream s;
						s.str(line);
						float t = 0.0;
						col_num = -1;
						while (s >> t)
						{
							if (!s.good())
							{
								break;
							}
							else
							{
								if (col_num==-1)
								{
									label = BClassType(t);
								}
								else
								{
									if (col_num>profile_dimension_-1)
									{
										std::cout<<"Wrong data line: "<<col_num<<" "<<row_num<<" "<<t<<std::endl;
									}
									data_single[0][col_num] = t;
								}
								col_num++;
							}
						}

						this->push_back(label, data_single, 0);
						
						row_num++;
					}
				}
				myfile.close();
			}
			else
			{
				std::cout << "Unable to open file: " << filename << std::endl;
				return;
			}

			std::cout<<"Load from sample file done"<<std::endl;

			delete[] data_single.ptr();
		}

		void autoWeight()
		{
			if ( isClustered() )
			{
				std::cout<<"Has been clustered. Automatic weighting can only be executed before clustering."<<std::endl;
				return;
			}
			
			KNNUnit* defaultUnit = this->getKNNUnitByClusterLabel(0);
			
			itk::VariableLengthVector<double> weights;
			itk::VariableLengthVector<double> means;
			weights.SetSize(defaultUnit->getProfileDimension());
			weights.Fill(1.0);
			means.SetSize(defaultUnit->getProfileDimension());
			means.Fill(0.0);
			double maxMean = 0.0;
			
			for (int c=0;c<defaultUnit->getProfileDimension();c++)
			{
				for (unsigned long r=0;r<defaultUnit->getPointsCount();r++)
				{
					means[c] += std::abs(defaultUnit->dataset[r][c]);
				}
				means[c] /= defaultUnit->getPointsCount();
				if (means[c]>maxMean)
				{
					maxMean = means[c];
				}
			}

			for (int c=0;c<defaultUnit->getProfileDimension();c++)
			{
				weights[c] = maxMean / means[c];
				if (weights[c] != 1.0)
				{
					for (unsigned long r=0;r<defaultUnit->getPointsCount();r++)
					{
						defaultUnit->dataset[r][c] *= weights[c];
					}
				}
			}

			std::cout<<"[AutoWeight] Means: "<<means<<std::endl;
			std::cout<<"[AutoWeight] Weights: "<<weights<<std::endl;
		}

		void load( const char* filename )
		{
			if (is_initilized)
			{
				std::cout<<"Has been loaded.."<<std::endl;
				return;
			}

			try
			{
				flann::load_from_file(cluster_labels,  filename, "cluster_labels");
				flann::load_from_file(cluster_labels_set,  filename, "cluster_labels_set");

				flann::load_from_file(shape_number_matrix,  filename, "shape_number");
				this->shape_number = shape_number_matrix[0][0];
				
				flann::load_from_file(shape_points_number_matrix,  filename, "shape_points_number");
				this->shape_points_number = shape_points_number_matrix[0][0];

				//FIXEME Need to check if the number of cluster is valid.
				for (int i=0;i<cluster_labels_set.rows;i++)
				{
					DatasetMatrixType datasetTmp;
					LabelMatrixType   labelsTmp;

					int cl = cluster_labels_set[i][0];

					std::stringstream sstm_dataset;
					sstm_dataset << "dataset_" << cl;
					flann::load_from_file(datasetTmp, filename, sstm_dataset.str());

					std::stringstream sstm_labels;
					sstm_labels << "labels_" << cl;
					flann::load_from_file(labelsTmp, filename, sstm_labels.str());

					KNNUnit* unit = new KNNUnit( datasetTmp, labelsTmp, cl );
					knnUnitMap.insert( KNNUnitPairType(cl, unit) );

					if (profile_dimension == 0)
					{
						profile_dimension = datasetTmp.cols;
					}
				}
			}
			catch(...)
			{
				std::cout<<"WARNING: exception thrown when load from HDF5 file."<<std::endl;
			}

			is_initilized = true;

			this->buildIndex();
		}

		void save( const char* filename )
		{
			std::cout<<"Start to save profiles..."<<std::endl;
			KNNUnitMapType::iterator it = knnUnitMap.begin();
			KNNUnitMapType::iterator itEnd = knnUnitMap.end();

			while(it!=itEnd)
			{
				int cl = it->first;

				DatasetMatrixType datasetTmp;
				LabelMatrixType   labelsTmp;
				this->getDatasetAndLabels( datasetTmp, labelsTmp, cl );

				std::stringstream sstm_dataset;
				sstm_dataset << "dataset_" << cl;
				flann::save_to_file( datasetTmp , filename, sstm_dataset.str() );

				std::stringstream sstm_labels;
				sstm_labels << "labels_" << cl;
				flann::save_to_file( labelsTmp, filename, sstm_labels.str() );

				it++;
			}

			flann::save_to_file(cluster_labels, filename, "cluster_labels");
			flann::save_to_file(cluster_labels_set, filename, "cluster_labels_set");
			this->shape_number_matrix[0][0] = this->shape_number;
			flann::save_to_file(shape_number_matrix, filename, "shape_number");
			this->shape_points_number_matrix[0][0] = this->shape_points_number;
			flann::save_to_file(shape_points_number_matrix, filename, "shape_points_number");
		}

		void test( double ratio_testing )
		{
			if (!is_initilized)
			{
				std::cout<<"Has not loaded."<<std::endl;
				return;
			}

			std::cout<<"Start to test classifying of ratio "<<ratio_testing<<"..."<<std::endl;

			double ratio_weighted = 0.0;
			unsigned long pt_count_sum = 0;

			for (int i=0;i<cluster_labels_set.rows;i++)
			{
				int cl = cluster_labels_set[i][0];
				KNNUnit* unit = this->getKNNUnitByClusterLabel(cl);

				double success_ratio = unit->test( ratio_testing, this->max_neighbour );

				std::cout<<"Label->"<<cl<<"; Success Ratio->"<<success_ratio<<std::endl;

				pt_count_sum += unit->getPointsCount();
				ratio_weighted += success_ratio*unit->getPointsCount();
			}

			ratio_weighted /= pt_count_sum;

			std::cout<<"Success Ratio Weighted: "<<ratio_weighted<<std::endl;
			std::cout<<"Test done!!!"<<std::endl;
		}

		void classify( PointMatrixType& point, int ptidx, ProbabilityContainer & probabilities, DistanceContainer & dists )
		{
			KNNUnit* unit = this->getKNNUnitByPointIndex( ptidx );

			unit->classify(point, probabilities, dists, this->max_neighbour);
		}

		void classify( PointMatrixType& point, int ptidx, BClassType& result_class, ProbabilityContainer & probabilities, DistanceContainer & dists )
		{
			this->classify(point,ptidx,probabilities,dists);

			BClassType classifiedLabel;
			if (probabilities[OPClass]>=probabilities[BPClass] && probabilities[OPClass]>=probabilities[IPClass])
			{
				classifiedLabel = OPClass;
			}
			else if (probabilities[BPClass]>=probabilities[OPClass] && probabilities[BPClass]>=probabilities[IPClass])
			{
				classifiedLabel = BPClass;
			}
			else
			{
				classifiedLabel = IPClass;
			}

			result_class = classifiedLabel;
		}

		void cluster( int cluster_number_request )
		{
			if (!is_initilized)
			{
				std::cout<<"Has not initialized. Cannot be clustered."<<std::endl;
				return;
			}

			if ( isClustered() )
			{
				std::cout<<"Has been clustered. No need to cluster again."<<std::endl;
				return;
			}

			KNNUnit* unitDefault = knnUnitMap[0];

			int numberOfSamplesPerLandmark = (unitDefault->getPointsCount()/shape_number)/shape_points_number;
			std::cout<<"number of samples per landmark: "<<numberOfSamplesPerLandmark<<std::endl;

			const unsigned int measureLength = profile_dimension*shape_number*numberOfSamplesPerLandmark;
			const unsigned int measureCount = shape_points_number;

			std::cout<<"measure length: "<<measureLength<<std::endl;
			std::cout<<"measure count:"<<measureCount<<std::endl;

			typedef itk::VariableLengthVector< double> MeasurementVectorType;
			typedef itk::Statistics::ListSample< MeasurementVectorType > SampleType;
			SampleType::Pointer sample = SampleType::New();
			sample->SetMeasurementVectorSize( measureLength );

			MeasurementVectorType mv;
			mv.SetSize( measureLength );
			mv.Fill( 0 );
			for (int i=0;i<shape_points_number;i++)
			{
				for (int j=0;j<shape_number;j++)
				{
					for (int k=0;k<numberOfSamplesPerLandmark;k++)
					{
						for (int t=0;t<profile_dimension;t++)
						{
							mv[ j*numberOfSamplesPerLandmark*profile_dimension + k*profile_dimension + t ] = unitDefault->dataset[j*shape_points_number*numberOfSamplesPerLandmark + i*numberOfSamplesPerLandmark + k][t];
						}
					}
				}

				//std::cout<<mv<<std::endl;

				//std::cout<<mv<<std::endl;
				sample->PushBack( mv );
			}

			typedef itk::Statistics::WeightedCentroidKdTreeGenerator< SampleType > TreeGeneratorType;
			TreeGeneratorType::Pointer treeGenerator = TreeGeneratorType::New();

			//std::cout << "TTTT" << std::endl;

			treeGenerator->SetSample( sample );
			treeGenerator->SetBucketSize( 16 );
			treeGenerator->Update();

			typedef TreeGeneratorType::KdTreeType TreeType;
			typedef itk::Statistics::KdTreeBasedKmeansEstimator<TreeType> EstimatorType;
			EstimatorType::Pointer estimator = EstimatorType::New();

			EstimatorType::ParametersType initialMeans( measureLength * cluster_number_request );
			initialMeans.Fill(0.0);

			//std::cout << "TTTT" << std::endl;

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

			//std::cout << "TTTT" << std::endl;

			int index = 0;
			for ( unsigned int i = 0 ; i < cluster_number_request ; i++ )
			{
				MembershipFunctionType::Pointer membershipFunction = MembershipFunctionType::New();
				MembershipFunctionType::CentroidType centroid( sample->GetMeasurementVectorSize() );
				for ( unsigned int j = 0 ; j < sample->GetMeasurementVectorSize(); j++ )
				{
					centroid[j] = estimatedMeans[index++];
				}
				membershipFunction->SetCentroid( centroid );
				membershipFunctionVector.push_back( membershipFunction.GetPointer() );
			}
			classifier->SetMembershipFunctions( membershipFunctionVectorObject );
			classifier->Update();

			//std::cout << "TTTT" << std::endl;

			const ClassifierType::MembershipSampleType* membershipSample = classifier->GetOutput();
			ClassifierType::MembershipSampleType::ConstIterator labelIt = membershipSample->Begin();
			ClassifierType::MembershipSampleType::ConstIterator labelItEnd = membershipSample->End();

			//Initial cluster map.
			typedef std::map<int, int> ClusterCountMapType;
			ClusterCountMapType clusterCountMap;//<cluster_label, point_count>
			for (int i=0;i<cluster_number_request;i++)
			{
				int cl = i;
				clusterCountMap.insert(std::pair<int, int>(cl, 0));
			}

			//Count points number of each cluster.
			int idx = 0;
			labelIt = membershipSample->Begin();
			while(labelIt!=labelItEnd)
			{
				int cl = static_cast<int>(labelIt.GetClassLabel());
				cluster_labels[idx++][0] = cl;
				clusterCountMap[cl]++;
				++labelIt;
			}
			
			//Shrink cluster map.
			for (int i=0;i<cluster_number_request;i++)
			{
				if (clusterCountMap[i] == 0) {
					clusterCountMap.erase ( i );
				}
			}

			int cluster_number_shrinked = clusterCountMap.size();
			std::cout<<"cluster number shrank: "<<cluster_number_shrinked<<std::endl;

			cluster_labels_set = IndexMatrixType(new int[cluster_number_shrinked], cluster_number_shrinked, 1);

			//Calculate profile number for each cluster for allocating.
			idx = 0;
			for (ClusterCountMapType::iterator it=clusterCountMap.begin(); it!=clusterCountMap.end(); ++it)
			{
				int cl = it->first;
				unsigned long ptcount = it->second;

				cluster_labels_set[idx++][0] = cl;

				std::cout<<"Profile count for cluster label "<<cl<<": "<<ptcount<<std::endl;

				unsigned long count = ptcount * this->shape_number * numberOfSamplesPerLandmark;

				DatasetMatrixType datasetTmp(new float[count*profile_dimension], count, profile_dimension);
				LabelMatrixType   labelsTmp(new int[count], count, 1);

				KNNUnit* unit = new KNNUnit( datasetTmp, labelsTmp, cl );
				knnUnitMap[cl] = unit;
			}

			PointMatrixType ptmatrix = PointMatrixType(new float[profile_dimension], 1, profile_dimension);

			for (unsigned long r=0;r<unitDefault->dataset.rows && r<unitDefault->labels.rows; r++)
			{
				int label = unitDefault->labels[r][0];
				for (int c=0;c<profile_dimension;c++)
				{
					ptmatrix[0][c] = unitDefault->dataset[r][c];
				}

				int cl = this->cluster_labels[r%(shape_points_number*numberOfSamplesPerLandmark)/numberOfSamplesPerLandmark][0];
				this->push_back( label, ptmatrix, cl );
			}

			//this->buildIndex();

			//DELETE_NULLABLE_PRT(ptmatrix);
			delete unitDefault;

			std::cout<<"Cluster done!!"<<std::endl;
		}

		void copyCluster(km::KNNProfileClassifier & classiferSource)
		{
			const int cluster_number = classiferSource.cluster_labels_set.rows;

			//Initial cluster map.
			typedef std::map<int, int> ClusterCountMapType;
			ClusterCountMapType clusterCountMap;//<cluster_label, point_count>

			cluster_labels_set = IndexMatrixType(new int[cluster_number], cluster_number, 1);

			//Copy cluster label set.
			for (int i=0;i<cluster_number;i++)
			{
				int cl = classiferSource.cluster_labels_set[i][0];
				this->cluster_labels_set[i][0] = cl;
				clusterCountMap.insert(std::pair<int, int>(cl, 0));
			}

			//Copy cluster labels.
			for (int i=0;i<classiferSource.cluster_labels.rows;i++)
			{
				int cl = classiferSource.cluster_labels[i][0];
				this->cluster_labels[i][0] = cl;
				clusterCountMap[cl]++;
			}

			KNNUnit* unitDefault = knnUnitMap[0];

			int numberOfSamplesPerLandmark = (unitDefault->getPointsCount()/shape_number)/shape_points_number;

			//Calculate profile number for each cluster for allocating.
			for ( int i=0;i<cluster_number;i++ )
			{
				int cl = classiferSource.cluster_labels_set[i][0]; //cluster label

				unsigned long profileCount = clusterCountMap[cl] * this->shape_number * numberOfSamplesPerLandmark;;

				std::cout<<"Profile count for cluster label "<<cl<<": "<<profileCount<<std::endl;

				DatasetMatrixType datasetTmp(new float[profileCount*profile_dimension], profileCount, profile_dimension);
				LabelMatrixType   labelsTmp(new int[profileCount], profileCount, 1);

				KNNUnit* unit = new KNNUnit( datasetTmp, labelsTmp, cl );
				knnUnitMap[cl] = unit;
			}

			PointMatrixType ptmatrix = PointMatrixType(new float[profile_dimension], 1, profile_dimension);

			for (unsigned long r=0;r<unitDefault->dataset.rows && r<unitDefault->labels.rows; r++)
			{
				int label = unitDefault->labels[r][0];
				for (int c=0;c<profile_dimension;c++)
				{
					ptmatrix[0][c] = unitDefault->dataset[r][c];
				}

				int cl = this->cluster_labels[r%(shape_points_number*numberOfSamplesPerLandmark)/numberOfSamplesPerLandmark][0];
				this->push_back( label, ptmatrix, cl );
			}

			//this->buildIndex();

			//DELETE_NULLABLE_PRT(ptmatrix);
			delete unitDefault;

			std::cout<<"Copy cluster done!!"<<std::endl;
		}

		void print()
		{
			std::cout<<"*********************KNNProfileClassifier***************************"<<std::endl;
			std::cout<<"Number of shape: "<<this->shape_number<<std::endl;
			std::cout<<"Number of points on each shape: "<<this->shape_points_number<<std::endl;
			std::cout<<"K-Means clustered: "<< (this->isClustered()?"YES":"NO") <<std::endl;
			std::cout<<"Maximum neighbor for KNN search: "<<this->max_neighbour<<std::endl;
			std::cout<<"********************************************************************"<<std::endl;
		}
	};
	
	/************************************************************************************************/
	// Adaboost profile classifier
	/************************************************************************************************/
	class ProfileClassifier
	{
	public:
		class AdaboostUnit
		{
		public:
			int cluster_label; //label for k means clustering.
			PROFILE_CATEGORY profile_category;

			vector<int> sign;
			vector<FLOATTYPE> alpha;
			vector<FLOATTYPE> threshold;
			vector<int> featureID;
			vector<FLOATTYPE> error;

			DatasetMatrixType alpha_set;
			IntegerMatrixType sign_set;
			DatasetMatrixType threshold_set;
			IndexMatrixType   featureID_set;
			DatasetMatrixType error_set;

			AdaboostUnit( PROFILE_CATEGORY category )
			{
				cluster_label = 0;
				profile_category = category;
			}

			~AdaboostUnit()
			{
				std::cout<<"De-constructor of AdaboostUnit.."<<std::endl;
				DELETE_NULLABLE_PRT(alpha_set);
				DELETE_NULLABLE_PRT(sign_set);
				DELETE_NULLABLE_PRT(threshold_set);
				DELETE_NULLABLE_PRT(featureID_set);
				std::cout<<"De-constructor end.."<<std::endl;
			}

			int label2Binary( BClassType c )
			{
				if (profile_category == BOUNDARY)
				{
					if (c == BPClass)
					{
						return 1;
					}
					else
					{
						return 0;
					}
				}
				else if(profile_category == LIVER)
				{
					if (c == IPClass)
					{
						return 1;
					}
					else
					{
						return 0;
					}
				}
				else
				{
					if (c == OPClass)
					{
						return 1;
					}
					else
					{
						return 0;
					}
				}
			}

			void FillMatrixToVector()
			{
				for(int i=0;i<alpha_set.rows;i++)
				{
					this->alpha.push_back(this->alpha_set[i][0]);
					this->sign.push_back(this->sign_set[i][0]);
					this->threshold.push_back(this->threshold_set[i][0]);
					this->featureID.push_back(this->featureID_set[i][0]);
					this->error.push_back(this->error_set[i][0]);
				}
			}

			void FillVectorToMatrix()
			{
				for(int i=0;i<alpha.size();i++)
				{
					this->alpha_set[i][0] = this->alpha[i];
					this->sign_set[i][0] = this->sign[i];
					this->threshold_set[i][0] = this->threshold[i];
					this->featureID_set[i][0] = this->featureID[i];
					this->error_set[i][0] = this->error[i];
				}
			}

			//Return 1: positive case. The classification success;
			//Return 0: negtive case. The classification failed.
			int check_positive( DatasetMatrixType & samples, LabelMatrixType & labels, LONGTYPE row_id)
			{
				BClassType label = (BClassType)labels[row_id][0];
				FLOATTYPE p = this->classify(samples, row_id);

				if ((p>=0.5 && label2Binary(label)==1) || (p<0.5 && label2Binary(label)==0) )
				{
					return 1;
				}
				else
				{
					return 0;
				}
			}

			double test(DatasetMatrixType & samples, LabelMatrixType & labels)
			{
				assert( samples.rows == labels.rows );

				const LONGTYPE NSample = samples.rows;
				LONGTYPE NSuccess = 0;

				for (LONGTYPE r=0;r<samples.rows;r++)
				{
					int positive_count = check_positive(samples, labels, r);
					NSuccess += positive_count;
				}

				double successRatio = static_cast<double>(NSuccess)/NSample;

				return successRatio;
			}

			FLOATTYPE classify( const std::vector<FLOATTYPE> & sample )
			{
				const int NFeature = sample.size();
				int LC = this->alpha.size();

				std::vector<FLOATTYPE> X;
				for (int i=0;i<sample.size();i++)
				{
					X.push_back(sample[i]);
				}

				FLOATTYPE H,cH;
				AdaBoostClassify(&X[0], 1, NFeature, LC, &featureID[0], &alpha[0], &sign[0], &threshold[0], &H, &cH);

				return H;
			}			

			FLOATTYPE classify( DatasetMatrixType & sample, LONGTYPE row_id )
			{
				std::vector<FLOATTYPE> sample_vec;
				for (int c=0;c<sample.cols;c++)
				{
					sample_vec.push_back( sample[row_id][c] );
				}

				FLOATTYPE H = this->classify(sample_vec);
				return H;
			}
		};

		typedef std::vector<FLOATTYPE> AdaFeatureContainer;
		typedef std::vector<char>      AdaLabelContainer;
		typedef std::pair<int, AdaboostUnit*> AdaboostUnitPairType;
		typedef std::map<int, AdaboostUnit*> AdaboostUnitMapType;
		typedef std::map<int, float> IntFloatMapType;
		
		PROFILE_CATEGORY profile_category;

		IndexMatrixType cluster_labels;
		IntegerMatrixType cluster_labels_set;
		
		AdaboostUnitMapType adaboostUnitMap;

		ProfileClassifier( PROFILE_CATEGORY category )
		{
			profile_category = category;

			cluster_labels_set    = IndexMatrixType(new int[1], 1, 1); //By default, there is only 1 cluster.
			fillFlannMatrix<IndexMatrixType>( cluster_labels_set, 0 );
		}

		~ProfileClassifier()
		{
			std::cout<<"De-constructor of ProfileClassifier.."<<std::endl;
			//DELETE_NULLABLE_PRT(cluster_labels);
			//DELETE_NULLABLE_PRT(cluster_labels_set);
			std::cout<<"De-constructor end.."<<std::endl;
		}

		int getNumberOfClusters()
		{
			if (cluster_labels_set.rows > 0)
			{
				return cluster_labels_set.rows;
			}
			else
			{
				return 0;
			}
		}

		AdaboostUnit* getAdaboostUnitByPointIndex( int ptidx )
		{
			int cl = 0;

			if (this->cluster_labels_set.rows <= 1)
			{
				cl = 0;
			}
			else if (ptidx<0 || ptidx>=this->cluster_labels.rows)
			{
				std::cerr<<"Cannot find Adaboost unit for point index: "<<ptidx<<std::endl;
				return NULL;
			}
			else
			{
				cl = this->cluster_labels[ptidx][0];
			}

			return this->getAdaboostUnitByClusterLabel(cl);
		}

		AdaboostUnit* getAdaboostUnitByClusterLabel( int cluster_label )
		{
			AdaboostUnitMapType::iterator unitIt = adaboostUnitMap.find( cluster_label );

			if ( unitIt == adaboostUnitMap.end() )
			{
				std::cerr<<"Cannot find Adaboost unit for label "<<cluster_label<<std::endl;
				unitIt = adaboostUnitMap.find(0);
			}

			return unitIt->second;
		}

		double classify( const std::vector<FLOATTYPE> & sample, int ptidx = 0 )
		{
			AdaboostUnit* adaboostUnit = this->getAdaboostUnitByPointIndex(ptidx);
			if (adaboostUnit==NULL)
			{
				std::cerr<<"Cannnot find any Adaboost unit for point: "<<ptidx<<std::endl;
				return 0;
			}
			return adaboostUnit->classify(sample);
		}

		double classify( DatasetMatrixType & sample, LONGTYPE row_id, int ptidx = 0 )
		{
			AdaboostUnit* adaboostUnit = this->getAdaboostUnitByPointIndex(ptidx);
			if (adaboostUnit==NULL)
			{
				std::cerr<<"Cannnot find any Adaboost unit for point: "<<ptidx<<std::endl;
				return 0;
			}
			return adaboostUnit->classify(sample, row_id);
		}

		void train(KNNProfileClassifier & classifier, const char* outputdir = NULL)
		{
			std::cout<<"***Label set: ";
			for (int tt=0;tt<classifier.cluster_labels_set.rows;tt++)
			{
				std::cout<<classifier.cluster_labels_set[tt][0]<<" ";
			}
			std::cout<<std::endl;

			classifier.getClusterLabelAndLabelSet( cluster_labels, cluster_labels_set );

			std::map<int, std::string> adaboostfilenames;

			std::cout<<"Number of clusters: "<<cluster_labels_set.rows<<std::endl;

			for (int tt=0;tt<cluster_labels_set.rows;tt++)
			{
				int cl = cluster_labels_set[tt][0];
				std::cout<<"Start to train classifer for label "<<cl<<std::endl;
				AdaboostUnit * unit = new AdaboostUnit(profile_category);
				unit->cluster_label = cl;
				adaboostUnitMap[cl] = unit;

				DatasetMatrixType tmpdataset;
				LabelMatrixType tmplabelset;
				classifier.getDatasetAndLabels(tmpdataset, tmplabelset, cl);
				
				assert( tmpdataset.rows == tmplabelset.rows );

				int NFeature = tmpdataset.cols;
				unsigned long long NSample = tmpdataset.cols * tmpdataset.rows;

				AdaFeatureContainer X;
				AdaLabelContainer Y;
				for (unsigned long r=0;r<tmpdataset.rows;r++)
				{
					for (unsigned long c=0;c<tmpdataset.cols;c++)
					{
						X.push_back( static_cast<FLOATTYPE>(tmpdataset[r][c]) );
					}

					Y.push_back(unit->label2Binary(tmplabelset[r][0]));
				}

				char adaboostFile[1024];
				if (outputdir != NULL)
				{
					sprintf (adaboostFile, "%s/AdaBoostResult_%d",outputdir, cl);
				}
				else
				{
					sprintf (adaboostFile, "AdaBoostResult_%d", cl);
				}

				std::cout<<adaboostFile<<std::endl;

				try
				{
					//AdaBoostTrain(&X[0],&Y[0],Y.size(),X.size()/Y.size(), 10, ss.str().c_str());
					AdaBoostTrain(&X[0],&Y[0],Y.size(),X.size()/Y.size(), 100, adaboostFile);
				}
				catch(const exception &e)
				{
					std::cout<<e.what()<<std::endl;
				}

				this->readTrainFile(adaboostFile, unit, 50);
			}
		}

		void readTrainFile( const char* adaboostFile, AdaboostUnit* unit, int maxWeakClassifier = 100)
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
				unit->alpha.push_back(t);
				ifs >> t;
				ifs >> t;
				unit->featureID.push_back(int(t));
				ifs >> t;
				unit->sign.push_back(int(t));
				ifs >> t;
				unit->threshold.push_back(t);
				ifs >> t;
				unit->error.push_back(t);

				weak_classifier_number++;
			}
			ifs.close();
			LC = unit->alpha.size();
			cout<<"# weak learners: "<<LC<<endl;

			unit->alpha_set = DatasetMatrixType(new float[LC], LC, 1);
			unit->sign_set  = IntegerMatrixType(new int[LC], LC, 1);
			unit->threshold_set = DatasetMatrixType(new float[LC], LC, 1);
			unit->featureID_set = IndexMatrixType(new int[LC], LC, 1);
			unit->error_set = DatasetMatrixType(new float[LC], LC, 1);
		}

		void save(const char* filename)
		{
			AdaboostUnitMapType::iterator it = adaboostUnitMap.begin();
			AdaboostUnitMapType::iterator itEnd = adaboostUnitMap.end();

			while(it!=itEnd)
			{
				int cl = it->first;
				AdaboostUnit* unit = it->second;

				unit->FillVectorToMatrix();

				std::stringstream ss;
				ss << "alpha_set_" << cl;
				flann::save_to_file( unit->alpha_set , filename, ss.str() );

				ss.str("");
				ss << "sign_set_" << cl;
				flann::save_to_file( unit->sign_set, filename, ss.str() );

				ss.str("");
				ss << "threshold_set_" << cl;
				flann::save_to_file( unit->threshold_set, filename, ss.str() );

				ss.str("");
				ss << "featureID_set_" << cl;
				flann::save_to_file( unit->featureID_set, filename, ss.str() );

				ss.str("");
				ss << "error_set_" << cl;
				flann::save_to_file( unit->error_set, filename, ss.str() );

				it++;
			}

			flann::save_to_file(cluster_labels, filename, "cluster_labels");
			flann::save_to_file(cluster_labels_set, filename, "cluster_labels_set");

			std::cout<<"==========>Save adaboost classifier successfully!<============="<<std::endl;
		}

		void load(const char* filename)
		{
			try
			{
				flann::load_from_file(cluster_labels,  filename, "cluster_labels");
				flann::load_from_file(cluster_labels_set,  filename, "cluster_labels_set");

				//FIXEME Need to check if the number of cluster is valid.
				for (int i=0;i<cluster_labels_set.rows;i++)
				{
					int cl = cluster_labels_set[i][0];

					AdaboostUnit* unit = new AdaboostUnit(profile_category);
					unit->cluster_label = cl;

					std::stringstream ss;
					ss << "alpha_set_" << cl;
					flann::load_from_file(unit->alpha_set, filename, ss.str());

					ss.str("");
					ss << "sign_set_" << cl;
					flann::load_from_file(unit->sign_set, filename, ss.str());

					ss.str("");
					ss << "threshold_set_" << cl;
					flann::load_from_file(unit->threshold_set, filename, ss.str());

					ss.str("");
					ss << "featureID_set_" << cl;
					flann::load_from_file(unit->featureID_set, filename, ss.str());

					ss.str("");
					ss << "error_set_" << cl;
					flann::load_from_file(unit->error_set, filename, ss.str());

					unit->FillMatrixToVector();

					adaboostUnitMap[unit->cluster_label] = unit;
				}

				std::cout<<"==========>Load adaboost classifier successfully!<============="<<std::endl;
			}
			catch(...)
			{
				std::cout<<"WARNING: exception thrown when load from HDF5 file."<<std::endl;
			}
		}

		void test(KNNProfileClassifier & classifier, IntFloatMapType & errorMap)
		{
			double ratio_weighted = 0.0;
			unsigned long pt_count_sum = 0;

			//Initialize error map.
			for (int i=0;i<cluster_labels_set.rows;i++)
			{
				errorMap[cluster_labels_set[i][0]] = 0;
			}

			for (int i=0;i<cluster_labels_set.rows;i++)
			{
				int cl = cluster_labels_set[i][0];
				AdaboostUnit* adaboostUnit = this->getAdaboostUnitByClusterLabel(cl);

				KNNProfileClassifier::KNNUnit* knnUnit = classifier.getKNNUnitByClusterLabel(cl);

				double success_ratio = adaboostUnit->test( knnUnit->dataset, knnUnit->labels );

				std::cout<<"Label->"<<cl<<"; Success Ratio->"<<success_ratio<<std::endl;

				pt_count_sum += knnUnit->getPointsCount();
				ratio_weighted += success_ratio*knnUnit->getPointsCount();

				errorMap[cl] = success_ratio;
			}

			ratio_weighted /= pt_count_sum;

			std::cout<<"Success Ratio Weighted: "<<ratio_weighted<<std::endl;
			std::cout<<"Test done!!!"<<std::endl;	
		}
		
		void test(KNNProfileClassifier & classifier)
		{
			IntFloatMapType errorMap;
			test(classifier, errorMap);
		}
	};
}

#include "kmProfileClassifier.hxx"

#undef DELETE_NULLABLE_PRT
#undef MAX_NEIGHBOUR
#endif