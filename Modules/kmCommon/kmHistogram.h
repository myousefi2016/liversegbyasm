#ifndef __kmHistogram_h
#define __kmHistogram_h

namespace km
{
	struct HistogramNode 
	{
		HistogramNode():lower(0), upper(0), percentage(0), count(0)
		{
		};
		double lower;
		double upper;
		double percentage;
		unsigned count;
	};

	struct Histogram
	{
		double statisticMin;
		double statisticMax;
		unsigned int binNum;
		double binWidth;
		HistogramNode * nodes;
		Histogram(double _staMin, double _staMax, double _binWidth)
		{
			statisticMin = _staMin;
			statisticMax = _staMax;
			binWidth = _binWidth;
			binNum = static_cast<unsigned int>(std::ceil((statisticMax-statisticMin)/binWidth)) + 1;
			nodes = new HistogramNode[binNum];
		}

		~Histogram()
		{
			delete[] nodes;
		}

		void clear()
		{
			for (int i=0;i<binNum;i++){
				this->nodes[i].count = 0;
				this->nodes[i].percentage = 0.0;
			}
		}

		void calcThresholdOfRatio( double ratio, double & thresholdLower, double & thresholdUpper )
		{
			if (ratio<=0){
				return;
			}
			double totalPercentage = 0;
			double max=this->nodes[0].percentage;
			int binNumM;
			for(int i=0;i<this->binNum;i++){
				if (this->nodes[i].percentage>max){
					max=this->nodes[i].percentage;
					binNumM = i;
				}
			}
			double sum=max;
			int i=binNumM;
			int j=binNumM;
			int flag=0;
			bool foundflag = false;
			while(i>=0&&j<this->binNum){
				if (flag==0){
					i--;
					flag=1;
					sum+=this->nodes[i].percentage;
				}else{
					j++; 
					flag=0; 
					sum+=this->nodes[j].percentage;
				}
				
				if (sum>ratio){
					thresholdLower = this->nodes[i].lower;
					thresholdUpper = this->nodes[j].upper; 
					foundflag = true;
					break;
				}
			}
			if (!foundflag){
				if(i<0) i+=1;
				if(j>binNumM) j-=1;
				thresholdLower = this->nodes[i].lower;
				thresholdUpper = this->nodes[j].upper;
			}
		}

		double calcRatioOfThreshold( double thresholdLower, double thresholdUpper )
		{
			thresholdLower = std::max(thresholdLower, this->statisticMin);
			thresholdUpper = std::min(thresholdUpper, this->statisticMax);

			if(thresholdUpper < thresholdLower){
				return 0.0;
			}
			unsigned int binIdxLower, binIdxUpper;
			binIdxLower = static_cast<unsigned int>(std::floor( (thresholdLower - this->statisticMin)/binWidth ));
			binIdxUpper = static_cast<unsigned int>(std::ceil(  (thresholdUpper - this->statisticMin)/binWidth ));
			double ratio = 0;
			for (unsigned int i=binIdxLower;i<binIdxUpper;i++){
				ratio += this->nodes[i].percentage;
			}
			return ratio;
		}
	};
	
	//Calculate image histogram
	template<typename ImageType>
	void 
		calculateHistogram( const ImageType* image, Histogram* histogram )
	{
		histogram->clear();
		unsigned int totalCount = 0; 
		itk::ImageRegionConstIterator<ImageType> it( image, image->GetLargestPossibleRegion() );
		it.GoToBegin();
		while(!it.IsAtEnd()){
			double val = it.Get();
			if (val<histogram->statisticMin || val>histogram->statisticMax){
				++it;
				continue;
			}else{
				unsigned int binIndex = static_cast<unsigned int>(std::floor( (val - histogram->statisticMin)/histogram->binWidth ));
				KM_ASSERT( binIndex >= 0 && binIndex < histogram->binNum );
				histogram->nodes[ binIndex ].count++;
				totalCount++;
				++it;
			}
		}

		double totalPercentage = 0;
		for(int i=0;i<histogram->binNum;i++){
			histogram->nodes[i].lower = histogram->statisticMin+i*histogram->binWidth;
			histogram->nodes[i].upper = histogram->nodes[i].lower + histogram->binWidth - 1;
			histogram->nodes[i].percentage = static_cast<double>(histogram->nodes[i].count)/totalCount;
			totalPercentage += histogram->nodes[i].percentage;
		}
		std::cout<< "Total percentage: " << totalPercentage << std::endl;
	}
}

#endif