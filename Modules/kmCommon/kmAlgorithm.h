#ifndef KM_ALGORITHM_H
#define KM_ALGORITHM_H

namespace km
{	
	//VectorType has to be fixed length vector.
	template<class VectorType>
	void sort_insert(VectorType& vec, bool asc = true)
	{
		for (int i=1;i<VectorType::Dimension;i++){
			int k = i;
			if(!asc){
				while( k>0 && vec[k]>vec[k-1] ){
					std::swap(vec[k], vec[k-1]);
					k--;
				}
			}else{
				while( k>0 && vec[k]<vec[k-1] ){
					std::swap(vec[k], vec[k-1]);
					k--;
				}
			}
		}
	}
	
	template<class T>
	void sort_insert(T* array1, int length, bool asc = true)
	{
		for (int i=1;i<length;i++){
			int k = i;
			if(!asc){
				while( k>0 && array1[k]>array1[k-1] ){
					std::swap(array1[k], array1[k-1]);
					k--;
				}
			}else{
				while( k>0 && array1[k]<array1[k-1] ){
					std::swap(array1[k], array1[k-1]);
					k--;
				}
			}
		}
	}

	template<class T1, class T2>
	void sort_insert(T1* array1, T2* array2, int length, bool asc = true)
	{
		for (int i=1;i<length;i++){
			int k = i;
			if(!asc){
				while( k>0 && array1[k]>array1[k-1] ){
					std::swap(array1[k], array1[k-1]);
					std::swap(array2[k], array2[k-1]);
					k--;
				}
			}else{
				while( k>0 && array1[k]<array1[k-1] ){
					std::swap(array1[k], array1[k-1]);
					std::swap(array2[k], array2[k-1]);
					k--;
				}
			}
		}
	}
}

#endif