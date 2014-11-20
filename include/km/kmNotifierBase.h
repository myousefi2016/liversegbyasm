#ifndef KM_NOTIFIER_BASE_H
#define KM_NOTIFIER_BASE_H

namespace km
{
	template<class TImage, class TPointSet>
	class NotifierBase
	{
	public:
		virtual void notifyMesh( const TPointSet* pointset, const char* filename = NULL ) const
		{
		}
		
		virtual void notifyImage( const TImage* image, const char* filename = NULL ) const
		{
		}
	};
}

#endif