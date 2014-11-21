#ifndef KM_NOTIFIER_BASE_H
#define KM_NOTIFIER_BASE_H

namespace km
{
	template<class TPointSet>
	class NotifierBase
	{
	public:
		virtual void notify( const TPointSet * pointset )
		{
		}
	};
}

#endif