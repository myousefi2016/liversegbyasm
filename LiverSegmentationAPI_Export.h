
#ifndef LiverSegmentationAPI_EXPORT_H
#define LiverSegmentationAPI_EXPORT_H

#ifdef LiverSegmentationAPI_BUILT_AS_STATIC
#  define LiverSegmentationAPI_EXPORT __declspec(dllexport)
#  define LIVERSEGMENTATIONAPI_NO_EXPORT
#else
#  ifndef LiverSegmentationAPI_EXPORT
#    ifdef LiverSegmentationAPI_EXPORTS
        /* We are building this library */
#      define LiverSegmentationAPI_EXPORT  __declspec(dllexport)
#    else
        /* We are using this library */
#      define LiverSegmentationAPI_EXPORT 
#    endif
#  endif

#  ifndef LIVERSEGMENTATIONAPI_NO_EXPORT
#    define LIVERSEGMENTATIONAPI_NO_EXPORT 
#  endif
#endif

#ifndef LIVERSEGMENTATIONAPI_DEPRECATED
#  define LIVERSEGMENTATIONAPI_DEPRECATED __declspec(deprecated)
#  define LIVERSEGMENTATIONAPI_DEPRECATED_EXPORT LiverSegmentationAPI_EXPORT __declspec(deprecated)
#  define LIVERSEGMENTATIONAPI_DEPRECATED_NO_EXPORT LIVERSEGMENTATIONAPI_NO_EXPORT __declspec(deprecated)
#endif

#define DEFINE_NO_DEPRECATED 0
#if DEFINE_NO_DEPRECATED
# define LIVERSEGMENTATIONAPI_NO_DEPRECATED
#endif

#endif
