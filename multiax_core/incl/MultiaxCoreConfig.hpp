#ifndef __MultiaxCoreConfig_hpp__
#define __MultiaxCoreConfig_hpp__

#ifdef MULTIAXCORE_DLL
	#define MULTIAXCORE_DLLEXPORTS __declspec( dllexport )
#else
	#define MULTIAXCORE_DLLEXPORTS __declspec( dllimport )
#endif

#endif // __MultiaxCoreConfig_hpp__