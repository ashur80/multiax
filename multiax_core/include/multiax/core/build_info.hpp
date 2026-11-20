#ifndef __MXBuildInfo_hpp__
#define __MXBuildInfo_hpp__

#include <string>

namespace MultiaxCore
{
	// Version and build configuration of the core library (compiler, C++ standard, OpenMP).
	struct BuildInfo
	{
		std::string version;
		std::string compiler;
		long cplusplus;
		int openmp_version; // _OPENMP (yyyymm), 0 if built without OpenMP
		int max_threads;    // omp_get_max_threads(), 1 without OpenMP
	};

	BuildInfo buildInfo();
}

#endif //__MXBuildInfo_hpp__
