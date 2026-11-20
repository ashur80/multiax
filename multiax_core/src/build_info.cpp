#include "multiax/core/build_info.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef MULTIAX_VERSION
#define MULTIAX_VERSION "0.0.0"
#endif

namespace MultiaxCore
{
	BuildInfo buildInfo()
	{
		BuildInfo info;
		info.version = MULTIAX_VERSION;
#if defined(__clang__)
		info.compiler = std::string("clang ") + __clang_version__;
#elif defined(_MSC_VER)
		info.compiler = "msvc " + std::to_string(_MSC_VER);
#elif defined(__GNUC__)
		info.compiler = std::string("gcc ") + __VERSION__;
#else
		info.compiler = "unknown";
#endif
		info.cplusplus = __cplusplus;
#ifdef _OPENMP
		info.openmp_version = _OPENMP;
		info.max_threads = omp_get_max_threads();
#else
		info.openmp_version = 0;
		info.max_threads = 1;
#endif
		return info;
	}
}
