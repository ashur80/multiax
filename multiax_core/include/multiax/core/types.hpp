#ifndef __MXTypes_hpp__
#define __MXTypes_hpp__

// Basic numeric types and portability macros of the multiax core.

#include <cstddef>
#include <limits>

// Marks functions that must later also compile as CUDA device code.
// Kernels marked with MX_HOST_DEVICE: no std containers, no exceptions, no allocation.
#if defined(__CUDACC__)
#define MX_HOST_DEVICE __host__ __device__
#else
#define MX_HOST_DEVICE
#endif

namespace MultiaxCore
{

#define MX_USES_DOUBLES

#ifdef MX_USES_DOUBLES
	using MXScalar = double;
#else
	using MXScalar = float;
#endif

	using MXIndex = std::ptrdiff_t;

	static constexpr MXScalar MXMachineEpsilon = std::numeric_limits<MXScalar>::epsilon();

}
#endif //__MXTypes_hpp__
