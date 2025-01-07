#ifndef __MXTypes_hpp__
#define __MXTypes_hpp__

#include <limits>

namespace MultiaxCore
{

#define MX_USES_DOUBLES

#ifdef MX_USES_DOUBLES
	using MXScalar = double; 
#else
	using MXScalar = float;
#endif

	static constexpr MXScalar MXMachineEpsilon = std::numeric_limits<MXScalar>::epsilon();

}
#endif //__MXTypes_hpp__