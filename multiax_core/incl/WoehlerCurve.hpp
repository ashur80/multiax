#ifndef __WoehlerCurve_hpp__
#define __WoehlerCurve_hpp__

#include  "MultiaxCoreConfig.hpp"
#include "mxtypes.hpp"

#include <concepts>

namespace MultiaxCore
{
	template <std::floating_point TNum>
	class WoehlerCurveT
	{
	public:
		WoehlerCurveT(TNum i_k, TNum i_kStar, TNum i_lk, TNum i_NK);

		//compute load, if N is given
		constexpr TNum LofN(TNum i_N) const;
		
		//compute number of cycles, if load is given
		constexpr TNum NofL(TNum i_L) const;

	private:
		TNum k, kStar; //inclination before and after the knee point
		TNum lK, NK; //load and number of cycles at the knee point

	};

	using WoehlerCurve = WoehlerCurveT <MXScalar> ;
}
	//Template implementations

template<std::floating_point TNum>
MultiaxCore::WoehlerCurveT<TNum>::WoehlerCurveT(TNum i_k, TNum i_kStar, TNum i_lK, TNum i_NK) :
	k(i_k), kStar(i_kStar), lK(i_lK), NK(i_NK)
{
	if (NK <= MXMachineEpsilon)
	{
		NK = 1000000.; //10^6
	}
}

template<std::floating_point TNum>
inline constexpr TNum MultiaxCore::WoehlerCurveT<TNum>::LofN(TNum i_N) const
{
	return TNum();
}

template<std::floating_point TNum>
inline constexpr TNum MultiaxCore::WoehlerCurveT<TNum>::NofL(TNum i_L) const
{
	return TNum();
}


#endif //__WoehlerCurve_hpp__

