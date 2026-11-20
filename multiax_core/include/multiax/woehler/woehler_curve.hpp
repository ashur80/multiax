#ifndef __WoehlerCurve_hpp__
#define __WoehlerCurve_hpp__

#include "multiax/core/types.hpp"

#include <concepts>

namespace MultiaxCore
{
	// Bilinear S-N curve, Eq. (2.1): slope k above and kStar below the knee point (lK, NK).
	template <std::floating_point TNum>
	class WoehlerCurveT
	{
	public:
		WoehlerCurveT(TNum i_k, TNum i_kStar, TNum i_lK, TNum i_NK);

		//compute load, if N is given
		constexpr TNum LofN(TNum i_N) const;

		//compute number of cycles, if load is given
		constexpr TNum NofL(TNum i_L) const;

		constexpr TNum k() const { return m_k; }
		constexpr TNum kStar() const { return m_kStar; }
		constexpr TNum lK() const { return m_lK; }
		constexpr TNum NK() const { return m_NK; }

	private:
		TNum m_k, m_kStar; //inclination before and after the knee point
		TNum m_lK, m_NK; //load and number of cycles at the knee point

	};

	using WoehlerCurve = WoehlerCurveT <MXScalar> ;
}
	//Template implementations

template<std::floating_point TNum>
MultiaxCore::WoehlerCurveT<TNum>::WoehlerCurveT(TNum i_k, TNum i_kStar, TNum i_lK, TNum i_NK) :
	m_k(i_k), m_kStar(i_kStar), m_lK(i_lK), m_NK(i_NK)
{
	if (m_NK <= MXMachineEpsilon)
	{
		m_NK = 1000000.; //10^6
	}
}

template<std::floating_point TNum>
inline constexpr TNum MultiaxCore::WoehlerCurveT<TNum>::LofN(TNum i_N) const
{
	// TODO(E1): Eq. (2.1), L(N) = lK (NK/N)^(1/k) for N <= NK, else lK (NK/N)^(1/kStar)
	(void)i_N;
	return TNum();
}

template<std::floating_point TNum>
inline constexpr TNum MultiaxCore::WoehlerCurveT<TNum>::NofL(TNum i_L) const
{
	// TODO(E1): Eq. (2.1), N(L) = NK (lK/L)^k for L >= lK, else NK (lK/L)^kStar
	(void)i_L;
	return TNum();
}


#endif //__WoehlerCurve_hpp__
