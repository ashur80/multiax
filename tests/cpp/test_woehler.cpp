#include <doctest.h>

#include "multiax/woehler/woehler_curve.hpp"

using MultiaxCore::WoehlerCurve;

TEST_CASE("WoehlerCurve: parameters are stored")
{
	const WoehlerCurve wc(3.0, 22.0, 109.19, 5.0e5);
	CHECK(wc.k() == doctest::Approx(3.0));
	CHECK(wc.kStar() == doctest::Approx(22.0));
	CHECK(wc.lK() == doctest::Approx(109.19));
	CHECK(wc.NK() == doctest::Approx(5.0e5));
}

TEST_CASE("WoehlerCurve: default knee point N_K = 1e6 if not given")
{
	const WoehlerCurve wc(3.0, 22.0, 100.0, 0.0);
	CHECK(wc.NK() == doctest::Approx(1.0e6));
}

// Exercise E1 (docs/learning/E1_woehler_soa.md): LofN / NofL, Eq. (2.1).
// These tests are skipped until E1 is started; remove the skip decorator to make them active (red first).
TEST_CASE("E1: WoehlerCurve LofN / NofL, Eq. (2.1)" * doctest::skip(true))
{
	const WoehlerCurve wc(3.0, 22.0, 109.19, 5.0e5);
	CHECK(wc.NofL(109.19) == doctest::Approx(5.0e5));
	CHECK(wc.LofN(5.0e5) == doctest::Approx(109.19));
	CHECK(wc.NofL(2.0 * 109.19) == doctest::Approx(5.0e5 / 8.0));
	CHECK(wc.LofN(wc.NofL(80.0)) == doctest::Approx(80.0)); // below the knee: slope kStar
}
