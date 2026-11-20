#include <doctest.h>

#include "multiax/core/build_info.hpp"
#include "multiax/core/types.hpp"

#include <type_traits>

TEST_CASE("build info: C++20 and OpenMP")
{
	const auto info = MultiaxCore::buildInfo();
	CHECK(info.cplusplus >= 202002L);
	CHECK(info.openmp_version > 0);
	CHECK(info.max_threads >= 1);
	CHECK_FALSE(info.version.empty());
}

TEST_CASE("core types")
{
	static_assert(std::is_floating_point_v<MultiaxCore::MXScalar>);
	CHECK(MultiaxCore::MXMachineEpsilon > 0);
}
