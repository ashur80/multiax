// Benchmarks (nanobench). Run the Release build: build/release/bench/multiax_bench.exe
// Exercises add their benchmarks here (E1: AoS vs SoA, E5/E6: rainflow serial vs parallel).

#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "multiax/core/build_info.hpp"

#include <cmath>
#include <cstdio>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

int main()
{
	const auto info = MultiaxCore::buildInfo();
	std::printf("multiax %s | %s | C++ %ld | OpenMP %d, %d threads\n\n",
	            info.version.c_str(), info.compiler.c_str(), info.cplusplus, info.openmp_version, info.max_threads);

	// Baseline: serial vs OpenMP transform over 1e7 doubles (memory-bound reference point)
	const std::size_t n = 10'000'000;
	std::vector<double> x(n, 1.5), y(n);

	ankerl::nanobench::Bench b;
	b.title("baseline: y = sqrt(x) * 2").unit("element").batch(n).minEpochIterations(3);
	b.run("serial", [&] {
		for (std::size_t i = 0; i < n; ++i)
			y[i] = std::sqrt(x[i]) * 2.0;
		ankerl::nanobench::doNotOptimizeAway(y.data());
	});
	b.run("openmp", [&] {
		const auto m = static_cast<std::ptrdiff_t>(n);
#pragma omp parallel for schedule(static)
		for (std::ptrdiff_t i = 0; i < m; ++i)
			y[static_cast<std::size_t>(i)] = std::sqrt(x[static_cast<std::size_t>(i)]) * 2.0;
		ankerl::nanobench::doNotOptimizeAway(y.data());
	});
	return 0;
}
