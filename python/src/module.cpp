// Python bindings of the multiax core (module multiax._core).
// Convention: numeric inputs are NumPy arrays passed without copy (buffer protocol -> std::span).

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "multiax/core/build_info.hpp"
#include "multiax/woehler/woehler_curve.hpp"

namespace py = pybind11;
using namespace MultiaxCore;

PYBIND11_MODULE(_core, m)
{
	m.doc() = "multiax C++ core (Multiax-Invariant, Bolchoun 2018)";

	m.def("build_info", [] {
		const auto bi = buildInfo();
		py::dict d;
		d["version"] = bi.version;
		d["compiler"] = bi.compiler;
		d["cplusplus"] = bi.cplusplus;
		d["openmp_version"] = bi.openmp_version;
		d["max_threads"] = bi.max_threads;
		return d;
	}, "Version and build configuration (compiler, C++ standard, OpenMP).");

	py::class_<WoehlerCurve>(m, "WoehlerCurve", "Bilinear S-N curve, Eq. (2.1).")
		.def(py::init<MXScalar, MXScalar, MXScalar, MXScalar>(), py::arg("k"), py::arg("k_star"), py::arg("L_k"), py::arg("N_k"))
		.def_property_readonly("k", &WoehlerCurve::k)
		.def_property_readonly("k_star", &WoehlerCurve::kStar)
		.def_property_readonly("L_k", &WoehlerCurve::lK)
		.def_property_readonly("N_k", &WoehlerCurve::NK)
		.def("L_of_N", &WoehlerCurve::LofN, py::arg("N"))
		.def("N_of_L", &WoehlerCurve::NofL, py::arg("L"))
		.def("__repr__", [](const WoehlerCurve& w) {
			return "WoehlerCurve(k=" + std::to_string(w.k()) + ", k_star=" + std::to_string(w.kStar()) +
			       ", L_k=" + std::to_string(w.lK()) + ", N_k=" + std::to_string(w.NK()) + ")";
		});
}
