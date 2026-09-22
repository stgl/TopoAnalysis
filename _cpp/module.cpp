// module.cpp -- pybind11 bindings for the TopoAnalysis C++ kernels.
//
// Everything exposed here has a pure-NumPy counterpart in
// ``TopoAnalysis.fastops``; the extension is an accelerator, never a hard
// requirement.  Arrays are accepted as C-contiguous float64/uint8/int64.
//
// One rule governs every function below: **resolve every buffer and validate
// every size while the GIL is still held, then release it and run nothing but
// plain C++**.  ``py::array::request()`` calls ``PyObject_GetBuffer``, which
// touches reference counts; calling it inside a ``gil_scoped_release`` block
// corrupts the interpreter as soon as two threads do it at once.  The
// ``Resolved*`` helpers exist so that rule is enforced by construction rather
// than remembered at each call site.

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include "priority_flood.hpp"
#include "flow_routing.hpp"

namespace py = pybind11;
using namespace topoanalysis;

namespace {

using F64Array = py::array_t<double, py::array::c_style | py::array::forcecast>;
using U8Array = py::array_t<std::uint8_t, py::array::c_style | py::array::forcecast>;
using I64Array = py::array_t<std::int64_t, py::array::c_style | py::array::forcecast>;
using I32Array = py::array_t<std::int32_t, py::array::c_style>;

// ---------------------------------------------------------------------------
// GIL-safe buffer resolution
// ---------------------------------------------------------------------------

// Holds an array alive together with the raw pointer into it.  Construct one
// while the GIL is held; use ``.data`` after releasing it.
template <typename Array, typename T>
struct Resolved {
    Array array;
    T* data = nullptr;
    std::size_t size = 0;
    std::vector<py::ssize_t> shape;

    Resolved() = default;

    explicit Resolved(Array a) : array(std::move(a)) {
        auto info = array.request();
        data = static_cast<T*>(info.ptr);
        size = static_cast<std::size_t>(info.size);
        shape.assign(info.shape.begin(), info.shape.end());
    }

    bool present() const { return data != nullptr; }

    void require_2d(const char* name) const {
        if (shape.size() != 2) {
            throw std::invalid_argument(std::string(name) + " must be a 2-D array");
        }
    }

    void require_size(std::size_t expected, const char* name) const {
        if (size != expected) {
            throw std::invalid_argument(
                std::string(name) + " has " + std::to_string(size) +
                " elements; expected " + std::to_string(expected));
        }
    }

    std::size_t ny() const { return static_cast<std::size_t>(shape[0]); }
    std::size_t nx() const { return static_cast<std::size_t>(shape[1]); }
};

using ResolvedF64 = Resolved<F64Array, double>;
using ResolvedConstF64 = Resolved<F64Array, const double>;
using ResolvedU8 = Resolved<U8Array, std::uint8_t>;
using ResolvedI64 = Resolved<I64Array, std::int64_t>;
using ResolvedI32 = Resolved<I32Array, std::int32_t>;

// Optional arguments that may be None.
ResolvedU8 optional_u8(const py::object& obj) {
    if (obj.is_none()) return ResolvedU8();
    return ResolvedU8(obj.cast<U8Array>());
}

ResolvedI64 optional_i64(const py::object& obj) {
    if (obj.is_none()) return ResolvedI64();
    return ResolvedI64(obj.cast<I64Array>());
}

ResolvedF64 optional_f64(const py::object& obj) {
    if (obj.is_none()) return ResolvedF64();
    return ResolvedF64(obj.cast<F64Array>());
}

// An output array of the given 2-D shape, resolved up front.
ResolvedF64 make_f64(std::size_t ny, std::size_t nx) {
    return ResolvedF64(F64Array({static_cast<py::ssize_t>(ny), static_cast<py::ssize_t>(nx)}));
}

ResolvedU8 make_u8(std::size_t ny, std::size_t nx) {
    return ResolvedU8(U8Array({static_cast<py::ssize_t>(ny), static_cast<py::ssize_t>(nx)}));
}

ResolvedI32 make_i32(std::size_t ny, std::size_t nx) {
    return ResolvedI32(I32Array({static_cast<py::ssize_t>(ny), static_cast<py::ssize_t>(nx)}));
}

// ---------------------------------------------------------------------------
// Depression filling
// ---------------------------------------------------------------------------

py::dict py_priority_flood(F64Array elevations,
                           py::object closed_obj,
                           py::object seeds_obj,
                           const std::string& mode,
                           double epsilon,
                           double cellsize,
                           double max_pit_depth,
                           bool track_visited,
                           bool reference_algorithm) {
    ResolvedF64 z(std::move(elevations));
    z.require_2d("elevations");
    const std::size_t ny = z.ny(), nx = z.nx(), n = ny * nx;

    std::vector<std::uint8_t> closed(n, 0);
    ResolvedU8 closed_arg = optional_u8(closed_obj);
    if (closed_arg.present()) {
        closed_arg.require_size(n, "closed");
        closed.assign(closed_arg.data, closed_arg.data + n);
    }
    for (std::size_t k = 0; k < n; ++k) {
        if (std::isnan(z.data[k])) closed[k] = 1;
    }

    std::vector<std::int64_t> seeds;
    ResolvedI64 seeds_arg = optional_i64(seeds_obj);
    if (!seeds_arg.present()) {
        seeds = default_seeds(z.data, ny, nx);
        // Honour a caller-supplied mask: a seed that is masked out would
        // otherwise silently disappear and leave the grid unflooded.
        std::vector<std::int64_t> kept;
        kept.reserve(seeds.size());
        for (std::int64_t k : seeds) {
            if (!closed[static_cast<std::size_t>(k)]) kept.push_back(k);
        }
        if (closed_arg.present() && kept.empty()) {
            // The mask excluded the whole grid perimeter.  Fall back to the
            // boundary of the masked region itself.
            for (std::size_t i = 0; i < ny; ++i) {
                for (std::size_t j = 0; j < nx; ++j) {
                    const std::size_t k = i * nx + j;
                    if (closed[k]) continue;
                    bool boundary = (i == 0 || j == 0 || i + 1 == ny || j + 1 == nx);
                    for (int d = 0; d < 8 && !boundary; ++d) {
                        const std::int64_t ni = static_cast<std::int64_t>(i) + kD8Di[d];
                        const std::int64_t nj = static_cast<std::int64_t>(j) + kD8Dj[d];
                        if (ni < 0 || nj < 0 || static_cast<std::size_t>(ni) >= ny ||
                            static_cast<std::size_t>(nj) >= nx) {
                            continue;
                        }
                        if (closed[static_cast<std::size_t>(ni) * nx +
                                   static_cast<std::size_t>(nj)]) {
                            boundary = true;
                        }
                    }
                    if (boundary) kept.push_back(static_cast<std::int64_t>(k));
                }
            }
        }
        seeds.swap(kept);
    } else {
        seeds.assign(seeds_arg.data, seeds_arg.data + seeds_arg.size);
    }

    FloodOptions options;
    if (mode == "epsilon") {
        options.mode = FillMode::Epsilon;
    } else if (mode == "flat") {
        options.mode = FillMode::Flat;
    } else {
        throw std::invalid_argument("mode must be 'flat' or 'epsilon'");
    }
    options.epsilon = epsilon;
    options.cellsize = cellsize;
    options.max_pit_depth = max_pit_depth;
    options.track_visited = track_visited;
    options.track_unfilled = max_pit_depth > 0.0;

    FloodResult result;
    {
        py::gil_scoped_release release;
        result = reference_algorithm
                     ? priority_flood_original(z.data, ny, nx, closed, seeds, options)
                     : priority_flood(z.data, ny, nx, closed, seeds, options);
    }

    py::dict out;
    out["cells_filled"] = result.cells_filled;
    out["max_fill_depth"] = result.max_fill_depth;
    if (track_visited) {
        ResolvedU8 visited = make_u8(ny, nx);
        std::copy(result.visited.begin(), result.visited.end(), visited.data);
        out["visited"] = visited.array;
    } else {
        out["visited"] = py::none();
    }
    if (options.track_unfilled) {
        ResolvedU8 unfilled = make_u8(ny, nx);
        std::copy(result.unfilled.begin(), result.unfilled.end(), unfilled.data);
        out["unfilled"] = unfilled.array;
    } else {
        out["unfilled"] = py::none();
    }
    return out;
}

// ---------------------------------------------------------------------------
// Flow routing
// ---------------------------------------------------------------------------

U8Array py_flow_directions(F64Array elevations, py::object cellsize_obj,
                           double uniform_cellsize) {
    ResolvedF64 z(std::move(elevations));
    z.require_2d("elevations");
    const std::size_t ny = z.ny(), nx = z.nx();

    ResolvedF64 cellsize = optional_f64(cellsize_obj);
    if (cellsize.present()) cellsize.require_size(ny * nx, "cellsize_grid");

    ResolvedU8 out = make_u8(ny, nx);
    {
        py::gil_scoped_release release;
        d8_flow_directions(z.data, ny, nx, cellsize.data, uniform_cellsize, out.data);
    }
    return out.array;
}

I64Array py_receivers(U8Array codes) {
    ResolvedU8 c(std::move(codes));
    c.require_2d("codes");
    const std::size_t ny = c.ny(), nx = c.nx();

    ResolvedI64 out(I64Array({static_cast<py::ssize_t>(ny), static_cast<py::ssize_t>(nx)}));
    {
        py::gil_scoped_release release;
        d8_receivers(c.data, ny, nx, out.data);
    }
    return out.array;
}

py::tuple py_topological_order(I64Array receivers) {
    ResolvedI64 recv(std::move(receivers));
    const std::size_t n = recv.size;

    ResolvedI64 out(I64Array(static_cast<py::ssize_t>(n)));
    std::size_t cycles = 0;
    {
        py::gil_scoped_release release;
        cycles = d8_topological_order(recv.data, n, out.data);
    }
    return py::make_tuple(out.array, cycles);
}

F64Array py_accumulate(I64Array receivers, I64Array order, F64Array weights,
                       py::object gate_obj) {
    ResolvedI64 recv(std::move(receivers));
    ResolvedI64 ord(std::move(order));
    ResolvedF64 w(std::move(weights));
    const std::size_t n = recv.size;
    w.require_size(n, "weights");
    ord.require_size(n, "order");

    ResolvedU8 gate = optional_u8(gate_obj);
    if (gate.present()) gate.require_size(n, "gate");

    ResolvedF64 out(F64Array(w.shape));
    {
        py::gil_scoped_release release;
        d8_accumulate(recv.data, ord.data, w.data, gate.data, n, out.data);
    }
    return out.array;
}

py::tuple py_flow_length(I64Array receivers, I64Array order, F64Array step_length) {
    ResolvedI64 recv(std::move(receivers));
    recv.require_2d("receivers");
    const std::size_t ny = recv.ny(), nx = recv.nx(), n = ny * nx;

    ResolvedI64 ord(std::move(order));
    ord.require_size(n, "order");
    ResolvedF64 step(std::move(step_length));
    step.require_size(n, "step_length");

    ResolvedF64 length = make_f64(ny, nx);
    ResolvedU8 from_codes = make_u8(ny, nx);
    {
        py::gil_scoped_release release;
        d8_flow_length(recv.data, ord.data, step.data, ny, nx, length.data, from_codes.data);
    }
    return py::make_tuple(length.array, from_codes.array);
}

F64Array py_propagate_along_main_stem(I64Array receivers, I64Array order,
                                      U8Array main_stem_from_codes, py::object gate_obj,
                                      F64Array values, const std::string& mode) {
    if (mode != "sum" && mode != "carry") {
        throw std::invalid_argument("mode must be 'carry' or 'sum'");
    }

    ResolvedI64 recv(std::move(receivers));
    recv.require_2d("receivers");
    const std::size_t ny = recv.ny(), nx = recv.nx(), n = ny * nx;

    ResolvedI64 ord(std::move(order));
    ord.require_size(n, "order");
    ResolvedU8 codes(std::move(main_stem_from_codes));
    codes.require_size(n, "main_stem_from_codes");
    ResolvedF64 vals(std::move(values));
    vals.require_size(n, "values");
    ResolvedU8 gate = optional_u8(gate_obj);
    if (gate.present()) gate.require_size(n, "gate");

    ResolvedF64 out = make_f64(ny, nx);
    std::copy(vals.data, vals.data + n, out.data);
    {
        py::gil_scoped_release release;
        d8_propagate_along_main_stem(recv.data, ord.data, codes.data, gate.data, ny, nx,
                                     mode == "sum" ? PropagateMode::Sum : PropagateMode::Carry,
                                     out.data);
    }
    return out.array;
}

U8Array py_upstream_mask(I64Array receivers, I64Array outlets) {
    ResolvedI64 recv(std::move(receivers));
    recv.require_2d("receivers");
    const std::size_t ny = recv.ny(), nx = recv.nx();
    ResolvedI64 out_idx(std::move(outlets));

    ResolvedU8 out = make_u8(ny, nx);
    {
        py::gil_scoped_release release;
        d8_upstream_mask(recv.data, ny * nx, out_idx.data, out_idx.size, out.data);
    }
    return out.array;
}

py::tuple py_chi(I64Array receivers, I64Array order, F64Array area, F64Array step_length,
                 I64Array outlets, double A0, double theta, bool trapezoid, double max_length,
                 py::object mask_obj) {
    ResolvedI64 recv(std::move(receivers));
    recv.require_2d("receivers");
    const std::size_t ny = recv.ny(), nx = recv.nx(), n = ny * nx;

    ResolvedI64 ord(std::move(order));
    ord.require_size(n, "order");
    ResolvedF64 a(std::move(area));
    a.require_size(n, "area");
    ResolvedF64 step(std::move(step_length));
    step.require_size(n, "step_length");
    ResolvedI64 out_idx(std::move(outlets));
    ResolvedU8 mask = optional_u8(mask_obj);
    if (mask.present()) mask.require_size(n, "mask");

    ResolvedF64 chi_out = make_f64(ny, nx);
    ResolvedF64 distance = make_f64(ny, nx);
    {
        py::gil_scoped_release release;
        d8_chi(recv.data, ord.data, a.data, step.data, n, out_idx.data, out_idx.size,
               A0, theta, trapezoid, max_length, mask.data, chi_out.data, distance.data);
    }
    return py::make_tuple(chi_out.array, distance.array);
}

F64Array py_downstream_distance(I64Array receivers, I64Array order, F64Array step_length) {
    ResolvedI64 recv(std::move(receivers));
    recv.require_2d("receivers");
    const std::size_t ny = recv.ny(), nx = recv.nx(), n = ny * nx;

    ResolvedI64 ord(std::move(order));
    ord.require_size(n, "order");
    ResolvedF64 step(std::move(step_length));
    step.require_size(n, "step_length");

    ResolvedF64 out = make_f64(ny, nx);
    {
        py::gil_scoped_release release;
        d8_downstream_distance(recv.data, ord.data, step.data, n, out.data);
    }
    return out.array;
}

void py_imposemin(I64Array receivers, I64Array order, F64Array step_length, double sl,
                  py::array_t<double, py::array::c_style> elevations) {
    ResolvedI64 recv(std::move(receivers));
    const std::size_t n = recv.size;
    ResolvedI64 ord(std::move(order));
    ord.require_size(n, "order");
    ResolvedF64 step(std::move(step_length));
    step.require_size(n, "step_length");

    auto einfo = elevations.request();
    if (static_cast<std::size_t>(einfo.size) != n) {
        throw std::invalid_argument("elevations must match the receiver grid");
    }
    double* z = static_cast<double*>(einfo.ptr);

    py::gil_scoped_release release;
    d8_imposemin(recv.data, ord.data, step.data, n, sl, z);
}

I32Array py_stream_order(I64Array receivers, I64Array order, U8Array is_stream,
                         const std::string& kind) {
    if (kind != "strahler" && kind != "shreve") {
        throw std::invalid_argument("kind must be 'strahler' or 'shreve'");
    }

    ResolvedI64 recv(std::move(receivers));
    recv.require_2d("receivers");
    const std::size_t ny = recv.ny(), nx = recv.nx(), n = ny * nx;

    ResolvedI64 ord(std::move(order));
    ord.require_size(n, "order");
    ResolvedU8 streams(std::move(is_stream));
    streams.require_size(n, "is_stream");

    ResolvedI32 out = make_i32(ny, nx);
    {
        py::gil_scoped_release release;
        if (kind == "shreve") {
            d8_shreve(recv.data, ord.data, streams.data, n, out.data);
        } else {
            d8_strahler(recv.data, ord.data, streams.data, n, out.data);
        }
    }
    return out.array;
}

I32Array py_drainage_basins(I64Array receivers, I64Array order, py::object outlets_obj,
                            py::object valid_obj) {
    ResolvedI64 recv(std::move(receivers));
    recv.require_2d("receivers");
    const std::size_t ny = recv.ny(), nx = recv.nx(), n = ny * nx;

    ResolvedI64 ord(std::move(order));
    ord.require_size(n, "order");
    ResolvedU8 valid = optional_u8(valid_obj);
    if (valid.present()) valid.require_size(n, "valid");
    ResolvedI64 outlets = optional_i64(outlets_obj);

    ResolvedI32 out = make_i32(ny, nx);
    {
        py::gil_scoped_release release;
        if (!outlets.present()) {
            d8_drainage_basins_all(recv.data, ord.data, valid.data, n, out.data);
        } else {
            // `valid` applies here too; dropping it made the two backends
            // disagree whenever both arguments were supplied.
            d8_drainage_basins(recv.data, ord.data, n, outlets.data, outlets.size,
                               valid.data, out.data);
        }
    }
    return out.array;
}

}  // namespace

PYBIND11_MODULE(_topoanalysis, m) {
    m.doc() =
        "Compiled kernels for TopoAnalysis.\n\n"
        "Depression filling implements the Priority-Flood family of algorithms of\n"
        "Barnes, R., Lehman, C. and Mulla, D. (2014), 'Priority-flood: An optimal\n"
        "depression-filling and watershed-labeling algorithm for digital elevation\n"
        "models', Computers & Geosciences 62, 117-127.";

    m.attr("__algorithm_reference__") =
        "Barnes, R., Lehman, C., Mulla, D. (2014). Priority-flood: An optimal "
        "depression-filling and watershed-labeling algorithm for digital elevation "
        "models. Computers & Geosciences 62, 117-127. doi:10.1016/j.cageo.2013.04.024";

    m.def("priority_flood", &py_priority_flood, py::arg("elevations").noconvert(),
          py::arg("closed") = py::none(), py::arg("seeds") = py::none(),
          py::arg("mode") = "flat", py::arg("epsilon") = 0.0, py::arg("cellsize") = 1.0,
          py::arg("max_pit_depth") = 0.0, py::arg("track_visited") = false,
          py::arg("reference_algorithm") = false,
          "Fill depressions in place using Priority-Flood (Barnes et al. 2014).");

    m.def("flow_directions", &py_flow_directions, py::arg("elevations"),
          py::arg("cellsize_grid") = py::none(), py::arg("cellsize") = 1.0,
          "Steepest-descent D8 flow directions as ArcGIS codes.");

    m.def("receivers", &py_receivers, py::arg("codes"),
          "Flat receiver index for each cell (-1 where flow leaves the grid).");

    m.def("topological_order", &py_topological_order, py::arg("receivers"),
          "Donors-before-receivers ordering; returns (order, n_cells_on_cycles).");

    m.def("accumulate", &py_accumulate, py::arg("receivers"), py::arg("order"),
          py::arg("weights"), py::arg("gate") = py::none(),
          "Accumulate weights downstream along the D8 network.");

    m.def("flow_length", &py_flow_length, py::arg("receivers"), py::arg("order"),
          py::arg("step_length"),
          "Longest upstream flow distance and the main-stem donor direction codes.");

    m.def("propagate_along_main_stem", &py_propagate_along_main_stem, py::arg("receivers"),
          py::arg("order"), py::arg("main_stem_from_codes"), py::arg("gate"), py::arg("values"),
          py::arg("mode"), "Carry or sum values downstream along the longest flow path.");

    m.def("upstream_mask", &py_upstream_mask, py::arg("receivers"), py::arg("outlets"),
          "Mask of every cell draining to any of the given outlets.");

    m.def("chi", &py_chi, py::arg("receivers"), py::arg("order"), py::arg("area"),
          py::arg("step_length"), py::arg("outlets"), py::arg("A0"), py::arg("theta"),
          py::arg("trapezoid") = false, py::arg("max_length") = 0.0,
          py::arg("mask") = py::none(), "Integrate chi upstream from the given outlets.");

    m.def("drainage_basins", &py_drainage_basins, py::arg("receivers"), py::arg("order"),
          py::arg("outlets") = py::none(), py::arg("valid") = py::none(),
          "Label drainage basins, either for given outlets or for every terminal cell.");

    m.def("downstream_distance", &py_downstream_distance, py::arg("receivers"), py::arg("order"),
          py::arg("step_length"),
          "Along-flow distance from each cell to the end of its flow path.");

    m.def("imposemin", &py_imposemin, py::arg("receivers"), py::arg("order"),
          py::arg("step_length"), py::arg("sl"), py::arg("elevations").noconvert(),
          "Carve elevations in place so flow paths descend at gradient >= sl.");

    m.def("stream_order", &py_stream_order, py::arg("receivers"), py::arg("order"),
          py::arg("is_stream"), py::arg("kind") = "strahler",
          "Strahler or Shreve stream order over the marked stream cells.");
}
