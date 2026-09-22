// flow_routing.hpp -- D8 flow routing kernels
//
// These replace the per-cell Python loops and the deeply recursive traversals
// that used to dominate TopoAnalysis' runtime.  Everything here is iterative,
// so there is no recursion limit to raise and no stack to overflow on large
// basins.
//
// Flow-direction codes follow the ArcGIS convention used throughout
// TopoAnalysis, with row 0 at the north edge of the grid:
//
//     | 32  64 128 |      | (i-1,j-1) (i-1,j) (i-1,j+1) |
//     | 16   X   1 |  ==  | (i  ,j-1)    X    (i  ,j+1) |
//     |  8   4   2 |      | (i+1,j-1) (i+1,j) (i+1,j+1) |
//
// A code of 0 means "no downstream neighbour" (a pit, an outlet, or no-data).

#ifndef TOPOANALYSIS_FLOW_ROUTING_HPP
#define TOPOANALYSIS_FLOW_ROUTING_HPP

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstddef>
#include <limits>
#include <vector>

#include "priority_flood.hpp"

namespace topoanalysis {

// ---------------------------------------------------------------------------
// Flow direction
// ---------------------------------------------------------------------------

// Steepest-descent D8 flow directions.
//
// The drop to each neighbour is divided by the centre-to-centre distance, so
// diagonal neighbours compete on gradient rather than on raw elevation drop.
// This is the standard D8 rule (O'Callaghan & Mark 1984) and matches
// TopoToolbox's ``gradient8``/FLOWobj steepest-neighbour selection.
//
// ``cellsize`` may be a single value (uniform grid) or a per-cell array of
// mean pixel dimensions, which is what geographic (lat/lon) grids need.
// Ties are broken by neighbour order (E, SE, S, SW, W, NW, N, NE), which is
// deterministic and reproducible.
inline void d8_flow_directions(const double* elevations,
                               std::size_t ny,
                               std::size_t nx,
                               const double* cellsize,  // ny*nx, or nullptr
                               double uniform_cellsize,
                               std::uint8_t* out_codes) {
    const std::int64_t rows = static_cast<std::int64_t>(ny);
    const std::int64_t cols = static_cast<std::int64_t>(nx);

    for (std::int64_t i = 0; i < rows; ++i) {
        for (std::int64_t j = 0; j < cols; ++j) {
            const std::int64_t k = i * cols + j;
            const double z = elevations[k];
            if (std::isnan(z)) {
                out_codes[k] = 0;
                continue;
            }
            const double de = cellsize ? cellsize[k] : uniform_cellsize;

            double best_slope = 0.0;
            int best_dir = -1;
            for (int d = 0; d < 8; ++d) {
                const std::int64_t ni = i + kD8Di[d];
                const std::int64_t nj = j + kD8Dj[d];
                if (ni < 0 || nj < 0 || ni >= rows || nj >= cols) continue;
                const double zn = elevations[ni * cols + nj];
                if (std::isnan(zn)) continue;
                const double slope = (z - zn) / (de * kD8Dist[d]);
                if (slope > best_slope) {
                    best_slope = slope;
                    best_dir = d;
                }
            }
            out_codes[k] = (best_dir < 0) ? 0 : static_cast<std::uint8_t>(1u << best_dir);
        }
    }
}

// Translate a flow-direction code grid into an explicit receiver index per
// cell (-1 where flow leaves the grid or the cell has no receiver).
inline void d8_receivers(const std::uint8_t* codes,
                         std::size_t ny,
                         std::size_t nx,
                         std::int64_t* out_receivers) {
    const std::int64_t rows = static_cast<std::int64_t>(ny);
    const std::int64_t cols = static_cast<std::int64_t>(nx);

    for (std::int64_t i = 0; i < rows; ++i) {
        for (std::int64_t j = 0; j < cols; ++j) {
            const std::int64_t k = i * cols + j;
            const std::uint8_t code = codes[k];
            out_receivers[k] = -1;
            if (code == 0) continue;
            for (int d = 0; d < 8; ++d) {
                if (code != static_cast<std::uint8_t>(1u << d)) continue;
                const std::int64_t ni = i + kD8Di[d];
                const std::int64_t nj = j + kD8Dj[d];
                if (ni < 0 || nj < 0 || ni >= rows || nj >= cols) break;
                out_receivers[k] = ni * cols + nj;
                break;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Topological ordering
// ---------------------------------------------------------------------------

// Order the cells so that every cell appears before the cell it drains into
// (Kahn's algorithm on the donor counts).  This is the ordering flow
// accumulation needs, and unlike sorting by elevation it is derived from the
// flow network itself -- so it stays correct for flow-direction grids that
// were loaded from disk, hand-edited, or derived from a different DEM.
//
// Returns the number of cells that could NOT be ordered because they lie on a
// cycle.  Those cells are appended at the end in index order so that callers
// still see every cell exactly once.
inline std::size_t d8_topological_order(const std::int64_t* receivers,
                                        std::size_t n,
                                        std::int64_t* out_order) {
    std::vector<std::int32_t> indegree(n, 0);
    for (std::size_t k = 0; k < n; ++k) {
        const std::int64_t r = receivers[k];
        if (r >= 0 && static_cast<std::size_t>(r) < n && static_cast<std::size_t>(r) != k) {
            ++indegree[static_cast<std::size_t>(r)];
        }
    }

    std::vector<std::int64_t> stack;
    stack.reserve(n);
    for (std::size_t k = 0; k < n; ++k) {
        if (indegree[k] == 0) stack.push_back(static_cast<std::int64_t>(k));
    }

    std::size_t written = 0;
    while (!stack.empty()) {
        const std::int64_t k = stack.back();
        stack.pop_back();
        out_order[written++] = k;
        const std::int64_t r = receivers[k];
        if (r >= 0 && static_cast<std::size_t>(r) < n && r != k) {
            if (--indegree[static_cast<std::size_t>(r)] == 0) stack.push_back(r);
        }
    }

    const std::size_t unordered = n - written;
    if (unordered > 0) {
        // Cells on a cycle: emit them in index order so the output is still a
        // permutation of all cells.
        std::vector<std::uint8_t> emitted(n, 0);
        for (std::size_t p = 0; p < written; ++p) emitted[static_cast<std::size_t>(out_order[p])] = 1;
        for (std::size_t k = 0; k < n; ++k) {
            if (!emitted[k]) out_order[written++] = static_cast<std::int64_t>(k);
        }
    }
    return unordered;
}

// ---------------------------------------------------------------------------
// Flow accumulation
// ---------------------------------------------------------------------------

// Accumulate ``weights`` downstream.  ``weights`` is normally the per-cell
// area (dx*dy, or the true spherical cell area for geographic grids); pass an
// array of ones to get the number of contributing cells, which is what
// TopoToolbox's ``flowacc`` returns by default.
//
// ``accumulate_from`` is an optional 0/1 mask: cells with a zero entry
// contribute nothing and do not pass their accumulated value downstream.
// This reproduces the ``evaluate_at`` option of the original Python code.
inline void d8_accumulate(const std::int64_t* receivers,
                          const std::int64_t* order,
                          const double* weights,
                          const std::uint8_t* accumulate_from,  // may be nullptr
                          std::size_t n,
                          double* out_accumulation) {
    for (std::size_t k = 0; k < n; ++k) {
        const bool active = (accumulate_from == nullptr) || (accumulate_from[k] != 0);
        out_accumulation[k] = active ? weights[k] : 0.0;
    }
    for (std::size_t p = 0; p < n; ++p) {
        const std::int64_t k = order[p];
        if (accumulate_from != nullptr && accumulate_from[k] == 0) continue;
        const std::int64_t r = receivers[k];
        if (r < 0 || static_cast<std::size_t>(r) >= n || r == k) continue;
        out_accumulation[r] += out_accumulation[k];
    }
}

// ---------------------------------------------------------------------------
// Flow length
// ---------------------------------------------------------------------------

// Longest upstream flow distance reaching each cell, plus the code of the
// upstream neighbour that supplied it.  ``step_length`` is the length of the
// step from each cell to its receiver.
//
// ``out_from_codes`` records, for each cell, the direction code *pointing back
// at* the donor that lies on the longest path.  ``FlowLength`` uses it to walk
// the main stem back upstream.
//
// When two donors tie on length the one with the lower flat index wins.  That
// tie-break matters: without it the "main stem" would depend on which valid
// topological order happened to be used, and two runs of the same analysis
// could disagree.
inline void d8_flow_length(const std::int64_t* receivers,
                           const std::int64_t* order,
                           const double* step_length,
                           std::size_t ny,
                           std::size_t nx,
                           double* out_length,
                           std::uint8_t* out_from_codes) {
    const std::size_t n = ny * nx;
    const std::int64_t cols = static_cast<std::int64_t>(nx);

    std::vector<std::int64_t> from_index(n, -1);
    for (std::size_t k = 0; k < n; ++k) out_length[k] = 0.0;

    for (std::size_t p = 0; p < n; ++p) {
        const std::int64_t k = order[p];
        const std::int64_t r = receivers[k];
        if (r < 0 || static_cast<std::size_t>(r) >= n || r == k) continue;
        const double candidate = out_length[k] + step_length[k];
        const std::int64_t incumbent = from_index[static_cast<std::size_t>(r)];
        const bool better = candidate > out_length[r] ||
                            (candidate == out_length[r] && (incumbent < 0 || k < incumbent));
        if (better) {
            out_length[r] = candidate;
            from_index[static_cast<std::size_t>(r)] = k;
        }
    }

    if (!out_from_codes) return;
    for (std::size_t r = 0; r < n; ++r) {
        out_from_codes[r] = 0;
        const std::int64_t k = from_index[r];
        if (k < 0) continue;
        const std::int64_t di = (k / cols) - (static_cast<std::int64_t>(r) / cols);
        const std::int64_t dj = (k % cols) - (static_cast<std::int64_t>(r) % cols);
        for (int d = 0; d < 8; ++d) {
            if (kD8Di[d] == di && kD8Dj[d] == dj) {
                out_from_codes[r] = static_cast<std::uint8_t>(1u << d);
                break;
            }
        }
    }
}

// Propagate a value downstream, but only along the longest-flow-length path.
// ``Relief`` uses mode=Carry (each cell inherits the value of its main-stem
// donor) and ``Ksi`` uses mode=Sum (each cell adds its main-stem donor's
// value).  ``gate`` optionally suppresses the transfer for cells where the
// gate is zero, which is how ``Relief`` implements its ``Ao`` cut-off.
enum class PropagateMode { Carry = 0, Sum = 1 };

inline void d8_propagate_along_main_stem(const std::int64_t* receivers,
                                         const std::int64_t* order,
                                         const std::uint8_t* main_stem_from_codes,
                                         const std::uint8_t* gate,  // may be nullptr
                                         std::size_t ny,
                                         std::size_t nx,
                                         PropagateMode mode,
                                         double* values) {
    const std::size_t n = ny * nx;
    const std::int64_t cols = static_cast<std::int64_t>(nx);

    for (std::size_t p = 0; p < n; ++p) {
        const std::int64_t k = order[p];
        const std::int64_t r = receivers[k];
        if (r < 0 || static_cast<std::size_t>(r) >= n || r == k) continue;

        // Is k the donor that defines r's longest flow path?
        const std::int64_t di = (k / cols) - (r / cols);
        const std::int64_t dj = (k % cols) - (r % cols);
        std::uint8_t code = 0;
        for (int d = 0; d < 8; ++d) {
            if (kD8Di[d] == di && kD8Dj[d] == dj) {
                code = static_cast<std::uint8_t>(1u << d);
                break;
            }
        }
        if (code == 0 || main_stem_from_codes[r] != code) continue;
        if (gate != nullptr && gate[k] == 0) continue;

        if (mode == PropagateMode::Sum) {
            values[r] += values[k];
        } else {
            values[r] = values[k];
        }
    }
}

// ---------------------------------------------------------------------------
// Upstream traversal
// ---------------------------------------------------------------------------

// Mark every cell that drains (directly or indirectly) into any of
// ``outlets``.  Iterative, so basin size is limited only by memory.
inline void d8_upstream_mask(const std::int64_t* receivers,
                             std::size_t n,
                             const std::int64_t* outlets,
                             std::size_t n_outlets,
                             std::uint8_t* out_mask) {
    // Build the donor adjacency in CSR form so the traversal is O(n).
    std::vector<std::int64_t> counts(n + 1, 0);
    for (std::size_t k = 0; k < n; ++k) {
        const std::int64_t r = receivers[k];
        if (r >= 0 && static_cast<std::size_t>(r) < n && static_cast<std::size_t>(r) != k) {
            ++counts[static_cast<std::size_t>(r) + 1];
        }
    }
    for (std::size_t k = 0; k < n; ++k) counts[k + 1] += counts[k];
    std::vector<std::int64_t> donors(static_cast<std::size_t>(counts[n]));
    std::vector<std::int64_t> cursor(counts.begin(), counts.end() - 1);
    for (std::size_t k = 0; k < n; ++k) {
        const std::int64_t r = receivers[k];
        if (r >= 0 && static_cast<std::size_t>(r) < n && static_cast<std::size_t>(r) != k) {
            donors[static_cast<std::size_t>(cursor[static_cast<std::size_t>(r)]++)] =
                static_cast<std::int64_t>(k);
        }
    }

    std::fill(out_mask, out_mask + n, static_cast<std::uint8_t>(0));
    std::vector<std::int64_t> stack;
    for (std::size_t o = 0; o < n_outlets; ++o) {
        const std::int64_t k = outlets[o];
        if (k < 0 || static_cast<std::size_t>(k) >= n) continue;
        if (out_mask[k]) continue;
        out_mask[k] = 1;
        stack.push_back(k);
    }
    while (!stack.empty()) {
        const std::int64_t k = stack.back();
        stack.pop_back();
        for (std::int64_t p = counts[static_cast<std::size_t>(k)];
             p < counts[static_cast<std::size_t>(k) + 1]; ++p) {
            const std::int64_t donor = donors[static_cast<std::size_t>(p)];
            if (out_mask[donor]) continue;
            out_mask[donor] = 1;
            stack.push_back(donor);
        }
    }
}

// ---------------------------------------------------------------------------
// Chi
// ---------------------------------------------------------------------------

// Integrate the chi coordinate upstream from a set of outlets:
//
//     chi(cell) = chi(receiver) + integrand * step_length(cell)
//
// with ``integrand = (A0 / A)^theta``.  ``trapezoid`` switches from the
// left-endpoint rule used by the original TopoAnalysis recursion to the
// trapezoidal rule used by TopoToolbox's ``chitransform``:
//
//     chi(cell) = chi(receiver)
//               + 0.5 * ((A0/A_cell)^theta + (A0/A_recv)^theta) * step_length
//
// Cells outside the outlets' basins, and cells where the drainage area is
// not positive, are left as NaN.  ``max_length`` (if > 0) stops the
// integration once the along-stream distance from the outlet exceeds it.
inline void d8_chi(const std::int64_t* receivers,
                   const std::int64_t* order,
                   const double* area,
                   const double* step_length,
                   std::size_t n,
                   const std::int64_t* outlets,
                   std::size_t n_outlets,
                   double A0,
                   double theta,
                   bool trapezoid,
                   double max_length,
                   const std::uint8_t* mask,  // may be nullptr
                   double* out_chi,
                   double* out_distance) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (std::size_t k = 0; k < n; ++k) {
        out_chi[k] = nan;
        if (out_distance) out_distance[k] = nan;
    }
    for (std::size_t o = 0; o < n_outlets; ++o) {
        const std::int64_t k = outlets[o];
        if (k < 0 || static_cast<std::size_t>(k) >= n) continue;
        if (mask != nullptr && mask[k] == 0) continue;
        out_chi[k] = 0.0;
        if (out_distance) out_distance[k] = 0.0;
    }

    // ``order`` runs donors-before-receivers, so walking it backwards visits
    // every receiver before its donors.
    for (std::size_t p = n; p-- > 0;) {
        const std::int64_t k = order[p];
        if (!std::isnan(out_chi[k])) continue;  // an outlet, already seeded
        const std::int64_t r = receivers[k];
        if (r < 0 || static_cast<std::size_t>(r) >= n || r == k) continue;
        if (std::isnan(out_chi[r])) continue;
        if (mask != nullptr && mask[k] == 0) continue;

        const double a = area[k];
        if (!(a > 0.0)) continue;

        const double dl = step_length[k];
        if (max_length > 0.0 && out_distance) {
            if (out_distance[r] + dl > max_length) continue;
        }

        double integrand = std::pow(A0 / a, theta);
        if (trapezoid) {
            const double ar = area[r];
            const double integrand_r = (ar > 0.0) ? std::pow(A0 / ar, theta) : integrand;
            integrand = 0.5 * (integrand + integrand_r);
        }
        out_chi[k] = out_chi[r] + integrand * dl;
        if (out_distance) out_distance[k] = out_distance[r] + dl;
    }
}

// ---------------------------------------------------------------------------
// Distance to the outlet
// ---------------------------------------------------------------------------

// Along-flow distance from each cell down to the end of its flow path.
// This is TopoToolbox's ``flowdistance(FD,'upstream')``: 0 at every outlet
// and at every cell with no receiver, increasing upstream.
inline void d8_downstream_distance(const std::int64_t* receivers,
                                   const std::int64_t* order,
                                   const double* step_length,
                                   std::size_t n,
                                   double* out_distance) {
    std::fill(out_distance, out_distance + n, 0.0);
    // Receivers before donors, so a cell's receiver is final when it is read.
    for (std::size_t p = n; p-- > 0;) {
        const std::int64_t k = order[p];
        const std::int64_t r = receivers[k];
        if (r < 0 || static_cast<std::size_t>(r) >= n || r == k) continue;
        out_distance[k] = out_distance[r] + step_length[k];
    }
}

// ---------------------------------------------------------------------------
// Minima imposition (carving)
// ---------------------------------------------------------------------------

// Lower each receiver until the flow path descends at a gradient of at least
// ``sl``.  This is TopoToolbox's ``imposemin``: unlike depression filling it
// only ever *lowers* cells, so the channel network is carved through
// obstructions rather than drowned by them.
//
// ``elevations`` is modified in place.
inline void d8_imposemin(const std::int64_t* receivers,
                         const std::int64_t* order,
                         const double* step_length,
                         std::size_t n,
                         double sl,
                         double* elevations) {
    for (std::size_t p = 0; p < n; ++p) {
        const std::int64_t k = order[p];
        const std::int64_t r = receivers[k];
        if (r < 0 || static_cast<std::size_t>(r) >= n || r == k) continue;
        const double candidate = elevations[k] - sl * step_length[k];
        if (candidate < elevations[r]) elevations[r] = candidate;
    }
}

// ---------------------------------------------------------------------------
// Stream order
// ---------------------------------------------------------------------------

// Strahler stream order over the cells marked in ``is_stream``.
//
// Order 1 at channel heads; where two or more streams of the highest incoming
// order meet, the order increases by one; a smaller tributary joining a
// larger stream leaves the order unchanged.  Non-stream cells stay 0.
//
// The rule is applied to each junction as a whole -- the largest incoming
// order and how many donors carry it -- rather than donor by donor.  The
// incremental form gives a different answer depending on the sequence the
// donors arrive in (orders {2, 2, 3} yield 4 or 3 according to whether the 3
// is seen first), which would make the result depend on which valid
// topological order the backend produced.
inline void d8_strahler(const std::int64_t* receivers,
                        const std::int64_t* order,
                        const std::uint8_t* is_stream,
                        std::size_t n,
                        std::int32_t* out_order) {
    std::fill(out_order, out_order + n, 0);

    std::vector<std::int32_t> best(n, 0);   // largest incoming order
    std::vector<std::int32_t> ties(n, 0);   // how many donors carry it

    for (std::size_t p = 0; p < n; ++p) {
        const std::int64_t k = order[p];
        if (!is_stream[k]) continue;

        // Donors are all settled by now, so this cell's own order is final.
        out_order[k] = (ties[static_cast<std::size_t>(k)] >= 2)
                           ? best[static_cast<std::size_t>(k)] + 1
                           : std::max(best[static_cast<std::size_t>(k)], 1);

        const std::int64_t r = receivers[k];
        if (r < 0 || static_cast<std::size_t>(r) >= n || r == k) continue;
        if (!is_stream[r]) continue;

        const std::size_t ru = static_cast<std::size_t>(r);
        if (out_order[k] > best[ru]) {
            best[ru] = out_order[k];
            ties[ru] = 1;
        } else if (out_order[k] == best[ru]) {
            ++ties[ru];
        }
    }
}

// Shreve magnitude: the number of channel heads upstream of each cell.
inline void d8_shreve(const std::int64_t* receivers,
                      const std::int64_t* order,
                      const std::uint8_t* is_stream,
                      std::size_t n,
                      std::int32_t* out_order) {
    std::vector<std::uint8_t> has_donor(n, 0);
    for (std::size_t k = 0; k < n; ++k) {
        if (!is_stream[k]) continue;
        const std::int64_t r = receivers[k];
        if (r >= 0 && static_cast<std::size_t>(r) < n && r != static_cast<std::int64_t>(k)
            && is_stream[r]) {
            has_donor[static_cast<std::size_t>(r)] = 1;
        }
    }
    for (std::size_t k = 0; k < n; ++k) {
        out_order[k] = (is_stream[k] && !has_donor[k]) ? 1 : 0;
    }
    for (std::size_t p = 0; p < n; ++p) {
        const std::int64_t k = order[p];
        if (!is_stream[k]) continue;
        const std::int64_t r = receivers[k];
        if (r < 0 || static_cast<std::size_t>(r) >= n || r == k) continue;
        if (!is_stream[r]) continue;
        out_order[r] += out_order[k];
    }
}

// ---------------------------------------------------------------------------
// Drainage basins
// ---------------------------------------------------------------------------

// Label every cell with the index (1-based) of the outlet it drains to.
// Cells that reach no listed outlet are labelled 0.
inline void d8_drainage_basins(const std::int64_t* receivers,
                               const std::int64_t* order,
                               std::size_t n,
                               const std::int64_t* outlets,
                               std::size_t n_outlets,
                               const std::uint8_t* valid,  // may be nullptr
                               std::int32_t* out_labels) {
    std::fill(out_labels, out_labels + n, 0);
    for (std::size_t o = 0; o < n_outlets; ++o) {
        const std::int64_t k = outlets[o];
        if (k < 0 || static_cast<std::size_t>(k) >= n) continue;
        if (valid != nullptr && valid[k] == 0) continue;
        out_labels[k] = static_cast<std::int32_t>(o + 1);
    }
    for (std::size_t p = n; p-- > 0;) {
        const std::int64_t k = order[p];
        if (out_labels[k] != 0) continue;
        if (valid != nullptr && valid[k] == 0) continue;
        const std::int64_t r = receivers[k];
        if (r < 0 || static_cast<std::size_t>(r) >= n || r == k) continue;
        out_labels[k] = out_labels[r];
    }
}

// Label each cell with the id of the terminal cell (pit or edge outlet) that
// its flow path ends at.  This is the ``drainagebasins(FD)`` behaviour of
// TopoToolbox when no outlets are supplied.
//
// A terminal cell only seeds a basin if something drains into it.  A cell
// that is neither a giver nor a receiver -- an isolated no-data cell, say --
// stays 0, which is also what TopoToolbox does (such a cell never appears in
// its edge list).
inline void d8_drainage_basins_all(const std::int64_t* receivers,
                                   const std::int64_t* order,
                                   const std::uint8_t* valid,  // may be nullptr
                                   std::size_t n,
                                   std::int32_t* out_labels) {
    std::fill(out_labels, out_labels + n, 0);

    std::vector<std::uint8_t> has_donor(n, 0);
    for (std::size_t k = 0; k < n; ++k) {
        if (valid != nullptr && valid[k] == 0) continue;
        const std::int64_t r = receivers[k];
        if (r >= 0 && static_cast<std::size_t>(r) < n && r != static_cast<std::int64_t>(k)) {
            has_donor[static_cast<std::size_t>(r)] = 1;
        }
    }

    std::int32_t next_label = 0;
    for (std::size_t k = 0; k < n; ++k) {
        if (valid != nullptr && valid[k] == 0) continue;
        const std::int64_t r = receivers[k];
        const bool terminal =
            (r < 0 || static_cast<std::size_t>(r) >= n || r == static_cast<std::int64_t>(k));
        if (terminal && has_donor[k]) {
            out_labels[k] = ++next_label;
        }
    }
    for (std::size_t p = n; p-- > 0;) {
        const std::int64_t k = order[p];
        if (out_labels[k] != 0) continue;
        if (valid != nullptr && valid[k] == 0) continue;
        const std::int64_t r = receivers[k];
        if (r < 0 || static_cast<std::size_t>(r) >= n || r == k) continue;
        out_labels[k] = out_labels[r];
    }
}

}  // namespace topoanalysis

#endif  // TOPOANALYSIS_FLOW_ROUTING_HPP
