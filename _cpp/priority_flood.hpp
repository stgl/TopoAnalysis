// priority_flood.hpp -- Priority-Flood depression filling
//
// Implements the Priority-Flood family of depression-filling algorithms
// described in
//
//     Barnes, R., Lehman, C., Mulla, D. (2014). "Priority-flood: An optimal
//     depression-filling and watershed-labeling algorithm for digital
//     elevation models." Computers & Geosciences 62, 117-127.
//     doi:10.1016/j.cageo.2013.04.024
//
// Three variants from that paper are provided:
//
//   * Algorithm 1 -- Priority-Flood.  Every cell passes through the priority
//     queue.  O(n log n).  Kept as a reference implementation used to test
//     the faster variants.
//   * Algorithm 2 -- Improved Priority-Flood.  Cells that need no filling are
//     handled through a FIFO ("pit") queue instead of the priority queue,
//     which reduces the practical cost to O(n) for most DEMs.  This is the
//     default.
//   * Algorithm 4 -- Priority-Flood+epsilon.  As Algorithm 2, but filled
//     cells are raised by a small increment so that the filled surface has a
//     strictly monotonic downhill path and contains no flat areas.
//
// The implementation is self-contained (standard library only) and templated
// on the elevation type so it can be instantiated for float and double.
//
// Conventions
// -----------
// * Grids are row-major, ``ny`` rows by ``nx`` columns; index ``k = i*nx + j``.
// * ``NaN`` marks no-data.  No-data cells are never filled and never enqueued.
// * Neighbourhood is D8 (8-connected), matching Barnes et al. and TopoToolbox.

#ifndef TOPOANALYSIS_PRIORITY_FLOOD_HPP
#define TOPOANALYSIS_PRIORITY_FLOOD_HPP

#include <cmath>
#include <cstdint>
#include <cstddef>
#include <limits>
#include <queue>
#include <stdexcept>
#include <vector>

namespace topoanalysis {

// D8 neighbour offsets.  Order matches the ArcGIS code order used elsewhere in
// TopoAnalysis: E, SE, S, SW, W, NW, N, NE.
constexpr int kD8Di[8] = {0, 1, 1, 1, 0, -1, -1, -1};
constexpr int kD8Dj[8] = {1, 1, 0, -1, -1, -1, 0, 1};
// Centre-to-centre distance multiples for the same ordering.
const double kD8Dist[8] = {1.0, 1.4142135623730951, 1.0, 1.4142135623730951,
                           1.0, 1.4142135623730951, 1.0, 1.4142135623730951};

enum class FillMode {
    Flat = 0,     // Algorithm 2: fill depressions to a flat surface
    Epsilon = 1,  // Algorithm 4: fill with a small upslope increment
};

// A priority-queue entry.  ``order`` makes the ordering total, so ties are
// broken first-in-first-out.  Without it std::priority_queue would resolve
// ties arbitrarily and the filled surface would not be reproducible.
template <typename T>
struct PqCell {
    T z;
    std::int64_t order;
    std::int64_t idx;
};

template <typename T>
struct PqGreater {
    bool operator()(const PqCell<T>& a, const PqCell<T>& b) const {
        if (a.z != b.z) return a.z > b.z;
        return a.order > b.order;
    }
};

// Options controlling a Priority-Flood run.
struct FloodOptions {
    FillMode mode = FillMode::Flat;

    // Epsilon increment.  Interpreted as a *slope*: a cell is raised to
    //     parent_elevation + epsilon * distance_to_parent
    // where the distance is in the same units as ``cellsize``.  A value of 0
    // with ``mode == Epsilon`` falls back to the smallest representable
    // increment (std::nextafter), which is what Barnes et al. Algorithm 4
    // prescribes.
    double epsilon = 0.0;
    double cellsize = 1.0;

    // Depressions deeper than this are left unfilled (0 or negative disables
    // the check).  Cells inside such a depression keep their original
    // elevation but are still traversed, so the rest of the grid is reached.
    double max_pit_depth = 0.0;

    // Record which cells the flood reached.
    bool track_visited = false;

    // Record which cells a depth limit left unfilled.  Without this the
    // caller cannot tell a cell the flood never reached from one it reached
    // but deliberately declined to raise.
    bool track_unfilled = false;
};

struct FloodResult {
    std::vector<std::uint8_t> visited;   // empty unless options.track_visited
    std::vector<std::uint8_t> unfilled;  // empty unless options.track_unfilled
    std::size_t cells_filled = 0;
    double max_fill_depth = 0.0;
};

namespace detail {

template <typename T>
inline bool is_nodata(T v) {
    return std::isnan(v);
}

// Raise ``value`` by the smallest amount that is representable at that
// magnitude.  Used by Algorithm 4 when no explicit epsilon is supplied.
template <typename T>
inline T next_up(T value) {
    return std::nextafter(value, std::numeric_limits<T>::infinity());
}

}  // namespace detail

// ---------------------------------------------------------------------------
// Algorithm 2 of Barnes et al. (2014), with the epsilon variant of Algorithm 4
// folded in.  ``elevations`` is modified in place.
//
// ``closed`` must be a pre-initialised ny*nx byte array: a non-zero entry
// marks a cell that must never be processed (no-data, or outside a caller
// supplied mask).  ``seeds`` lists the flat indices the flood starts from;
// they are usually the grid perimeter plus every valid cell adjacent to
// no-data, but callers may supply outlet cells instead.
// ---------------------------------------------------------------------------
template <typename T>
FloodResult priority_flood(T* elevations,
                           std::size_t ny,
                           std::size_t nx,
                           std::vector<std::uint8_t>& closed,
                           const std::vector<std::int64_t>& seeds,
                           const FloodOptions& options) {
    if (ny == 0 || nx == 0) return FloodResult{};
    const std::size_t n = ny * nx;
    if (closed.size() != n) {
        throw std::invalid_argument("closed mask has the wrong size");
    }

    FloodResult result;
    if (options.track_visited) result.visited.assign(n, 0);
    if (options.track_unfilled) result.unfilled.assign(n, 0);

    std::priority_queue<PqCell<T>, std::vector<PqCell<T>>, PqGreater<T>> open;
    std::queue<std::int64_t> pit;

    std::int64_t order = 0;
    for (std::int64_t k : seeds) {
        if (k < 0 || static_cast<std::size_t>(k) >= n) continue;
        if (closed[static_cast<std::size_t>(k)]) continue;
        if (detail::is_nodata(elevations[k])) continue;
        closed[static_cast<std::size_t>(k)] = 1;
        open.push(PqCell<T>{elevations[k], order++, k});
    }

    const bool use_epsilon = (options.mode == FillMode::Epsilon);
    const bool limit_depth = options.max_pit_depth > 0.0;
    // The FIFO shortcut is only sound when every cell it holds sits at exactly
    // the current flood level.  Two options break that invariant: an epsilon
    // fill raises each cell by its own increment, and a depth limit re-queues
    // unfilled cells *below* the current level.  In either case the priority
    // queue handles everything, i.e. the run degrades to Algorithm 1.
    const bool use_pit = !use_epsilon && !limit_depth;

    while (!open.empty() || !pit.empty()) {
        std::int64_t cur;
        T cur_z;
        // Draining the FIFO first is what makes Algorithm 2 faster than
        // Algorithm 1 while producing an identical surface: cells already at
        // or below the current flood level cannot be reordered by the
        // priority queue.  When both queues hold the same elevation the
        // priority queue wins, preserving Algorithm 1's visit order.
        if (!pit.empty() && !open.empty() && open.top().z <= elevations[pit.front()]) {
            cur = open.top().idx;
            cur_z = open.top().z;
            open.pop();
        } else if (!pit.empty()) {
            cur = pit.front();
            pit.pop();
            cur_z = elevations[cur];
        } else {
            cur = open.top().idx;
            cur_z = open.top().z;
            open.pop();
        }

        if (options.track_visited) result.visited[static_cast<std::size_t>(cur)] = 1;

        const std::size_t i = static_cast<std::size_t>(cur) / nx;
        const std::size_t j = static_cast<std::size_t>(cur) % nx;

        for (int d = 0; d < 8; ++d) {
            const std::int64_t ni = static_cast<std::int64_t>(i) + kD8Di[d];
            const std::int64_t nj = static_cast<std::int64_t>(j) + kD8Dj[d];
            if (ni < 0 || nj < 0 || static_cast<std::size_t>(ni) >= ny ||
                static_cast<std::size_t>(nj) >= nx) {
                continue;
            }
            const std::int64_t nk = ni * static_cast<std::int64_t>(nx) + nj;
            if (closed[static_cast<std::size_t>(nk)]) continue;
            if (detail::is_nodata(elevations[nk])) {
                closed[static_cast<std::size_t>(nk)] = 1;
                continue;
            }
            closed[static_cast<std::size_t>(nk)] = 1;

            T target;
            if (use_epsilon) {
                if (options.epsilon > 0.0) {
                    target = static_cast<T>(cur_z + options.epsilon * options.cellsize * kD8Dist[d]);
                    if (target <= cur_z) target = detail::next_up(cur_z);
                } else {
                    target = detail::next_up(cur_z);
                }
            } else {
                target = cur_z;
            }

            if (elevations[nk] < target) {
                const double depth = static_cast<double>(target) - static_cast<double>(elevations[nk]);
                if (limit_depth && depth > options.max_pit_depth) {
                    // Leave this depression alone, but keep traversing so the
                    // terrain beyond it is still reached.
                    if (options.track_unfilled) result.unfilled[static_cast<std::size_t>(nk)] = 1;
                    open.push(PqCell<T>{elevations[nk], order++, nk});
                    continue;
                }
                elevations[nk] = target;
                ++result.cells_filled;
                if (depth > result.max_fill_depth) result.max_fill_depth = depth;
                if (use_pit) {
                    pit.push(nk);
                } else {
                    open.push(PqCell<T>{elevations[nk], order++, nk});
                }
            } else {
                open.push(PqCell<T>{elevations[nk], order++, nk});
            }
        }
    }

    return result;
}

// ---------------------------------------------------------------------------
// Algorithm 1 of Barnes et al. (2014).  Every cell goes through the priority
// queue.  Retained because it is the simplest statement of the method and
// therefore a good oracle for testing the improved variant.
// ---------------------------------------------------------------------------
template <typename T>
FloodResult priority_flood_original(T* elevations,
                                    std::size_t ny,
                                    std::size_t nx,
                                    std::vector<std::uint8_t>& closed,
                                    const std::vector<std::int64_t>& seeds,
                                    const FloodOptions& options) {
    if (ny == 0 || nx == 0) return FloodResult{};
    const std::size_t n = ny * nx;
    if (closed.size() != n) {
        throw std::invalid_argument("closed mask has the wrong size");
    }

    FloodResult result;
    if (options.track_visited) result.visited.assign(n, 0);
    if (options.track_unfilled) result.unfilled.assign(n, 0);

    std::priority_queue<PqCell<T>, std::vector<PqCell<T>>, PqGreater<T>> open;
    std::int64_t order = 0;
    for (std::int64_t k : seeds) {
        if (k < 0 || static_cast<std::size_t>(k) >= n) continue;
        if (closed[static_cast<std::size_t>(k)]) continue;
        if (detail::is_nodata(elevations[k])) continue;
        closed[static_cast<std::size_t>(k)] = 1;
        open.push(PqCell<T>{elevations[k], order++, k});
    }

    const bool use_epsilon = (options.mode == FillMode::Epsilon);
    const bool limit_depth = options.max_pit_depth > 0.0;

    while (!open.empty()) {
        const PqCell<T> c = open.top();
        open.pop();
        if (options.track_visited) result.visited[static_cast<std::size_t>(c.idx)] = 1;

        const std::size_t i = static_cast<std::size_t>(c.idx) / nx;
        const std::size_t j = static_cast<std::size_t>(c.idx) % nx;

        for (int d = 0; d < 8; ++d) {
            const std::int64_t ni = static_cast<std::int64_t>(i) + kD8Di[d];
            const std::int64_t nj = static_cast<std::int64_t>(j) + kD8Dj[d];
            if (ni < 0 || nj < 0 || static_cast<std::size_t>(ni) >= ny ||
                static_cast<std::size_t>(nj) >= nx) {
                continue;
            }
            const std::int64_t nk = ni * static_cast<std::int64_t>(nx) + nj;
            if (closed[static_cast<std::size_t>(nk)]) continue;
            if (detail::is_nodata(elevations[nk])) {
                closed[static_cast<std::size_t>(nk)] = 1;
                continue;
            }
            closed[static_cast<std::size_t>(nk)] = 1;

            T target;
            if (use_epsilon) {
                if (options.epsilon > 0.0) {
                    target = static_cast<T>(c.z + options.epsilon * options.cellsize * kD8Dist[d]);
                    if (target <= c.z) target = detail::next_up(c.z);
                } else {
                    target = detail::next_up(c.z);
                }
            } else {
                target = c.z;
            }

            if (elevations[nk] < target) {
                const double depth = static_cast<double>(target) - static_cast<double>(elevations[nk]);
                if (limit_depth && depth > options.max_pit_depth) {
                    if (options.track_unfilled) result.unfilled[static_cast<std::size_t>(nk)] = 1;
                    open.push(PqCell<T>{elevations[nk], order++, nk});
                    continue;
                }
                elevations[nk] = target;
                ++result.cells_filled;
                if (depth > result.max_fill_depth) result.max_fill_depth = depth;
            }
            open.push(PqCell<T>{elevations[nk], order++, nk});
        }
    }

    return result;
}

// ---------------------------------------------------------------------------
// Default seed set: the grid perimeter plus every valid cell that touches a
// no-data cell.  This matches TopoToolbox's ``fillsinks``, whose marker image
// keeps the original elevation on the perimeter and around NaN holes.
// ---------------------------------------------------------------------------
template <typename T>
std::vector<std::int64_t> default_seeds(const T* elevations,
                                        std::size_t ny,
                                        std::size_t nx) {
    std::vector<std::int64_t> seeds;
    if (ny == 0 || nx == 0) return seeds;

    const std::int64_t rows = static_cast<std::int64_t>(ny);
    const std::int64_t cols = static_cast<std::int64_t>(nx);

    for (std::int64_t i = 0; i < rows; ++i) {
        for (std::int64_t j = 0; j < cols; ++j) {
            const std::int64_t k = i * cols + j;
            if (detail::is_nodata(elevations[k])) continue;

            bool is_seed = (i == 0 || j == 0 || i + 1 == rows || j + 1 == cols);
            for (int d = 0; d < 8 && !is_seed; ++d) {
                const std::int64_t ni = i + kD8Di[d];
                const std::int64_t nj = j + kD8Dj[d];
                if (ni < 0 || nj < 0 || ni >= rows || nj >= cols) continue;
                if (detail::is_nodata(elevations[ni * cols + nj])) is_seed = true;
            }
            if (is_seed) seeds.push_back(k);
        }
    }
    return seeds;
}

}  // namespace topoanalysis

#endif  // TOPOANALYSIS_PRIORITY_FLOOD_HPP
