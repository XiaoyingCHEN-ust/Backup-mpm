#define CATCH_CONFIG_NO_POSIX_SIGNALS
#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstddef>
#include <string>
#include <thread>
#include <vector>

#include "parallelism.h"

#include <tbb/parallel_for.h>

#include "io.h"

namespace {

void update_peak(std::atomic<unsigned>& peak, unsigned value) {
  unsigned observed = peak.load(std::memory_order_relaxed);
  while (observed < value &&
         !peak.compare_exchange_weak(observed, value,
                                     std::memory_order_relaxed)) {
  }
}

}  // namespace

TEST_CASE("-p maps through IO to the TBB parallelism setting") {
  std::vector<std::string> arguments{
      "mpm", "-f", MPM_TEST_FIXTURE_DIR, "-i", "minimal_io.json", "-p",
      "4"};
  std::vector<char*> argv;
  argv.reserve(arguments.size());
  for (auto& argument : arguments) argv.push_back(&argument[0]);

  const mpm::IO io(static_cast<int>(argv.size()), argv.data());
  REQUIRE(io.nthreads() == 4);

  const auto parameter = tbb::global_control::max_allowed_parallelism;
  const std::size_t previous_limit =
      tbb::global_control::active_value(parameter);
  mpm::ParallelismGuard parallelism(io.nthreads());
  REQUIRE(tbb::global_control::active_value(parameter) ==
          std::min<std::size_t>(previous_limit, 4));
}

TEST_CASE("ParallelismGuard limits and restores TBB parallelism") {
  constexpr unsigned requested_threads = 4;
  const auto parameter = tbb::global_control::max_allowed_parallelism;
  const std::size_t previous_limit =
      tbb::global_control::active_value(parameter);
  const std::size_t expected_limit =
      std::min<std::size_t>(previous_limit, requested_threads);

#ifdef _OPENMP
  const int previous_openmp_threads = omp_get_max_threads();
#endif

  std::atomic<unsigned> active{0};
  std::atomic<unsigned> peak{0};
  {
    mpm::ParallelismGuard parallelism(requested_threads);
    REQUIRE(tbb::global_control::active_value(parameter) == expected_limit);
#ifdef _OPENMP
    REQUIRE(omp_get_max_threads() == static_cast<int>(requested_threads));
#endif

    // Sleeping while counted as active makes oversubscription observable
    // without relying on process thread counts, CPU sampling, or timing.
    tbb::parallel_for(std::size_t{0}, std::size_t{128}, [&](std::size_t) {
      const unsigned now = active.fetch_add(1, std::memory_order_relaxed) + 1;
      update_peak(peak, now);
      std::this_thread::sleep_for(std::chrono::milliseconds(5));
      active.fetch_sub(1, std::memory_order_relaxed);
    });

    REQUIRE(peak.load(std::memory_order_relaxed) >= 1);
    REQUIRE(peak.load(std::memory_order_relaxed) <= expected_limit);
  }

  REQUIRE(tbb::global_control::active_value(parameter) == previous_limit);
#ifdef _OPENMP
  REQUIRE(omp_get_max_threads() == previous_openmp_threads);
#endif
}
