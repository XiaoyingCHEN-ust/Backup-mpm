#ifndef MPM_PARALLELISM_H_
#define MPM_PARALLELISM_H_

#include <memory>

// global_control was a preview API in older Intel TBB releases.  Defining the
// preview switch is harmless with oneTBB and keeps the compatibility include
// usable on releases that shipped the header before making it fully supported.
#ifndef TBB_PREVIEW_GLOBAL_CONTROL
#define TBB_PREVIEW_GLOBAL_CONTROL 1
#endif
#include <tbb/global_control.h>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace mpm {

//! Scope the solver's OpenMP and TBB parallelism settings.
//! \brief A positive limit controls the OpenMP team size selected by the
//! calling thread and the process-wide TBB scheduler limit.  Independent
//! third-party thread pools are outside this guard.
class ParallelismGuard {
 public:
  explicit ParallelismGuard(unsigned nthreads) {
    if (nthreads == 0) return;

    // Construct this first so a failed TBB allocation cannot leave OpenMP in a
    // partially updated state.
    tbb_parallelism_ = std::make_unique<tbb::global_control>(
        tbb::global_control::max_allowed_parallelism, nthreads);

#ifdef _OPENMP
    previous_openmp_threads_ = omp_get_max_threads();
    omp_set_num_threads(static_cast<int>(nthreads));
    restore_openmp_ = true;
#endif
  }

  ~ParallelismGuard() {
#ifdef _OPENMP
    if (restore_openmp_) omp_set_num_threads(previous_openmp_threads_);
#endif
  }

  ParallelismGuard(const ParallelismGuard&) = delete;
  ParallelismGuard& operator=(const ParallelismGuard&) = delete;

 private:
  std::unique_ptr<tbb::global_control> tbb_parallelism_;
#ifdef _OPENMP
  int previous_openmp_threads_{0};
  bool restore_openmp_{false};
#endif
};

}  // namespace mpm

#endif  // MPM_PARALLELISM_H_
