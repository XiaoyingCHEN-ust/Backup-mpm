#include <iostream>
#include <memory>

#ifdef USE_MPI
#include "mpi.h"
#endif
#include "spdlog/spdlog.h"

#include "parallelism.h"
#include "io.h"
#include "mpm.h"

int main(int argc, char** argv) {

#ifdef USE_MPI
  // Initialise MPI
  MPI_Init(&argc, &argv);
  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  // Get number of MPI ranks
  int mpi_size;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  // Allocate enough space to issue the buffered send
  int mpi_buffer_size = 2000000000;
  void* mpi_buffer = malloc(mpi_buffer_size);
  // Pass the buffer allocated to MPI so it uses it when we issue MPI_Bsend
  MPI_Buffer_attach(mpi_buffer, mpi_buffer_size);

#endif

  try {
    // Logger level (trace, debug, info, warn, error, critical, off)
    spdlog::set_level(spdlog::level::trace);

    // Initialise logger
    auto console = spdlog::stdout_color_mt("main");

    // Create an IO object
    auto io = std::make_shared<mpm::IO>(argc, argv);

    // Keep both runtime controls alive for the entire solve.  This bounds the
    // OpenMP and TBB work used by MPM; unrelated third-party pools may still
    // own idle/helper threads.
    mpm::ParallelismGuard parallelism(io->nthreads());

    // Get analysis type
    const std::string analysis = io->analysis_type();

    // Create an MPM analysis
    auto mpm =
        Factory<mpm::MPM, const std::shared_ptr<mpm::IO>&>::instance()->create(
            analysis, std::move(io));
    // Solve
    mpm->solve();

  } catch (std::exception& exception) {
    std::cerr << "MPM main: " << exception.what() << std::endl;
#ifdef USE_MPI
    free(mpi_buffer);
    MPI_Buffer_detach(&mpi_buffer, &mpi_buffer_size);
    MPI_Abort(MPI_COMM_WORLD, 1);
#endif
    std::terminate();
  }

#ifdef USE_MPI
  free(mpi_buffer);
  MPI_Buffer_detach(&mpi_buffer, &mpi_buffer_size);
  MPI_Finalize();
#endif
}
