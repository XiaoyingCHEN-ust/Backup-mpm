#include <boost/python.hpp>
#include <iostream>
#include <memory>

#include <boost/python/extract.hpp>
#include <boost/python/list.hpp>
#include <sstream>
#include <string>

#ifdef USE_MPI
#include "mpi.h"
#endif
#include "spdlog/spdlog.h"
#include "tbb/task_scheduler_init.h"

#include "io.h"
#include "solvers/mpm.h"

int my(boost::python::list msgs) {

  std::cout << "hello world\n";

  int argc = len(msgs);
  char** argv;
  argv = new char*[argc];

  std::cout << "argc:" << argc << "\n";
  std::string a;
  int l;
  for (int i = 0; i < argc; i++) {

    a = boost::python::extract<std::string>(msgs[i]);
    std::cout << a << "\n";
    char* item = new char[strlen(a.c_str()) + 1];
    strcpy(item, (a).c_str());
    argv[i] = item;
  }

  std::cout << *argv << "\n";
  // #ifdef USE_MPI
  //   // Initialise MPI
  //   MPI_Init(&argc, &argv);
  //   int mpi_rank;
  //   MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  //   // Get number of MPI ranks
  //   int mpi_size;
  //   MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  // #endif

  try {
    // Logger level (trace, debug, info, warn, error, critical, off)
    spdlog::set_level(spdlog::level::trace);

    // Initialise logger
    auto console = spdlog::stdout_color_mt("main");

    // Create an IO object
    auto io = std::make_shared<mpm::IO>(argc, argv);

    // Set TBB threads
    unsigned nthreads = tbb::task_scheduler_init::default_num_threads();
    // If number of TBB threads are positive set to nthreads
    if (io->nthreads() > 0) nthreads = io->nthreads();
    tbb::task_scheduler_init init(nthreads);
  
    std::cout << "TBB nthreads = "<<nthreads <<"\n";
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
    // #ifdef USE_MPI
    //     MPI_Abort(MPI_COMM_WORLD, 1);
    // #endif
    std::terminate();
  }

#ifdef USE_MPI
  MPI_Finalize();
#endif
  return 1;
}
