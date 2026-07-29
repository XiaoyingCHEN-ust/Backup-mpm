#ifndef MPM_MPM_H_
#define MPM_MPM_H_

#include <chrono>
#include <memory>
#include <vector>
#include <functional>

#include <boost/lexical_cast.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

#include "io.h"
#include "io_mesh.h"
#include "io_mesh_ascii.h"
#include "mesh.h"

#ifdef USE_VTK
#include "vtk_writer.h"
#endif

#ifdef USE_PARTIO
#include "partio_writer.h"
#endif

namespace mpm {
//! MPM class
//! \brief MPM class calls solver and algorithm
//! \details MPM class: implicit and explicit MPM
class MPM {
 public:
  //! Constructor
  MPM(const std::shared_ptr<IO>& io) : io_(io) {

    analysis_ = io_->analysis();

    // Unique id
    if (analysis_.find("uuid") != analysis_.end())
      uuid_ = analysis_["uuid"].template get<std::string>();

    if (uuid_.empty())
      uuid_ =
          boost::lexical_cast<std::string>(boost::uuids::random_generator()());
  }

  // Initialise mesh and particles
  virtual bool initialise_mesh() = 0;

  // Initialise particles
  virtual bool initialise_particles() = 0;

  // Initialise materials
  virtual bool initialise_materials() = 0;

  // Initialise external loads
  virtual bool initialise_loads() = 0;

  // Initialise math functions
  virtual bool initialise_math_functions(const Json&) = 0;

  // Initialise vtk output for single phase
  virtual bool initialise_vtk() = 0;

  // Initialise vtk output for two phase
  virtual bool initialise_vtk_twophase() = 0;

  // Solve
  virtual bool solve() = 0;

  // Check point restart
  virtual bool checkpoint_resume() = 0;

  //! Write HDF5 files
  virtual void write_hdf5(mpm::Index step, mpm::Index max_steps) = 0;

  //! pre processing
  virtual bool pre_process() = 0;

  //! Get analysis information
  virtual void get_info(unsigned& dim, bool& resume,
                        unsigned& checkpoint_step) = 0;

  //! get step, time
  virtual void get_status(double& dt, unsigned& step, unsigned& nsteps,
                          unsigned& output_steps) = 0;

  virtual bool get_deformation_task() = 0;

  //! send deformation task
  virtual bool send_deformation_task(
      std::vector<unsigned>& id,
      std::vector<Eigen::MatrixXd>& displacement_gradients) = 0;

  //! send temperature task
  virtual bool send_temperature_task(
      std::vector<unsigned>& id,
      std::vector<double>& particle_temperature) = 0;        

  //! send saturation task
  virtual bool send_saturation_task(
      std::vector<unsigned>& id,
      std::vector<double>& hydrate_saturation) = 0;  

  //! Set particle stress
  virtual bool set_stress_task(const Eigen::MatrixXd& stresses,
                               bool increment) = 0;

  //! Set particle porosity
  virtual bool set_porosity_task(const Eigen::MatrixXd& porosities) = 0;

  // Set particle fabric
  virtual bool set_fabric_task(std::string fabric_type,
                               const Eigen::MatrixXd& fabrics) = 0;

  // Set particle rotation
  virtual bool set_rotation_task(const Eigen::MatrixXd& rotations) = 0;

  // Update particle state eg. velocity and position
  virtual bool update_state_task() = 0;

  //! Write reaction force
  virtual void write_reaction_force(bool overwrite, mpm::Index step,
                                    mpm::Index max_steps) = 0;

#ifdef USE_VTK
  //! Write VTK files
  virtual void write_vtk(mpm::Index step, mpm::Index max_steps) = 0;
#endif

#ifdef USE_PARTIO
  //! Write PARTIO files
  virtual void write_partio(mpm::Index step, mpm::Index max_steps) = 0;
#endif

 protected:
  //! A unique id for the analysis
  std::string uuid_;
  //! Time step size
  double dt_{std::numeric_limits<double>::max()};
  //! Current step
  mpm::Index step_{0};
  //! Number of steps
  double current_time_{0};
  //! Number of steps  
  mpm::Index nsteps_{std::numeric_limits<mpm::Index>::max()};
  //! Output steps
  mpm::Index output_steps_{std::numeric_limits<mpm::Index>::max()};
  //! A shared ptr to IO object
  std::shared_ptr<mpm::IO> io_;
  //! JSON analysis object
  Json analysis_;
  //! JSON post-process object
  Json post_process_;  
  //! Logger
  std::shared_ptr<spdlog::logger> console_;
};
}  // namespace mpm

#endif  // MPM_MPM_H_
