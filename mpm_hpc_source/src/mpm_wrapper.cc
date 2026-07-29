#include <boost/python.hpp>
#include <iostream>
#include <memory>

#include <boost/python/extract.hpp>
#include <boost/python/list.hpp>
#include <sstream>
#include <string>

#include "spdlog/spdlog.h"
#include "tbb/task_scheduler_init.h"

#include "boost/python/stl_iterator.hpp"
#include "io.h"
#include "solvers/mpm.h"
#include <Eigen/Dense>
#include <boost/python.hpp>

class pyMPM {
 private:
  // number of rves
  unsigned nrve_;
  unsigned dim_;
  std::shared_ptr<mpm::MPM> mpm_;
  std::vector<unsigned> ids_;
  std::vector<Eigen::MatrixXd> displacement_gradients_;
  std::vector<double> particle_temperatures_;
  std::vector<double> hydrate_saturations_;    
  Eigen::MatrixXd stresses_;
  Eigen::MatrixXd porosities_;
  Eigen::MatrixXd fabrics_CN_;
  Eigen::MatrixXd fabrics_PO_;
  Eigen::MatrixXd rotations_;

 public:
  pyMPM() { std::cout << "Welcome to MPM!\n"; }

  void initialize(boost::python::list msgs) {
    int argc = len(msgs);
    char** argv;
    argv = new char*[argc];
    std::string a;
    for (int i = 0; i < argc; i++) {
      a = boost::python::extract<std::string>(msgs[i]);
      char* item = new char[strlen(a.c_str()) + 1];
      strcpy(item, (a).c_str());
      argv[i] = item;
    }

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

      std::cout << "TBB nthreads = " << nthreads << "\n";
      // Get analysis type
      const std::string analysis = io->analysis_type();

      // Create an MPM analysis
      mpm_ = Factory<mpm::MPM, const std::shared_ptr<mpm::IO>&>::instance()
                 ->create(analysis, std::move(io));
    } catch (std::exception& exception) {
      std::cerr << "MPM main: " << exception.what() << std::endl;
      std::terminate();
    };
  }

  //! Solve the problem by pure MPM
  void solve() { mpm_->solve(); }

  //! Prepreocessing for mpm (multiscale modeling)
  void pre_process() { mpm_->pre_process(); }

  //! Get information
  boost::python::list get_info() {
    boost::python::list list;
    unsigned dim;
    bool resume;
    unsigned checkpoint_step;
    mpm_->get_info(dim, resume, checkpoint_step);

    dim_ = dim;
    list.append(dim);
    list.append(resume);
    list.append(checkpoint_step);

    return list;
  }

  //! Get state
  boost::python::list get_status() {
    double dt;
    unsigned step, nsteps, output_steps;
    boost::python::list list;
    mpm_->get_status(dt, step, nsteps, output_steps);

    list.append(dt);
    list.append(step);
    list.append(nsteps);
    list.append(output_steps);

    return list;
  }

  //! Get deformation information for dem
  void get_deformation_task() { mpm_->get_deformation_task(); }

  //! Send particle id to dem
  boost::python::list send_ids_task() {
    boost::python::list list;
    ids_.clear();
    displacement_gradients_.clear();

    // ids_, displacement_gradients used to store result
    mpm_->send_deformation_task(ids_, displacement_gradients_);

    nrve_ = ids_.size();

    for (auto iter = ids_.begin(); iter != ids_.end(); ++iter) {
      list.append(*iter);
    }
    return list;
  }

  //! Send deformation information to dem
  boost::python::list send_deformations_task() {
    boost::python::list list;
    unsigned dim = displacement_gradients_.at(0).rows();

    for (unsigned i = 0; i < nrve_; i++) {
      for (unsigned j = 0; j < dim; j++) {
        for (unsigned k = 0; k < dim; k++) {
          list.append(displacement_gradients_.at(i)(j, k));
        }
      }
    }
    return list;
  }

  //! Send temperature information to dem
  boost::python::list send_temperatures_task() {
    boost::python::list list;

    particle_temperatures_.clear();
    mpm_->send_temperature_task(ids_, particle_temperatures_);

    for (unsigned i = 0; i < nrve_; i++) {
          list.append(particle_temperatures_.at(i));
        }
    return list;
  }

  //! Send hydrate saturation information to dem
  boost::python::list send_saturations_task() {
    boost::python::list list;

    hydrate_saturations_.clear();
    mpm_->send_saturation_task(ids_, hydrate_saturations_);

    for (unsigned i = 0; i < nrve_; i++) {
          list.append(hydrate_saturations_.at(i));
        }
    return list;
  }  

  //! Set stress for each particle
  void set_stress_task(const boost::python::list& py_stresses, bool increment) {
    stresses_.resize(6 * nrve_, 1);
    stresses_.setZero();
    for (unsigned i = 0; i < 6 * nrve_; i++)
      stresses_(i) = boost::python::extract<double>(py_stresses[i]);

    stresses_.resize(6, nrve_);
    mpm_->set_stress_task(stresses_, increment);
  }

  //! Set porosity for each particle
  void set_porosity_task(const boost::python::list& py_porosities) {
    porosities_.resize(nrve_, 1);
    porosities_.setZero();
    for (unsigned i = 0; i < nrve_; i++)
      porosities_(i) = boost::python::extract<double>(py_porosities[i]);

    mpm_->set_porosity_task(porosities_);
  }

  //! Set fabric for each particle
  void set_fabric_CN_task(const boost::python::list& py_fabrics) {
    fabrics_CN_.resize(nrve_ * dim_ * dim_, 1);
    fabrics_CN_.setZero();
    for (unsigned i = 0; i < nrve_ * dim_ * dim_; i++)
      fabrics_CN_(i) = boost::python::extract<double>(py_fabrics[i]);

    fabrics_CN_.resize(dim_ * dim_, nrve_);
    std::string fabric_type = "CN";
    mpm_->set_fabric_task(fabric_type, fabrics_CN_);
  }

  //! Set fabric for each particle
  void set_fabric_PO_task(const boost::python::list& py_fabrics) {
    fabrics_PO_.resize(nrve_ * dim_ * dim_, 1);
    fabrics_PO_.setZero();
    for (unsigned i = 0; i < nrve_ * dim_ * dim_; i++)
      fabrics_PO_(i) = boost::python::extract<double>(py_fabrics[i]);

    fabrics_PO_.resize(dim_ * dim_, nrve_);
    std::string fabric_type = "PO";
    mpm_->set_fabric_task(fabric_type, fabrics_PO_);
  }

  //! Set rotation for each particle
  void set_rotation_task(const boost::python::list& py_rotations) {
    rotations_.resize(3 * nrve_, 1);
    rotations_.setZero();
    for (unsigned i = 0; i < 3 * nrve_; i++)
      rotations_(i) = boost::python::extract<double>(py_rotations[i]);

    // dimension
    rotations_.resize(3, nrve_);
    mpm_->set_rotation_task(rotations_);
  }

  //! update particle position and so on
  bool update_state_task() {
    bool status = true;
    try {
      mpm_->update_state_task();
    } catch (std::exception& exception) {
      std::cerr << exception.what() << "\n";
      status = false;
    }
    return status;
  }
};

//! Python wrapper for mpm
BOOST_PYTHON_MODULE(lmpm) {
  boost::python::class_<pyMPM>("mpm")
      .def("initialize", &pyMPM::initialize, "Initialise MPM\n")
      .def("solve", &pyMPM::solve)
      .def("pre_process", &pyMPM::pre_process)
      .def("get_info", &pyMPM::get_info)
      .def("get_status", &pyMPM::get_status)
      .def("get_deformation_task", &pyMPM::get_deformation_task)
      .def("send_temperatures_task", &pyMPM::send_temperatures_task)
      .def("send_saturations_task", &pyMPM::send_saturations_task)
      .def("send_ids_task", &pyMPM::send_ids_task)
      .def("send_deformations_task", &pyMPM::send_deformations_task)
      .def("set_stress_task", &pyMPM::set_stress_task,
           boost::python::arg("increment") = false)
      .def("set_porosity_task", &pyMPM::set_porosity_task)
      .def("set_fabric_CN_task", &pyMPM::set_fabric_CN_task)
      .def("set_fabric_PO_task", &pyMPM::set_fabric_PO_task)
      .def("set_rotation_task", &pyMPM::set_rotation_task)
      .def("update_state_task", &pyMPM::update_state_task);
};