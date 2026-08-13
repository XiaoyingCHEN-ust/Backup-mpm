#ifndef MPM_MATERIAL_MATERIAL_H_
#define MPM_MATERIAL_MATERIAL_H_

#include <limits>

#include "Eigen/Dense"
#include "json.hpp"

#include "factory.h"
#include "logger.h"
#include "map.h"
#include "particle.h"
#include "particles/particle_base.h"

// JSON
using Json = nlohmann::json;

namespace mpm {

// Forward declaration of ParticleBase
template <unsigned Tdim>
class ParticleBase;

//! Material base class
//! \brief Base class that stores the information about materials
//! \details Material class stresses and strains
//! \tparam Tdim Dimension
template <unsigned Tdim>
class Material {
 public:
  //! Define a vector of 6 dof
  using Vector6d = Eigen::Matrix<double, 6, 1>;
  //! Define a Matrix of 6 x 6
  using Matrix6x6 = Eigen::Matrix<double, 6, 6>;

  // Constructor with id
  //! \param[in] id Material id
  Material(unsigned id, const Json& material_properties) : id_{id} {
    //! Logger
    std::string logger = "material::" + std::to_string(id);
    console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);
  }

  //! Destructor
  virtual ~Material(){};

  //! Delete copy constructor
  Material(const Material&) = delete;

  //! Delete assignement operator
  Material& operator=(const Material&) = delete;

  //! Return id of the material
  unsigned id() const { return id_; }

  //! Get material property
  //! \tparam Ttype Return type for proerpty
  //! \param[in] key Material property key
  //! \retval result Value of material property
  template <typename Ttype>
  Ttype property(const std::string& key);

  //! Get optional material property
  //! \tparam Ttype Return type for property
  //! \param[in] key Material property key
  //! \param[in] default_value Value returned when key is absent or null
  //! \retval result Value of material property or default value
  template <typename Ttype>
  Ttype property_or(const std::string& key, const Ttype& default_value);

  //! Initialise history variables
  virtual mpm::dense_map initialise_state_variables() = 0;

  //! Initialise history variables from restored particle data
  //! \param[in] porosity Restored particle porosity
  //! \details Materials must explicitly opt in to cross-model checkpoint
  //! reinitialisation. The default rejects the operation.
  virtual mpm::dense_map initialise_state_variables_from_particle(
      double porosity) {
    static_cast<void>(porosity);
    throw std::runtime_error(
        "Material does not support state reinitialisation from particle data");
  }

  //! Initialise history variables from restored particle data and stress
  //! \param[in] porosity Restored particle porosity
  //! \param[in] stress Restored effective stress in MPM Voigt ordering
  //! \details The stress-aware overload lets pressure-dependent materials
  //! place their internal yield surface consistently around a checkpointed
  //! stress state. Materials that do not need the stress retain the legacy
  //! porosity-only behaviour.
  virtual mpm::dense_map initialise_state_variables_from_particle(
      double porosity, const Vector6d& stress) {
    static_cast<void>(stress);
    return this->initialise_state_variables_from_particle(porosity);
  }

  //! Compute stress
  //! \param[in] stress Stress
  //! \param[in] dstrain Strain
  //! \param[in] particle Constant point to particle base
  //! \param[in] state_vars History-dependent state variables
  //! \retval updated_stress Updated value of stress
  virtual Eigen::Matrix<double, 6, 1> compute_stress(const Vector6d& stress,
                                  const Vector6d& dstrain,
                                  const ParticleBase<Tdim>* ptr,
                                  mpm::dense_map* state_vars) = 0;

  // Return plastic work
  virtual double plastic_work() const {
    throw std::runtime_error(
        "Calling the base class function (plastic_work) in "
        "ParticleBase:: illegal operation!");
    return 0;
  };                                  

  // Return drviatoric plastic strain increment 
  virtual double pdstrain() const {
    throw std::runtime_error(
        "Calling the base class function (pdstrain) in "
        "ParticleBase:: illegal operation!");
    return 0;
  };                                  

  // Return yield_type 
  virtual double yield_type() const {
    throw std::runtime_error(
        "Calling the base class function (yield_type) in "
        "ParticleBase:: illegal operation!");
    return 0;
  };  

  // Return cohesion 
  virtual double cohesion() const {
    throw std::runtime_error(
        "Calling the base class function (cohesion) in "
        "ParticleBase:: illegal operation!");
    return 0;
  };  

  // Return phi
  virtual double phi() const {
    throw std::runtime_error(
        "Calling the base class function (phi) in "
        "ParticleBase:: illegal operation!");
    return 0;
  }; 

 protected:
  //! material id
  unsigned id_{std::numeric_limits<unsigned>::max()};
  //! Material properties
  Json properties_;
  //! Logger
  std::unique_ptr<spdlog::logger> console_;
};  // Material class
}  // namespace mpm

#include "material.tcc"

#endif  // MPM_MATERIAL_MATERIAL_H_
