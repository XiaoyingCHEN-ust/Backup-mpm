//! Read material properties
template <unsigned Tdim>
mpm::Rigid<Tdim>::Rigid(unsigned id, const Json& material_properties)
    : Material<Tdim>(id, material_properties) {
  try {
    density_ = material_properties.at("density").template get<double>();
    constraint_dir_ =
        material_properties.at("constraint_dir").template get<std::string>();
    properties_ = material_properties;
  } catch (Json::exception& except) {
    console_->error("Material parameter not set: {} {}\n", except.what(),
                    except.id);
  }
}

//! Compute stress
template <unsigned Tdim>
Eigen::Matrix<double, 6, 1> mpm::Rigid<Tdim>::compute_stress(
    const Vector6d& stress, const Vector6d& dstrain,
    const ParticleBase<Tdim>* ptr, mpm::dense_map* state_vars) {
  Vector6d dstress;
  dstress.setZero();
  return dstress;
}
