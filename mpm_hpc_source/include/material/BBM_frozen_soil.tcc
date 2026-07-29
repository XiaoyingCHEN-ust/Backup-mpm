template <unsigned Tdim>
mpm::Modified_BBM_for_frozen_soil<Tdim>::MohrCoulomb(unsigned id,
                                    const Json& material_properties)
    : Material<Tdim>(id, material_properties) {
  try {
    // General parameters
    // Density
    density_ = material_properties.at("density").template get<double>();
    // Young's modulus
    youngs_modulus_ =
        material_properties.at("youngs_modulus").template get<double>();
    // Poisson ratio
    poisson_ratio_ =
        material_properties.at("poisson_ratio").template get<double>();   
    //  preconsolidation stress
    pn0_star_ = material_properties.at("pn0_star").template get<double>();  
    //  slope of critical state line
    M_ = material_properties.at("M").template get<double>();     
    // material constant controlling the increase of strength with 
    // cryogenic suction
    k = material_properties.at("k").template get<double>();     
    // elastic stiffness parameter for changes in mean net stress
    kappa_ = material_properties.at("kappa").template get<double>();    
    // reference pressure
    pc_ = material_properties.at("pc").template get<double>();
    // slope of virgin state for saturated conditions
    lambda0_ = material_properties.at("lambda0").template get<double>();
    // model constant that controls the rate of increase of soil stiffness
    beta_ = material_properties.at("beta").template get<double>();
    // model parameter that defines the minimum soil compressibility
    r_ = material_properties.at("r").template get<double>();
    // Bulk modulus
    bulk_modulus_ = youngs_modulus_ / (3.0 * (1. - 2. * poisson_ratio_));
    // Shear modulus
    shear_modulus_ = youngs_modulus_ / (2.0 * (1 + poisson_ratio_));   

    // Set elastic tensor
    this->compute_elastic_tensor();
  } catch (std::exception& except) {
    console_->error("Material parameter not set: {}\n", except.what());
  }
}

//! Initialise state variables
template <unsigned Tdim>
mpm::dense_map mpm::Modified_BBM_for_frozen_soil<Tdim>::initialise_state_variables() {
  mpm::dense_map state_vars = {// model parameters

                               // stress variables
                               // net stress
                               {"pn", 0},
                               // cryo suction
                               {"s_cryo", 0},
                               // deviatoric stress
                               {"q", 0}};
  return state_vars;
}

//! Compute elastic tensor
template <unsigned Tdim>
bool mpm::MohrCoulomb<Tdim>::compute_elastic_tensor() {
  // Shear modulus
  const double G = shear_modulus_;
  const double a1 = bulk_modulus_ + (4.0 / 3.0) * G;
  const double a2 = bulk_modulus_ - (2.0 / 3.0) * G;
  // compute elastic stiffness matrix
  // clang-format off
  de_(0,0)=a1;    de_(0,1)=a2;    de_(0,2)=a2;    de_(0,3)=0;    de_(0,4)=0;    de_(0,5)=0;
  de_(1,0)=a2;    de_(1,1)=a1;    de_(1,2)=a2;    de_(1,3)=0;    de_(1,4)=0;    de_(1,5)=0;
  de_(2,0)=a2;    de_(2,1)=a2;    de_(2,2)=a1;    de_(2,3)=0;    de_(2,4)=0;    de_(2,5)=0;
  de_(3,0)= 0;    de_(3,1)= 0;    de_(3,2)= 0;    de_(3,3)=G;    de_(3,4)=0;    de_(3,5)=0;
  de_(4,0)= 0;    de_(4,1)= 0;    de_(4,2)= 0;    de_(4,3)=0;    de_(4,4)=G;    de_(4,5)=0;
  de_(5,0)= 0;    de_(5,1)= 0;    de_(5,2)= 0;    de_(5,3)=0;    de_(5,4)=0;    de_(5,5)=G;
  // clang-format on

  return true;
}

//! Compute yield function and yield state
