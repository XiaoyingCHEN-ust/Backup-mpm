#include <cstdio>
#include <filesystem>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "hdf5_particle.h"

namespace {

mpm::HDF5Particle valid_particle(mpm::Index id) {
  mpm::HDF5Particle particle{};
  particle.id = id;
  particle.cell_id = id;
  particle.material_id = 0;
  particle.liquid_material_id = 1;
  particle.volume = 1.E-4 + 1.E-6 * id;
  particle.porosity = 0.45 + 0.01 * id;
  particle.coord_x = 0.1 * id;
  particle.coord_y = 0.2 + 0.1 * id;
  particle.coord_z = 0.;
  particle.displacement_x = 1.E-3 * id;
  particle.displacement_y = -2.E-3 * id;
  particle.displacement_z = 0.;
  particle.velocity_x = 3.E-4 * id;
  particle.velocity_y = -4.E-4 * id;
  particle.velocity_z = 0.;
  particle.stress_xx = -1000. - id;
  particle.stress_yy = -2000. - id;
  particle.stress_zz = -1500. - id;
  particle.tau_xy = 25. * id;
  particle.tau_yz = 0.;
  particle.tau_xz = 0.;
  particle.status = true;
  particle.nstate_vars = 2;
  particle.svars[0] = 0.1 * id;
  particle.svars[1] = -0.2 * id;
  return particle;
}

}  // namespace

int main(int argc, char** argv) {
  H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
  if (argc != 3) {
    std::cerr << "usage: hdf5_checkpoint_audit_fixture OUTPUT.h5 "
                 "valid|duplicate|nonfinite|bad_porosity|false_status|"
                 "bad_nstate|nonfinite_svar|legacy\n";
    return 2;
  }
  const std::filesystem::path output = argv[1];
  const std::string mode = argv[2];
  std::vector<mpm::HDF5Particle> particles{
      valid_particle(0), valid_particle(1), valid_particle(2)};
  hsize_t field_count = mpm::hdf5::particle::NFIELDS;
  if (mode == "duplicate") particles[2].id = 1;
  else if (mode == "nonfinite")
    particles[1].velocity_y = std::numeric_limits<double>::quiet_NaN();
  else if (mode == "bad_porosity")
    particles[1].porosity = 1.;
  else if (mode == "false_status")
    particles[1].status = false;
  else if (mode == "bad_nstate")
    particles[1].nstate_vars = 21;
  else if (mode == "nonfinite_svar")
    particles[1].svars[1] = std::numeric_limits<double>::quiet_NaN();
  else if (mode == "legacy")
    field_count = mpm::hdf5::particle::LEGACY_NFIELDS;
  else if (mode != "valid") {
    std::cerr << "unknown fixture mode: " << mode << '\n';
    return 2;
  }

  const hid_t file = H5Fcreate(output.string().c_str(), H5F_ACC_TRUNC,
                               H5P_DEFAULT, H5P_DEFAULT);
  if (file < 0) return 2;
  const herr_t status = H5TBmake_table(
      "checkpoint audit fixture", file, "table", field_count,
      particles.size(), mpm::hdf5::particle::dst_size,
      mpm::hdf5::particle::field_names, mpm::hdf5::particle::dst_offset,
      mpm::hdf5::particle::field_type, particles.size(), nullptr, 0,
      particles.data());
  const herr_t close_status = H5Fclose(file);
  if (status < 0 || close_status < 0) {
    std::remove(output.string().c_str());
    return 2;
  }
  return 0;
}
