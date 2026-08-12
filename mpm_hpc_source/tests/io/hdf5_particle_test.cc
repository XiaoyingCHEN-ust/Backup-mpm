#define CATCH_CONFIG_NO_POSIX_SIGNALS
#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include <chrono>
#include <cstddef>
#include <cstdio>
#include <string>

#include "hdf5_particle.h"

namespace {

class TemporaryFile {
 public:
  TemporaryFile() {
    const auto suffix =
        std::chrono::high_resolution_clock::now().time_since_epoch().count();
    path_ = "mpm_hdf5_particle_" + std::to_string(suffix) + ".h5";
  }

  ~TemporaryFile() { std::remove(path_.c_str()); }

  const std::string& path() const { return path_; }

 private:
  std::string path_;
};

class Hdf5File {
 public:
  explicit Hdf5File(hid_t id) : id_(id) {}
  ~Hdf5File() {
    if (id_ >= 0) H5Fclose(id_);
  }

  Hdf5File(const Hdf5File&) = delete;
  Hdf5File& operator=(const Hdf5File&) = delete;

  bool valid() const { return id_ >= 0; }
  hid_t get() const { return id_; }

  herr_t close() {
    if (id_ < 0) return 0;
    const herr_t status = H5Fclose(id_);
    id_ = -1;
    return status;
  }

 private:
  hid_t id_{-1};
};

std::size_t field_index(const std::string& name) {
  for (std::size_t i = 0; i < mpm::hdf5::particle::NFIELDS; ++i) {
    if (name == mpm::hdf5::particle::field_names[i]) return i;
  }
  return mpm::hdf5::particle::NFIELDS;
}

mpm::HDF5Particle legacy_round_trip(unsigned nstate_vars) {
  TemporaryFile temporary;
  const std::string& filename = temporary.path();

  mpm::HDF5Particle source{};
  source.id = 77;
  source.nstate_vars = nstate_vars;
  for (std::size_t state = 0; state < 20; ++state)
    source.svars[state] = 2000. + 0.25 * static_cast<double>(state);

  Hdf5File output(
      H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
  REQUIRE(output.valid());
  REQUIRE(H5TBmake_table(
              "legacy particle state", output.get(), "table",
              mpm::hdf5::particle::LEGACY_NFIELDS, 1,
              mpm::hdf5::particle::dst_size,
              mpm::hdf5::particle::field_names,
              mpm::hdf5::particle::dst_offset,
              mpm::hdf5::particle::field_type, 1, nullptr, 0, &source) >= 0);

  hsize_t fields = 0;
  hsize_t records = 0;
  REQUIRE(H5TBget_table_info(output.get(), "table", &fields, &records) >= 0);
  REQUIRE(fields == mpm::hdf5::particle::LEGACY_NFIELDS);
  REQUIRE(records == 1);
  REQUIRE_NOTHROW(mpm::hdf5::particle::validate_table_metadata(fields, records,
                                                               1));
  REQUIRE(output.close() >= 0);

  Hdf5File input(H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
  REQUIRE(input.valid());
  mpm::HDF5Particle restored{};
  REQUIRE(H5TBread_table(input.get(), "table",
                         mpm::hdf5::particle::dst_size,
                         mpm::hdf5::particle::dst_offset,
                         mpm::hdf5::particle::dst_sizes, &restored) >= 0);
  return restored;
}

}  // namespace

TEST_CASE("HDF5 particle datatype registers all twenty state variables") {
  REQUIRE(mpm::hdf5::particle::NFIELDS == 173);
  REQUIRE(field_index("gamma_b") == 158);
  REQUIRE(field_index("svars_6") == 159);

  for (std::size_t state = 0; state < 20; ++state) {
    const std::string name = "svars_" + std::to_string(state);
    const std::size_t index = field_index(name);
    CAPTURE(state);
    CAPTURE(name);
    REQUIRE(index < mpm::hdf5::particle::NFIELDS);
    REQUIRE(mpm::hdf5::particle::dst_offset[index] ==
            offsetof(mpm::HDF5Particle, svars) + state * sizeof(double));
    REQUIRE(mpm::hdf5::particle::dst_sizes[index] == sizeof(double));
    REQUIRE(mpm::hdf5::particle::field_type[index] == H5T_NATIVE_DOUBLE);
    if (state >= 6) REQUIRE(index >= 159);
  }
}

TEST_CASE("HDF5 particle round trip preserves all twenty state variables") {
  TemporaryFile temporary;
  const std::string& filename = temporary.path();

  mpm::HDF5Particle source{};
  source.id = 42;
  source.nstate_vars = 20;
  for (std::size_t state = 0; state < 20; ++state)
    source.svars[state] = 1000. + 0.125 * static_cast<double>(state);

  Hdf5File output(
      H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
  REQUIRE(output.valid());
  REQUIRE(H5TBmake_table(
              "particle state round trip", output.get(), "table",
              mpm::hdf5::particle::NFIELDS, 1,
              mpm::hdf5::particle::dst_size,
              mpm::hdf5::particle::field_names,
              mpm::hdf5::particle::dst_offset,
              mpm::hdf5::particle::field_type, 1, nullptr, 0, &source) >= 0);
  hsize_t fields = 0;
  hsize_t records = 0;
  REQUIRE(H5TBget_table_info(output.get(), "table", &fields, &records) >= 0);
  REQUIRE(fields == mpm::hdf5::particle::NFIELDS);
  REQUIRE(records == 1);
  REQUIRE(output.close() >= 0);

  Hdf5File input(
      H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
  REQUIRE(input.valid());
  mpm::HDF5Particle restored{};
  REQUIRE(H5TBread_table(input.get(), "table", mpm::hdf5::particle::dst_size,
                         mpm::hdf5::particle::dst_offset,
                         mpm::hdf5::particle::dst_sizes, &restored) >= 0);

  REQUIRE(restored.id == source.id);
  REQUIRE(restored.nstate_vars == source.nstate_vars);
  for (std::size_t state = 0; state < 20; ++state) {
    CAPTURE(state);
    REQUIRE(restored.svars[state] == Approx(source.svars[state]));
  }
  REQUIRE_NOTHROW(mpm::hdf5::particle::prepare_state_variables_for_restart(
      mpm::hdf5::particle::NFIELDS, &restored, 1));
}

TEST_CASE("Legacy 159-field particle tables cannot fake extended history") {
  SECTION("linear elastic checkpoint with no state is accepted") {
    const auto restored = legacy_round_trip(0);
    auto prepared = restored;
    REQUIRE_NOTHROW(mpm::hdf5::particle::prepare_state_variables_for_restart(
        mpm::hdf5::particle::LEGACY_NFIELDS, &prepared, 1));
    REQUIRE(prepared.nstate_vars == 0);
    for (std::size_t state = 0; state < 6; ++state)
      REQUIRE(prepared.svars[state] ==
              Approx(2000. + 0.25 * static_cast<double>(state)));
    for (std::size_t state = 6; state < 20; ++state)
      REQUIRE(prepared.svars[state] == Approx(0.0));
  }

  SECTION("SANISAND checkpoint declaring twenty states is rejected") {
    auto restored = legacy_round_trip(20);
    REQUIRE(restored.nstate_vars == 20);
    REQUIRE_THROWS_WITH(
        mpm::hdf5::particle::prepare_state_variables_for_restart(
            mpm::hdf5::particle::LEGACY_NFIELDS, &restored, 1),
        Catch::Matchers::Contains("refusing a lossy restart"));
  }
}

TEST_CASE("HDF5 particle metadata validation is strict") {
  REQUIRE_NOTHROW(mpm::hdf5::particle::validate_table_metadata(
      mpm::hdf5::particle::LEGACY_NFIELDS, 3, 3));
  REQUIRE_NOTHROW(mpm::hdf5::particle::validate_table_metadata(
      mpm::hdf5::particle::NFIELDS, 3, 3));
  REQUIRE_THROWS(mpm::hdf5::particle::validate_table_metadata(160, 3, 3));
  REQUIRE_THROWS(mpm::hdf5::particle::validate_table_metadata(
      mpm::hdf5::particle::NFIELDS, 2, 3));
}
