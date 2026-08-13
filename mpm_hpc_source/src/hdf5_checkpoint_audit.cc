#include <algorithm>
#include <cmath>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include "hdf5_particle.h"
#include "json.hpp"

// HDF5-only half of the checkpoint/VTP cross-format audit.  The production
// VTP files use VTK's compressed-binary XML encoding; the pipeline's existing
// Python reader decodes that format and joins these JSONL rows by particle ID.
namespace {

using Json = nlohmann::json;

class Hdf5File {
 public:
  explicit Hdf5File(const std::filesystem::path& path)
      : id_(H5Fopen(path.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT)) {
    if (id_ < 0) throw std::runtime_error("cannot open HDF5 file");
  }
  ~Hdf5File() {
    if (id_ >= 0) H5Fclose(id_);
  }
  Hdf5File(const Hdf5File&) = delete;
  Hdf5File& operator=(const Hdf5File&) = delete;
  hid_t get() const { return id_; }

 private:
  hid_t id_{-1};
};

struct Options {
  std::filesystem::path hdf5;
  bool jsonl{false};
  bool allow_legacy{false};
};

Options parse_options(int argc, char** argv) {
  Options options;
  for (int argument = 1; argument < argc; ++argument) {
    const std::string token(argv[argument]);
    if (token == "--hdf5") {
      if (++argument >= argc)
        throw std::invalid_argument("--hdf5 requires a path");
      options.hdf5 = argv[argument];
    } else if (token == "--jsonl") {
      options.jsonl = true;
    } else if (token == "--allow-legacy") {
      options.allow_legacy = true;
    } else if (token == "--help" || token == "-h") {
      std::cout
          << "Usage: hdf5_checkpoint_audit --hdf5 CHECKPOINT.h5 [--jsonl] "
             "[--allow-legacy]\n\n"
             "Default output is one JSON summary. --jsonl appends one compact "
             "particle object per line for a Python VTP cross-check. Formal "
             "new outputs require the 173-field schema; --allow-legacy accepts "
             "the historical 159-field schema for diagnostics only.\n";
      std::exit(0);
    } else {
      throw std::invalid_argument("unknown argument: " + token);
    }
  }
  if (options.hdf5.empty())
    throw std::invalid_argument("--hdf5 CHECKPOINT.h5 is required");
  return options;
}

template <std::size_t Size>
bool finite(const double (&values)[Size]) {
  return std::all_of(std::begin(values), std::end(values),
                     [](double value) { return std::isfinite(value); });
}

Json summary_json(const std::filesystem::path& path, hsize_t fields,
                  hsize_t count, mpm::Index minimum_id,
                  mpm::Index maximum_id, bool legacy) {
  return Json{{"schema", "mpm-hdf5-checkpoint-audit-v1"},
              {"passed", true},
              {"hdf5", std::filesystem::absolute(path).string()},
              {"table", "table"},
              {"table_fields", fields},
              {"particle_count", count},
              {"minimum_id", minimum_id},
              {"maximum_id", maximum_id},
              {"ids_contiguous_unique", true},
              {"formal_schema", !legacy},
              {"legacy_schema", legacy},
              {"compared_fields",
               {"coordinates", "velocities", "displacements", "stresses",
                "porosities", "volumes", "status"}}};
}

Json particle_json(const mpm::HDF5Particle& particle) {
  return Json{{"record_type", "particle"},
              {"id", particle.id},
              {"coordinates",
               {particle.coord_x, particle.coord_y, particle.coord_z}},
              {"displacements",
               {particle.displacement_x, particle.displacement_y,
                particle.displacement_z}},
              {"velocities",
               {particle.velocity_x, particle.velocity_y,
                particle.velocity_z}},
              {"stresses",
               {particle.stress_xx, particle.stress_yy, particle.stress_zz,
                particle.tau_xy, particle.tau_yz, particle.tau_xz}},
              {"porosity", particle.porosity},
              {"volume", particle.volume},
              {"status", particle.status}};
}

}  // namespace

int main(int argc, char** argv) {
  // Do not let routine expected failures print HDF5's low-level stack.  The
  // single stderr message below is stable enough for tests and batch logs.
  H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
  try {
    const Options options = parse_options(argc, argv);
    if (!std::filesystem::is_regular_file(options.hdf5))
      throw std::runtime_error("HDF5 path is not a regular file");
    Hdf5File file(options.hdf5);

    hsize_t field_count = 0;
    hsize_t record_count = 0;
    if (H5TBget_table_info(file.get(), "table", &field_count, &record_count) <
        0)
      throw std::runtime_error("cannot read HDF5 particle table metadata");
    if (record_count == 0)
      throw std::runtime_error("HDF5 particle table has no records");
    const bool legacy = field_count == mpm::hdf5::particle::LEGACY_NFIELDS;
    if (field_count != mpm::hdf5::particle::NFIELDS && !legacy)
      throw std::runtime_error("HDF5 particle table has unsupported field count " +
                               std::to_string(field_count));
    if (legacy && !options.allow_legacy)
      throw std::runtime_error(
          "formal checkpoint requires 173 fields; legacy table has 159 "
          "(pass --allow-legacy only for a diagnostic audit)");

    std::vector<mpm::HDF5Particle> particles(record_count);
    if (H5TBread_table(file.get(), "table", mpm::hdf5::particle::dst_size,
                       mpm::hdf5::particle::dst_offset,
                       mpm::hdf5::particle::dst_sizes, particles.data()) < 0)
      throw std::runtime_error("cannot read HDF5 particle table data");

    std::sort(particles.begin(), particles.end(),
              [](const auto& left, const auto& right) {
                return left.id < right.id;
              });
    for (std::size_t index = 0; index < particles.size(); ++index) {
      const auto& particle = particles[index];
      // The registered mesh uses the canonical dense ID space [0,N).  This
      // simultaneously detects duplicates, gaps, and out-of-range IDs.
      if (particle.id != static_cast<mpm::Index>(index))
        throw std::runtime_error(
            "particle IDs must be unique and contiguous from zero; record " +
            std::to_string(index) + " has ID " +
            std::to_string(particle.id));
      const double coordinates[3] = {particle.coord_x, particle.coord_y,
                                     particle.coord_z};
      const double displacements[3] = {
          particle.displacement_x, particle.displacement_y,
          particle.displacement_z};
      const double velocities[3] = {particle.velocity_x, particle.velocity_y,
                                    particle.velocity_z};
      const double stresses[6] = {particle.stress_xx, particle.stress_yy,
                                  particle.stress_zz, particle.tau_xy,
                                  particle.tau_yz, particle.tau_xz};
      if (!finite(coordinates) || !finite(displacements) ||
          !finite(velocities) || !finite(stresses) ||
          !std::isfinite(particle.porosity) ||
          !std::isfinite(particle.volume))
        throw std::runtime_error("particle " + std::to_string(particle.id) +
                                 " has a non-finite core field");
      if (particle.porosity <= 0. || particle.porosity >= 1.)
        throw std::runtime_error("particle " + std::to_string(particle.id) +
                                 " porosity is outside (0,1)");
      if (particle.volume <= 0.)
        throw std::runtime_error("particle " + std::to_string(particle.id) +
                                 " volume is not positive");
      if (!particle.status)
        throw std::runtime_error("particle " + std::to_string(particle.id) +
                                 " status is false");
      if (particle.nstate_vars > 20)
        throw std::runtime_error("particle " + std::to_string(particle.id) +
                                 " nstate_vars exceeds checkpoint capacity");
      for (unsigned state = 0; state < particle.nstate_vars; ++state)
        if (!std::isfinite(particle.svars[state]))
          throw std::runtime_error(
              "particle " + std::to_string(particle.id) +
              " has a non-finite active constitutive state variable");
    }

    Json summary = summary_json(options.hdf5, field_count, record_count,
                                particles.front().id, particles.back().id,
                                legacy);
    if (legacy)
      summary["warning"] =
          "legacy checkpoint accepted for diagnostics, not formal evidence";
    std::cout << summary.dump() << '\n';
    if (options.jsonl) {
      std::cout << std::setprecision(std::numeric_limits<double>::max_digits10);
      for (const auto& particle : particles)
        std::cout << particle_json(particle).dump() << '\n';
    }
    return 0;
  } catch (const std::exception& exception) {
    std::cerr << "hdf5_checkpoint_audit: " << exception.what() << '\n';
    return 2;
  }
}
