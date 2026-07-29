"""Apply the small source fixes required for reliable HDF5 checkpoint resume."""

from __future__ import annotations

import argparse
from pathlib import Path


PATCHES = (
    (
        Path("include/mesh/mesh.tcc"),
        """  std::vector<HDF5Particle> dst_buf;
  dst_buf.reserve(nparticles);
""",
        """  std::vector<HDF5Particle> dst_buf(nparticles);
""",
        "size the HDF5 read buffer",
    ),
    (
        Path("include/solvers/thm_mpm_explicit_threephase_new.tcc"),
        """  // Check point resume
  if (resume) {
    this->checkpoint_resume();
    this->current_time_ = analysis_["resume"]["current_time"].template get<double>();
    std::cout << "current_time" << this->current_time_ << "\\n";
  }
""",
        """  // Check point resume
  mpm::Index start_step = 0;
  if (resume) {
    if (!this->checkpoint_resume())
      throw std::runtime_error("Checkpoint resume failed");
    this->current_time_ = analysis_["resume"]["current_time"].template get<double>();
    start_step = this->step_;
    std::cout << "current_time" << this->current_time_ << "\\n";
  }
""",
        "retain the resumed global step",
    ),
    (
        Path("include/solvers/thm_mpm_explicit_threephase_new.tcc"),
        """  for (step_ = 0; step_ <= nsteps_; ++step_) {
""",
        """  for (step_ = start_step; step_ <= nsteps_; ++step_) {
""",
        "start the loop at the resumed global step",
    ),
)


def apply_exact_patch(
    path: Path, text: str, old: str, new: str, description: str
) -> tuple[str, bool]:
    if new in text:
        print(f"Already patched: {description} ({path})")
        return text, False
    occurrences = text.count(old)
    if occurrences != 1:
        raise RuntimeError(
            f"Cannot safely patch {description} in {path}: "
            f"expected one source pattern, found {occurrences}"
        )
    print(f"Will patch: {description} ({path})")
    return text.replace(old, new, 1), True


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mpm-source", type=Path, required=True)
    args = parser.parse_args()

    source = args.mpm_source.resolve()
    if not (source / "build").is_dir():
        parser.error(f"MPM build directory is missing: {source / 'build'}")

    original_text = {}
    patched_text = {}
    changed_paths = set()
    for relative_path, old, new, description in PATCHES:
        target = source / relative_path
        if not target.is_file():
            parser.error(f"MPM source file is missing: {target}")
        if target not in original_text:
            original_text[target] = target.read_text(encoding="utf-8")
            patched_text[target] = original_text[target]
        patched_text[target], changed = apply_exact_patch(
            target, patched_text[target], old, new, description
        )
        if changed:
            changed_paths.add(target)

    # Write only after every pattern has passed its preflight check, preventing
    # a partially modified source tree when a later pattern is unrecognised.
    for target in changed_paths:
        target.write_text(patched_text[target], encoding="utf-8")
        print(f"Updated checkpoint-resume source: {target}")


if __name__ == "__main__":
    main()
