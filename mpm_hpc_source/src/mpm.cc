#include <memory>

#include "factory.h"
#include "io.h"
#include "solvers/mpm.h"
#include "solvers/mpm_explicit.h"
#include "solvers/mpm_explicit_twophase.h"
#include "solvers/mpm_semi_implicit_twophase.h"
#include "solvers/tm_mpm_explicit.h"
#include "solvers/thm_mpm_explicit_twophase.h"
#include "solvers/thm_mpm_explicit_threephase.h"
#include "solvers/thm_mpm_explicit_threephase_new.h"
#include "solvers/thm_mpm_explicit_threephase_lag.h"
#include "solvers/thm_mpm_explicit_sat_frozen.h"
#include "solvers/thm_mpm_explicit_MHBS.h"
#include "solvers/thm_mpm_semi_implicit_twophase.h"
#include "solvers/thm_mpm_semi_implicit_sat_frozen.h"
#include "solvers/thm_mpm_semi_implicit_unsat_frozen.h"
#include "solvers/thermo_mpm_with_phase_change.h"
#include "solvers/thermo_mpm_implicit.h"
#include "solvers/hydro_mpm_explicit.h"
#include "solvers/hydro_mpm_explicit_MHBS.h"

namespace mpm {
// Stress update method
std::map<std::string, StressUpdate> stress_update = {
    {"usf", StressUpdate::USF},
    {"usl", StressUpdate::USL},
    {"musl", StressUpdate::MUSL}};
}  // namespace mpm

// 2D Explicit MPM
static Register<mpm::MPM, mpm::MPMExplicit<2>, const std::shared_ptr<mpm::IO>&>
    mpm_explicit_2d("MPMExplicit2D");

// 3D Explicit MPM
static Register<mpm::MPM, mpm::MPMExplicit<3>, const std::shared_ptr<mpm::IO>&>
    mpm_explicit_3d("MPMExplicit3D");

// 2D Implicit Thermo-mechnical MPM
static Register<mpm::MPM, mpm::ThermoMPMImplicit<2>, const std::shared_ptr<mpm::IO>&>
    thermo_mpm_implicit_2d("ThermoMPMImplicit2D");

// 3D Implicit Thermal MPM
static Register<mpm::MPM, mpm::ThermoMPMImplicit<3>, const std::shared_ptr<mpm::IO>&>
    thermo_mpm_implicit_3d("ThermoMPMImplicit3D");

// 2D Explicit Thermo-mechnical MPM
static Register<mpm::MPM, mpm::ThermoMPMExplicit<2>, const std::shared_ptr<mpm::IO>&>
    thermo_mpm_explicit_2d("ThermoMPMExplicit2D");

// 3D Explicit Thermo-mechnical MPM
static Register<mpm::MPM, mpm::ThermoMPMExplicit<3>, const std::shared_ptr<mpm::IO>&>
    thermo_mpm_explicit_3d("ThermoMPMExplicit3D");

// 2D Explicit Two Phase MPM
static Register<mpm::MPM, mpm::MPMExplicitTwoPhase<2>,
                const std::shared_ptr<mpm::IO>&>
    mpm_explicit_twophase_2d("MPMExplicitTwoPhase2D");

// 3D Explicit Two Phase MPM
static Register<mpm::MPM, mpm::MPMExplicitTwoPhase<3>,
                const std::shared_ptr<mpm::IO>&>
    mpm_explicit_twophase_3d("MPMExplicitTwoPhase3D");

// 2D Explicit Two Phase Thermo-Hydro-mechnical MPM
static Register<mpm::MPM, mpm::ThermoMPMExplicitTwoPhase<2>,
                const std::shared_ptr<mpm::IO>&>
    thermo_mpm_explicit_twophase_2d("ThermoMPMExplicitTwoPhase2D");

// 3D Explicit Two Phase Thermo-Hydro-mechnical MPM
static Register<mpm::MPM, mpm::ThermoMPMExplicitTwoPhase<3>,
                const std::shared_ptr<mpm::IO>&>
    thermo_mpm_explicit_twophase_3d("ThermoMPMExplicitTwoPhase3D");

// 2D Explicit Three Phase Thermo-Hydro-mechnical MPM
static Register<mpm::MPM, mpm::ThermoMPMExplicitThreePhase<2>,
                const std::shared_ptr<mpm::IO>&>
    thermo_mpm_explicit_threephase_2d("ThermoMPMExplicitThreePhase2D");

// 3D Explicit Three Phase Thermo-Hydro-mechnical MPM
static Register<mpm::MPM, mpm::ThermoMPMExplicitThreePhase<3>,
                const std::shared_ptr<mpm::IO>&>
    thermo_mpm_explicit_threephase_3d("ThermoMPMExplicitThreePhase3D");

// 2D Explicit Three Phase Thermo-Hydro-mechnical MPM
static Register<mpm::MPM, mpm::ThermoMPMExplicitThreePhaseNew<2>,
                const std::shared_ptr<mpm::IO>&>
    thermo_mpm_explicit_threephase_new_2d("ThermoMPMExplicitThreePhaseNew2D");

// 2D Explicit Three Phase Thermo-Hydro-mechnical MPM
static Register<mpm::MPM, mpm::ThermoMPMExplicitThreePhaseLag<2>,
                const std::shared_ptr<mpm::IO>&>
    thermo_mpm_explicit_threephase_lag_2d("ThermoMPMExplicitThreePhaseLag2D");

// 3D Explicit Three Phase Thermo-Hydro-mechnical MPM
static Register<mpm::MPM, mpm::ThermoMPMExplicitThreePhaseNew<3>,
                const std::shared_ptr<mpm::IO>&>
    thermo_mpm_explicit_threephase_new_3d("ThermoMPMExplicitThreePhaseNew3D");

// 3D Explicit Three Phase Thermo-Hydro-mechnical MPM
static Register<mpm::MPM, mpm::ThermoMPMExplicitThreePhaseLag<3>,
                const std::shared_ptr<mpm::IO>&>
    thermo_mpm_explicit_threephase_lag_3d("ThermoMPMExplicitThreePhaseLag3D");

// 2D Explicit THM-MPM for saturated frozen soil
static Register<mpm::MPM, mpm::THMMPMExplicitSatFrozen<2>,
                const std::shared_ptr<mpm::IO>&>
    thm_mpm_explicit_sat_frozen_2d("THMMPMExplicitSatFrozen2D");

// 3D Explicit THM-MPM for saturated frozen soil
static Register<mpm::MPM, mpm::THMMPMExplicitSatFrozen<3>,
                const std::shared_ptr<mpm::IO>&>
    thm_mpm_explicit_sat_frozen_3d("THMMPMExplicitSatFrozen3D");

// 2D SemiImplicit Two Phase MPM
static Register<mpm::MPM, mpm::MPMSemiImplicitTwoPhase<2>,
                const std::shared_ptr<mpm::IO>&>
    mpm_semi_implicit_twophase_2d("MPMSemiImplicitTwoPhase2D");

// 3D SemiImplicit Two Phase MPM
static Register<mpm::MPM, mpm::MPMSemiImplicitTwoPhase<3>,
                const std::shared_ptr<mpm::IO>&>
    mpm_semi_implicit_twophase_3d("MPMSemiImplicitTwoPhase3D");

// 2D SemiImplicit Two Phase Thermo-mechnical MPM
static Register<mpm::MPM, mpm::ThermoMPMSemiImplicitTwoPhase<2>,
                const std::shared_ptr<mpm::IO>&>
    thermo_mpm_semi_implicit_twophase_2d("ThermoMPMSemiImplicitTwoPhase2D");

// 3D SemiImplicit Two Phase Thermo-mechnical MPM
static Register<mpm::MPM, mpm::ThermoMPMSemiImplicitTwoPhase<3>,
                const std::shared_ptr<mpm::IO>&>
    thermo_mpm_semi_implicit_twophase_3d("ThermoMPMSemiImplicitTwoPhase3D");

// 2D SemiImplicit THM-MPM for saturated frozen soil
static Register<mpm::MPM, mpm::THMMPMSemiImplicitSatFrozen<2>,
                const std::shared_ptr<mpm::IO>&>
    thm_mpm_semi_implicit_sat_frozen_2d("THMMPMSemiImplicitSatFrozen2D");

// 3D SemiImplicit THM-MPM for saturated frozen soil
static Register<mpm::MPM, mpm::THMMPMSemiImplicitSatFrozen<3>,
                const std::shared_ptr<mpm::IO>&>
    thm_mpm_semi_implicit_sat_frozen_3d("THMMPMSemiImplicitSatFrozen3D");

// 2D SemiImplicit THM-MPM for unsaturated frozen soil
static Register<mpm::MPM, mpm::THMMPMSemiImplicitUnsatFrozen<2>,
                const std::shared_ptr<mpm::IO>&>
    thm_mpm_semi_implicit_unsat_frozen_2d("THMMPMSemiImplicitUnsatFrozen2D");

// 3D SemiImplicit THM-MPM for unsaturated frozen soil
static Register<mpm::MPM, mpm::THMMPMSemiImplicitUnsatFrozen<3>,
                const std::shared_ptr<mpm::IO>&>
    thm_mpm_semi_implicit_unsat_frozen_3d("THMMPMSemiImplicitUnsatFrozen3D");

// 2D SemiImplicit THM-MPM with phase change
static Register<mpm::MPM, mpm::ThermoMPMPhaseChange<2>,
                const std::shared_ptr<mpm::IO>&>
    thermo_mpm_phase_change_2d("ThermoMPMPhaseChange2D");

// 3D SemiImplicit THM-MPM with phase change
static Register<mpm::MPM, mpm::ThermoMPMPhaseChange<3>,
                const std::shared_ptr<mpm::IO>&>
    thermo_mpm_phase_change_3d("ThermoMPMPhaseChange3D");

// 2D Explicit THM-MPM for MHBS
static Register<mpm::MPM, mpm::THMMPMExplicitMHBS<2>,
                const std::shared_ptr<mpm::IO>&>
    thm_mpm_explicit_MHBS_2d("THMMPMExplicitMHBS2D");

// 3D  Explicit THM-MPM for MHBS
static Register<mpm::MPM, mpm::THMMPMExplicitMHBS<3>,
                const std::shared_ptr<mpm::IO>&>
    thm_mpm_explicit_MHBS_3d("THMMPMExplicitMHBS3D");

// 2D Explicit hydro MPM
static Register<mpm::MPM, mpm::HydroMPMExplicit<2>,
                const std::shared_ptr<mpm::IO>&>
    hydro_mpm_explicit_2d("HydroMPMExplicit2D");

// 3D  Explicit hydro MPM
static Register<mpm::MPM, mpm::HydroMPMExplicit<3>,
                const std::shared_ptr<mpm::IO>&>
    hydro_mpm_explicit_3d("HydroMPMExplicit3D");

// 2D Explicit hydro MHBS
static Register<mpm::MPM, mpm::HydroMPMExplicitMHBS<2>,
                const std::shared_ptr<mpm::IO>&>
    hydro_mpm_explicit_MHBS_2d("HydroMPMExplicitMHBS2D");

// 3D  Explicit hydro MHBS
static Register<mpm::MPM, mpm::HydroMPMExplicitMHBS<3>,
                const std::shared_ptr<mpm::IO>&>
    hydro_mpm_explicit_MHBS_3d("HydroMPMExplicitMHBS3D");