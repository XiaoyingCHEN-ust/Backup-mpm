#include "material/material.h"
#include "material/bingham.h"
#include "material/cam_clay.h"
#include "material/linear_elastic.h"
#include "material/mohr_coulomb.h"
#include "material/MC_frozen_soil.h"
#include "material/newtonian.h"
#include "material/norsand.h"
#include "material/rigid.h"
#include "material/sanisand.h"
// Bingham 2D
static Register<mpm::Material<2>, mpm::Bingham<2>, unsigned, const Json&>
    bingham_2d("Bingham2D");

// Bingham 3D
static Register<mpm::Material<3>, mpm::Bingham<3>, unsigned, const Json&>
    bingham_3d("Bingham3D");

// CamClay 2D
static Register<mpm::Material<2>, mpm::CamClay<2>, unsigned, const Json&>
    cam_clay_2d("CamClay2D");

// CamClay 3D
static Register<mpm::Material<3>, mpm::CamClay<3>, unsigned, const Json&>
    cam_clay_3d("CamClay3D");

// LinearElastic 2D
static Register<mpm::Material<2>, mpm::LinearElastic<2>, unsigned, const Json&>
    linear_elastic_2d("LinearElastic2D");

// LinearElastic 3D
static Register<mpm::Material<3>, mpm::LinearElastic<3>, unsigned, const Json&>
    linear_elastic_3d("LinearElastic3D");

// MohrCoulomb 2D
static Register<mpm::Material<2>, mpm::MohrCoulomb<2>, unsigned, const Json&>
    mohr_coulomb_2d("MohrCoulomb2D");

// MohrCoulomb 3D
static Register<mpm::Material<3>, mpm::MohrCoulomb<3>, unsigned, const Json&>
    mohr_coulomb_3d("MohrCoulomb3D");

// MohrCoulomb 2D
static Register<mpm::Material<2>, mpm::MohrCoulomb_FrozenSoil<2>, unsigned, const Json&>
    frozen_mohr_coulomb_2d("FrozenMohrCoulomb2D");

// MohrCoulomb 3D
static Register<mpm::Material<3>, mpm::MohrCoulomb_FrozenSoil<3>, unsigned, const Json&>
    frozen_mohr_coulomb_3d("FrozenMohrCoulomb3D");

// Newtonian 2D
static Register<mpm::Material<2>, mpm::Newtonian<2>, unsigned, const Json&>
    newtonian_2d("Newtonian2D");

// Newtonian 3D
static Register<mpm::Material<3>, mpm::Newtonian<3>, unsigned, const Json&>
    newtonian_3d("Newtonian3D");

// Norsand 2D
static Register<mpm::Material<2>, mpm::NorSand<2>, unsigned, const Json&>
    nor_sand_2d("NorSand2D");

// Norsand 3D
static Register<mpm::Material<3>, mpm::NorSand<3>, unsigned, const Json&>
    nor_sand_3d("NorSand3D");

// SANISAND 2D
static Register<mpm::Material<2>, mpm::Sanisand<2>, unsigned, const Json&>
    sanisand_2d("SANISAND2D");

// Rigid 2D
static Register<mpm::Material<2>, mpm::Rigid<2>, unsigned, const Json&>
    rigid_2d("Rigid2D");

// Rigid 3D
static Register<mpm::Material<3>, mpm::Rigid<3>, unsigned, const Json&>
    rigid_3d("Rigid3D");
